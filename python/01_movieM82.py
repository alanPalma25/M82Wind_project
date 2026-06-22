# Import third-party libraries
from csv import writer

import yt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import constants
import os
import scienceplots
import pandas as pd
import multiprocessing as mp
import glob
from PIL import Image 
import imageio.v2 as imageio

def output_dir_check(output_dir = None):
    
    # Create the output directory on the current direcctory if the output_dir is None
    if output_dir is None:
        output_dir = os.getcwd()

    # Create an output directory for the plots if it doesn't exist
    if os.path.isdir(output_dir + "/" + "output"):
        print("Directory already exists.")
    else:
        print("Directory has been created.")
        os.mkdir(output_dir + "/" + "output")
    
    return output_dir

class CreateMovie:
    """
    A class to create a movie from a series of FLASH HDF5 outputs.
    """

    def __init__(self, path=None, n_cpus=1, output_name="m82_movie.mp4", fps=10, step_imgs=None):

        # Deactivate the logging from yt in the workers
        yt.set_log_level(40)

        # Build the glob pattern
        if path is None:
            pattern = os.getcwd() + "/M82Wind_2D_hdf5_plt_cnt_????"
        else:
            pattern = path

        # Store the SORTED FILE LIST names
        self.file_list = sorted(glob.glob(pattern))
        if not self.file_list:
            raise FileNotFoundError(f"No files matched pattern: {pattern}")

        self.output_folder = output_dir_check() + "/output"
        os.makedirs(self.output_folder, exist_ok=True)

        self.n_cpus    = n_cpus
        self.step_imgs = step_imgs if step_imgs is not None else 1
        self.fps       = fps
        self.images_name  = "dens_M82.*.png"
        self.output_name  = output_name

        # plt.style.use(['science', 'notebook', 'no-latex'])

        # Compute global colour limits from the FIRST file only (with YT)
        ds0   = yt.load(self.file_list[0])
        dens0 = ds0.all_data()[("flash", "dens")]
        self.vmin = np.min(dens0).to("g/cm**3")
        self.vmax = np.max(dens0).to("g/cm**3")

    def plot_density(self, j):
        """
        Plot density slice for frame index j (into self.file_list).
        Each worker loads its own yt dataset — no shared state needed.
        """
        filepath = self.file_list[j]          # plain string
        ds       = yt.load(filepath)           # load inside the worker

        # Extract the time from the dataset
        time = ds.current_time.to("Myr")

        # Create the slice plot

        p_slice_r = yt.SlicePlot(ds, "z", ("flash", "dens"))
        p_slice_r.set_cmap(("flash", "dens"), "turbo")
        p_slice_r.set_unit(("flash", "dens"), "g/cm**3")

        p_slice_r.set_zlim(("flash", "dens"), self.vmin, self.vmax)
        p_slice_r.annotate_title("t = {:.3f} Myr".format(time.value))

        # Save the plot
        out = f"{self.output_folder}/dens_M82.{j:03d}.png"
        print(f"Saving frame {j} to {out} …")

        p_slice_r.save(out)

    def create_movie_mpi(self):
        """
        Same as create_movie, but parallelizes frame generation using MPI
        instead of multiprocessing.Pool. Each rank renders its own subset
        of frames; rank 0 assembles the final MP4 once all ranks are done.
        """
    
        from mpi4py import MPI

        world_comm = MPI.COMM_WORLD
        world_size = world_comm.Get_size()
        my_rank    = world_comm.Get_rank()

        # Build frame index list
        imgs_i = list(range(0, len(self.file_list), self.step_imgs))
        if imgs_i[-1] != len(self.file_list) - 1:
            imgs_i.append(len(self.file_list) - 1)
        N_frames = len(imgs_i)

        if my_rank == 0:
            print(f"Generating {N_frames} frames with {world_size} MPI ranks …")

        # Distribute frame indices across ranks
        # as the temperature example)
        workloads = [N_frames // world_size for _ in range(world_size)]
        for i in range(N_frames % world_size):
            workloads[i] += 1

        my_start = sum(workloads[:my_rank])
        my_end   = my_start + workloads[my_rank]
        my_imgs_i = imgs_i[my_start:my_end]

        print(f"Rank {my_rank} assigned frames: {my_imgs_i}")

        # Each rank renders its own frames
        for j in my_imgs_i:
            self.plot_density(j)

        # Completion signal to root
        if my_rank != 0:
            world_comm.send(True, dest=0, tag=77)
            print(f"Rank {my_rank} sent completion signal")
        else:
            for i in range(1, world_size):
                world_comm.recv(source=i, tag=77)
                print(f"Rank 0 received completion signal from rank {i}")

        # Barrier to be safe before any rank exits and root assembles video
        world_comm.Barrier()

        # Only rank 0 assembles the MP4
        if my_rank == 0:
            images_input = f"{self.output_folder}/{self.images_name}"
            output_video = f"{self.output_folder}/{self.output_name}"
            filenames    = sorted(glob.glob(images_input))
            with imageio.get_writer(output_video, fps=self.fps, codec='libx264') as writer:
                for fn in filenames:
                    writer.append_data(imageio.imread(fn))
            print("MP4 saved to:", output_video)