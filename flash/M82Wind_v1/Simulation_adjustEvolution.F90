!!****f* source/Simulation/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  Here we enforce a steady injection reservoir inside the starburst sphere.
!!
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***

subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

  use Simulation_data, only : sim_Rsph, sim_rhoSph, sim_TiSph, sim_TeSph, &
                              sim_eSphIon, sim_eSphEle
  use Grid_interface,   only : Grid_getBlkIndexLimits, &
                               Grid_getCellCoords,      &
                               Grid_putPointData,       &
                               Grid_getPointData
  use Eos_interface,    only : Eos_wrapped
  use Driver_interface, only : Driver_getMyPE

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real,    intent(in) :: dt
  real,    intent(in) :: stime

  ! Local variables
  integer :: blockID, lb
  integer, dimension(LOW:HIGH, MDIM) :: blkLimits, blkLimitsGC
  real, allocatable :: xCoord(:), yCoord(:), zCoord(:)
  integer :: nxGC, nyGC, nzGC
  integer :: i, j, k
  integer, dimension(MDIM) :: axis
  real :: xcell, ycell, zcell, dist
  integer :: myPE
  
  ! Check if the routine is executing
  call Driver_getMyPE(MESH_COMM, myPE)
  if (myPE == 0 .and. (nstep == 1 .or. mod(nstep, 10) == 0)) then
     write(*,*) 'Simulation_adjustEvolution: step=', nstep, ' time=', stime
  end if


  ! Loop over all blocks
  do lb = 1, blkcnt
     blockID = blklst(lb)

     ! Get block limits
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     nxGC = blkLimitsGC(HIGH, IAXIS) - blkLimitsGC(LOW, IAXIS) + 1
     nyGC = blkLimitsGC(HIGH, JAXIS) - blkLimitsGC(LOW, JAXIS) + 1
     nzGC = blkLimitsGC(HIGH, KAXIS) - blkLimitsGC(LOW, KAXIS) + 1

     allocate(xCoord(nxGC), yCoord(nyGC), zCoord(nzGC))

     call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCoord, nxGC)
     call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yCoord, nyGC)
     call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zCoord, nzGC)

     ! Loop over interior cells only
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        zcell = zCoord(k - blkLimitsGC(LOW, KAXIS) + 1)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           ycell = yCoord(j - blkLimitsGC(LOW, JAXIS) + 1)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              xcell = xCoord(i - blkLimitsGC(LOW, IAXIS) + 1)

              dist = sqrt(xcell**2 + ycell**2 + zcell**2)

              if (dist <= sim_Rsph) then
                 axis(IAXIS) = i
                 axis(JAXIS) = j
                 axis(KAXIS) = k

                 ! Reset to initial conditions
                 call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoSph)
                 call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, 0.0e0)
                 call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, 0.0e0)
                 call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, 0.0e0)
                 call Grid_putPointData(blockID, CENTER, TION_VAR, EXTERIOR, axis, sim_TiSph)
                 call Grid_putPointData(blockID, CENTER, TELE_VAR, EXTERIOR, axis, sim_TeSph)

                 ! Do NOT set internal energies; Eos_wrapped will recompute them.
              end if
           end do
        end do
     end do

     ! Call EOS to recompute internal energies, pressure, total energy
     call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blockID)

     deallocate(xCoord, yCoord, zCoord)

  end do  
