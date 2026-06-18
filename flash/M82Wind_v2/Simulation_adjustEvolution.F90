subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

  use Simulation_data, only : sim_Rsph, sim_volMassRate, sim_volEnergyRate, sim_Delta
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
  real :: xcell, ycell, zcell, dist, r_scaled
  real :: rho_old, vx_old, vy_old, vz_old, eion_old, eele_old, ener_old
  real :: rho_new, eion_new, eele_new, ener_new
  real :: mass_rate, energy_rate
  integer :: myPE

  call Driver_getMyPE(MESH_COMM, myPE)
  if (myPE == 0 .and. (nstep == 1 .or. mod(nstep, 10) == 0)) then
     write(*,*) 'Simulation_adjustEvolution: step=', nstep, ' time=', stime
  end if

  do lb = 1, blkcnt
     blockID = blklst(lb)

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

                 ! 1. Read current state
                 call Grid_getPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, rho_old)
                 call Grid_getPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, vx_old)
                 call Grid_getPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, vy_old)
                 call Grid_getPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, vz_old)
                 call Grid_getPointData(blockID, CENTER, EION_VAR, EXTERIOR, axis, eion_old)
                 call Grid_getPointData(blockID, CENTER, EELE_VAR, EXTERIOR, axis, eele_old)
                 call Grid_getPointData(blockID, CENTER, ENER_VAR, EXTERIOR, axis, ener_old)

                 ! 2. Local injection rates (avoid singularity at r=0)
                 r_scaled = max(dist, 1.0e-10 * sim_Rsph)
                 mass_rate = sim_volMassRate * (sim_Rsph / r_scaled)**sim_Delta
                 energy_rate = sim_volEnergyRate * (sim_Rsph / r_scaled)**sim_Delta

                 ! 3. Add mass to density
                 rho_new = rho_old + mass_rate * dt

                 ! 4. Add energy to ion fluid only (SNe heat ions)
                 !    new specific internal energy = (old energy density + added) / new density
                 eion_new = (rho_old * eion_old + energy_rate * dt) / rho_new

                 ! 5. Electron internal energy unchanged (will be updated by Coulomb coupling)
                 eele_new = eele_old

                 ! 6. Velocities unchanged (zero momentum source)
                 !    (vx, vy, vz remain as read)

                 ! 7. Update total energy density: rho*(eion+eele+0.5*v^2)
                 ener_new = rho_new * (eion_new + eele_new + 0.5*(vx_old**2 + vy_old**2 + vz_old**2))

                 ! 8. Write updated conserved variables
                 call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, rho_new)
                 call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, vx_old)
                 call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, vy_old)
                 call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, vz_old)
                 call Grid_putPointData(blockID, CENTER, EION_VAR, EXTERIOR, axis, eion_new)
                 call Grid_putPointData(blockID, CENTER, EELE_VAR, EXTERIOR, axis, eele_new)
                 call Grid_putPointData(blockID, CENTER, ENER_VAR, EXTERIOR, axis, ener_new)

                 ! Do NOT set TION_VAR or TELE_VAR here; Eos will compute them.
              end if
           end do
        end do
     end do

     ! Call EOS to compute pressure and temperatures from density and internal energies
     call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blockID)

     deallocate(xCoord, yCoord, zCoord)

  end do  ! end block loop

end subroutine Simulation_adjustEvolution
