!! Simulation_sourceTerms.F90
!!
!! Optional operator-split source-term hook called by Driver after each
!! hydrodynamic step. 


subroutine Simulation_sourceTerms(blockID, dt)

   use Simulation_data
   use Grid_interface,   only : Grid_getBlkIndexLimits,  &
                                Grid_getCellCoords,       &
                                Grid_getPointData,        &
                                Grid_putPointData

   implicit none

#include "constants.h"
#include "Flash.h"


   integer, intent(in) :: blockID
   real,    intent(in) :: dt         ! current hydrodynamic timestep

   integer, dimension(LOW:HIGH, MDIM) :: blkLimits, blkLimitsGC
   real, allocatable :: xCoord(:), yCoord(:), zCoord(:)
   integer :: nxGC, nyGC, nzGC
   integer :: i, j, k
   integer, dimension(MDIM) :: axis

   real :: xcell, ycell, zcell
   real :: dist, r_perp

   !! Current cell state (read back from UNK before possible reset)
   real :: dens_cell, tion_cell, tele_cell
   real :: velx_cell, vely_cell, velz_cell

   !! Diagnostic quantities
   real :: cs_ion, vtot, tratio_cell

   !! Get block index limits and coordinates
   call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

   nxGC = blkLimitsGC(HIGH, IAXIS) - blkLimitsGC(LOW, IAXIS) + 1
   nyGC = blkLimitsGC(HIGH, JAXIS) - blkLimitsGC(LOW, JAXIS) + 1
   nzGC = blkLimitsGC(HIGH, KAXIS) - blkLimitsGC(LOW, KAXIS) + 1

   allocate(xCoord(nxGC)); allocate(yCoord(nyGC)); allocate(zCoord(nzGC))

   call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCoord, nxGC)
   call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yCoord, nyGC)
   call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zCoord, nzGC)

   !! Loop over INTERIOR cells only

   do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      zcell = zCoord(k - blkLimitsGC(LOW, KAXIS) + 1)

      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
         ycell = yCoord(j - blkLimitsGC(LOW, JAXIS) + 1)

         do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
            xcell = xCoord(i - blkLimitsGC(LOW, IAXIS) + 1)

            dist   = sqrt(xcell**2 + ycell**2 + zcell**2)
            r_perp = sqrt(xcell**2 + ycell**2)

            !! Sphere replenishment
            !! Reset cells inside R_sph to injection state every timestep.

            if (dist <= sim_Rsph) then

               axis(IAXIS) = i
               axis(JAXIS) = j
               axis(KAXIS) = k

               !! Overwrite with injection state
               call Grid_putPointData(blockID, CENTER, DENS_VAR, INTERIOR, &
                                      axis, sim_rhoSph)
               call Grid_putPointData(blockID, CENTER, VELX_VAR, INTERIOR, &
                                      axis, 0.0e0)
               call Grid_putPointData(blockID, CENTER, VELY_VAR, INTERIOR, &
                                      axis, 0.0e0)
               call Grid_putPointData(blockID, CENTER, VELZ_VAR, INTERIOR, &
                                      axis, 0.0e0)
               call Grid_putPointData(blockID, CENTER, TION_VAR, INTERIOR, &
                                      axis, sim_TiSph)
               call Grid_putPointData(blockID, CENTER, TELE_VAR, INTERIOR, &
                                      axis, sim_TeSph)

            else

               !! Update diagnostic variables for non-sphere cells
               !! Read current UNK
               axis(IAXIS) = i
               axis(JAXIS) = j
               axis(KAXIS) = k

               call Grid_getPointData(blockID, CENTER, TION_VAR, INTERIOR, &
                                      axis, tion_cell)
               call Grid_getPointData(blockID, CENTER, TELE_VAR, INTERIOR, &
                                      axis, tele_cell)
               call Grid_getPointData(blockID, CENTER, VELX_VAR, INTERIOR, &
                                      axis, velx_cell)
               call Grid_getPointData(blockID, CENTER, VELY_VAR, INTERIOR, &
                                      axis, vely_cell)
               call Grid_getPointData(blockID, CENTER, VELZ_VAR, INTERIOR, &
                                      axis, velz_cell)

            end if

         end do  ! i
      end do  ! j
   end do  ! k

   deallocate(xCoord, yCoord, zCoord)

   return
end subroutine Simulation_sourceTerms
