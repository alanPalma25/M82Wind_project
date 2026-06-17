!! Simulation_initBlock.F90
!!
!! Initialises a single AMR block for the M82 Biconical Wind simulation.

subroutine Simulation_initBlock(blockID)

   use Simulation_data
   use Grid_interface,       only : Grid_getBlkIndexLimits,  &
                                    Grid_getCellCoords,       &
                                    Grid_putPointData
   use Eos_interface,        only : Eos_wrapped

   implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

   integer, intent(in) :: blockID


   !! Block index limits in the UNK array (interior + guard cells)
   integer, dimension(LOW:HIGH, MDIM) :: blkLimits 
   integer, dimension(LOW:HIGH, MDIM) :: blkLimitsGC
   integer, dimension(MDIM) :: axis

   !! Cell-centre coordinate arrays for each axis
   real, allocatable :: xCoord(:), yCoord(:), zCoord(:)

   !! Number of cells along each axis including guard cells
   integer :: nxGC, nyGC, nzGC

   !! Cell-loop indices
   integer :: i, j, k

   !! Cell centroid coordinates
   real :: xcell, ycell, zcell

   !! Distance from box centre (origin)
   real :: dist

   !! Perpendicular distance from z-axis: r_perp = sqrt(x^2 + y^2)
   real :: r_perp

   !! Flags
   logical :: inSphere    ! cell is inside the starburst sphere
   logical :: inConeFunnel! cell is inside the biconical funnel (not sphere)

   !! Scratchpad for a single cell's fluid state, passed to Eos_wrapped
   real :: tion_cell, tele_cell, dens_cell
   real :: velx_cell, vely_cell, velz_cell

   !! Diagnostic Mach number (for MACH_VAR diagnostic variable)
   real :: cs_ion, vtot, tratio_cell

   !! block's index ranges

   call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

   nxGC = blkLimitsGC(HIGH, IAXIS) - blkLimitsGC(LOW, IAXIS) + 1
   nyGC = blkLimitsGC(HIGH, JAXIS) - blkLimitsGC(LOW, JAXIS) + 1
   nzGC = blkLimitsGC(HIGH, KAXIS) - blkLimitsGC(LOW, KAXIS) + 1


   !!cell-centre coordinates (including guard cells)

   allocate(xCoord(nxGC))
   allocate(yCoord(nyGC))
   allocate(zCoord(nzGC))

   call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCoord, nxGC)
   call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yCoord, nyGC)
   call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zCoord, nzGC)

   !! Loop over all cells INCLUDING guard cells
   do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
      zcell = zCoord(k - blkLimitsGC(LOW, KAXIS) + 1)

      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
         ycell = yCoord(j - blkLimitsGC(LOW, JAXIS) + 1)

         do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
            xcell = xCoord(i - blkLimitsGC(LOW, IAXIS) + 1)

            !! Spherical distance from origin
            dist = sqrt(xcell**2 + ycell**2 + zcell**2)

            !! Cylindrical radius (perpendicular distance from z-axis)
            r_perp = sqrt(xcell**2 + ycell**2)

            !! nside the starburst sphere
            inSphere = (dist <= sim_Rsph)

            !! Region B: inside the biconical funnel but outside the sphere
            inConeFunnel = (.not. inSphere) .and. &
                           (r_perp <= sim_tanHalfAngle * abs(zcell))

            if (inSphere) then
               !! sphere

               velx_cell = 0.0e0
               vely_cell = 0.0e0
               velz_cell = 0.0e0
               tion_cell = sim_TiSph
               tele_cell = sim_TeSph

            else
               !! Ambient medium (and cone-wall cells).

               dens_cell = sim_rhoAmb
               velx_cell = 0.0e0
               vely_cell = 0.0e0
               velz_cell = 0.0e0
               tion_cell = sim_Tamb
               tele_cell = sim_Tamb

            end if
            
            axis(IAXIS) = i
            axis(JAXIS) = j
            axis(KAXIS) = k


            call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, dens_cell)
            call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, velx_cell)
            call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, vely_cell)
            call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, velz_cell)
            call Grid_putPointData(blockID, CENTER, TION_VAR, EXTERIOR, axis, tion_cell)
            call Grid_putPointData(blockID, CENTER, TELE_VAR, EXTERIOR, axis, tele_cell)


         end do  ! i
      end do  ! j
   end do  ! k

   call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blockID)

   deallocate(xCoord)
   deallocate(yCoord)
   deallocate(zCoord)

   return
end subroutine Simulation_initBlock
