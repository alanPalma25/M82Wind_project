!! Simulation_init.F90
!!
!! Called once at startup by FLASH's Driver unit before the main time loop.
subroutine Simulation_init()

   use Simulation_data
   use Driver_interface,           only : Driver_getMype
   use RuntimeParameters_interface, only : RuntimeParameters_get
   use PhysicalConstants_interface, only : PhysicalConstants_get

   implicit none

#include "constants.h"
#include "Flash.h"


   real :: theta_rad   ! converted to radians


   !! MPI rank
   call Driver_getMype(MESH_COMM, sim_meshMe)

   !! Read geometry / domain parameters
   call RuntimeParameters_get("sim_Lbox",          sim_Lbox)
   call RuntimeParameters_get("sim_Rsph",          sim_Rsph)
   call RuntimeParameters_get("sim_coneHalfAngle", sim_coneHalfAngle)

   !! Read initial-condition parameters
   call RuntimeParameters_get("sim_nAmb",  sim_nAmb)
   call RuntimeParameters_get("sim_Tamb",  sim_Tamb)
   call RuntimeParameters_get("sim_nSph",  sim_nSph)
   call RuntimeParameters_get("sim_TiSph", sim_TiSph)
   call RuntimeParameters_get("sim_TeSph", sim_TeSph)

   !! Read physical-constant parameters
   call RuntimeParameters_get("sim_kB",         sim_kB)
   call RuntimeParameters_get("sim_mp",         sim_mp)
   call RuntimeParameters_get("sim_mu",         sim_mu)
   call RuntimeParameters_get("sim_gamma",      sim_gamma)
   call RuntimeParameters_get("sim_coulombLog", sim_coulombLog)

   !! Read AMR refinement threshold
   call RuntimeParameters_get("sim_refDensCutoff", sim_refDensCutoff)
   call RuntimeParameters_get("sim_refTionCutoff", sim_refTionCutoff)

   !! Derived: cone tangent
   !! r_perp = tan(theta_half) * |z|  defines the cone surface.

   theta_rad          = sim_coneHalfAngle * sim_deg2rad
   sim_tanHalfAngle   = tan(theta_rad)

   !! Derived: ambient mass density
   sim_rhoAmb = sim_mu * sim_mp * sim_nAmb

   !! 8.  Derived: ambient specific internal energy (each fluid)
   !!     e = (kB * T) / ((gamma-1) * mu * m_p)
   !!     Both fluids start in thermal equilibrium for the ambient medium.
   sim_eAmbIon = sim_kB * sim_Tamb / ((sim_gamma - 1.0e0) * sim_mu * sim_mp)
   sim_eAmbEle = sim_eAmbIon   ! T_ion = T_ele = T_amb initially

   !! Derived: sphere mass density
   sim_rhoSph = sim_mu * sim_mp * sim_nSph

   !! Derived: sphere specific internal energy per fluid
   !!       t_eq ~ 252 yr * T_ion^(3/2) / (n_ele * lnLambda)
   !!            ~ 252 * (1.16e8)^1.5 / (100 * 35)  ~ 3 kyr
   !!     Reference: Spitzer (1962); Strickland & Heckman (2009).

   sim_eSphIon = sim_kB * sim_TiSph / ((sim_gamma - 1.0e0) * sim_mu * sim_mp)
   sim_eSphEle = sim_kB * sim_TeSph / ((sim_gamma - 1.0e0) * sim_mu * sim_mp)


   !! Sanity diagnostics (printed by rank 0 only)

   if (sim_meshMe == MASTER_PE) then
      write(*,*)
      write(*,'(A)') "========================================================"
      write(*,'(A)') "  M82 Biconical Wind -- Simulation_init summary"
      write(*,'(A,ES12.4,A)') "  Box half-length      sim_Lbox  = ", sim_Lbox,   " cm"
      write(*,'(A,ES12.4,A)') "  Sphere radius        sim_Rsph  = ", sim_Rsph,   " cm"
      write(*,'(A,F8.2,  A)') "  Cone half-angle               = ", sim_coneHalfAngle, " deg"
      write(*,'(A,ES12.4  )') "  tan(half-angle)               = ", sim_tanHalfAngle
      write(*,'(A)') "  --- Ambient ---"
      write(*,'(A,ES12.4,A)') "  n_amb                = ", sim_nAmb,    " cm^-3"
      write(*,'(A,ES12.4,A)') "  T_amb                = ", sim_Tamb,    " K"
      write(*,'(A,ES12.4,A)') "  rho_amb              = ", sim_rhoAmb,  " g/cm^3"
      write(*,'(A,ES12.4,A)') "  e_amb (ion=ele)      = ", sim_eAmbIon, " erg/g"
      write(*,'(A)') "  --- Sphere (starburst source) ---"
      write(*,'(A,ES12.4,A)') "  n_sph                = ", sim_nSph,    " cm^-3"
      write(*,'(A,ES12.4,A)') "  T_ion_sph            = ", sim_TiSph,   " K  (10 keV)"
      write(*,'(A,ES12.4,A)') "  T_ele_sph            = ", sim_TeSph,   " K  (0.1 keV)"
      write(*,'(A,ES12.4  )') "  T_ion/T_ele (init)   = ", sim_TiSph/sim_TeSph
      write(*,'(A,ES12.4,A)') "  rho_sph              = ", sim_rhoSph,  " g/cm^3"
      write(*,'(A,ES12.4,A)') "  e_ion_sph            = ", sim_eSphIon, " erg/g"
      write(*,'(A,ES12.4,A)') "  e_ele_sph            = ", sim_eSphEle, " erg/g"
      write(*,'(A,F8.1    )') "  Coulomb logarithm    = ", sim_coulombLog
      write(*,'(A)') "========================================================"
      write(*,*)
   end if

end subroutine Simulation_init
