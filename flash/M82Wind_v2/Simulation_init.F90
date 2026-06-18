!! Simulation_init.F90

subroutine Simulation_init()

   use Simulation_data
   use Driver_interface,           only : Driver_getMype
   use RuntimeParameters_interface, only : RuntimeParameters_get
   use PhysicalConstants_interface, only : PhysicalConstants_get

   implicit none

#include "constants.h"
#include "Flash.h"


   real :: theta_rad   ! sim_coneHalfAngle converted to radians

   !! MPI rank
   call Driver_getMype(MESH_COMM, sim_meshMe)

   call RuntimeParameters_get("sim_Lbox",          sim_Lbox)
   call RuntimeParameters_get("sim_Rsph",          sim_Rsph)
   call RuntimeParameters_get("sim_coneHalfAngle", sim_coneHalfAngle)

   !! initial-condition parameters
   call RuntimeParameters_get("sim_nAmb",  sim_nAmb)
   call RuntimeParameters_get("sim_Tamb",  sim_Tamb)
   call RuntimeParameters_get("sim_nSph",  sim_nSph)
   call RuntimeParameters_get("sim_TiSph", sim_TiSph)
   call RuntimeParameters_get("sim_TeSph", sim_TeSph)

   !! physical-constant parameters
   call RuntimeParameters_get("sim_kB",         sim_kB)
   call RuntimeParameters_get("sim_mp",         sim_mp)
   call RuntimeParameters_get("sim_mu",         sim_mu)
   call RuntimeParameters_get("sim_gamma",      sim_gamma)
   call RuntimeParameters_get("sim_coulombLog", sim_coulombLog)

   !! AMR refinement thresholds
   call RuntimeParameters_get("sim_refDensCutoff", sim_refDensCutoff)
   call RuntimeParameters_get("sim_refTionCutoff", sim_refTionCutoff)

   !! cone tangent
   theta_rad          = sim_coneHalfAngle * sim_deg2rad
   sim_tanHalfAngle   = tan(theta_rad)

   !! Ambient mass density
   sim_rhoAmb = sim_mu * sim_mp * sim_nAmb

   !! ambient specific internal energy (each fluid)
   sim_eAmbIon = sim_kB * sim_Tamb / ((sim_gamma - 1.0e0) * sim_mu * sim_mp)
   sim_eAmbEle = sim_eAmbIon   ! T_ion = T_ele = T_amb initially

   !! Sphere mass density
   sim_rhoSph = sim_mu * sim_mp * sim_nSph

   !! Sphere specific internal energy per fluid
   sim_eSphIon = sim_kB * sim_TiSph / ((sim_gamma - 1.0e0) * sim_mu * sim_mp)
   sim_eSphEle = sim_kB * sim_TeSph / ((sim_gamma - 1.0e0) * sim_mu * sim_mp)

   !! Injection source paremeters
   call RuntimeParameters_get("sim_Mdot_T", sim_Mdot_T)
   call RuntimeParameters_get("sim_Edot_T", sim_Edot_T)
   call RuntimeParameters_get("sim_Delta",  sim_Delta)
   
   ! Volumetric injection rates (CC85 normalisation)
   ! Q0 = (3-Delta)/3 * Mdot_T / ( (4/3) * pi * R^3 )
   ! q0 = (3-Delta)/3 * Edot_T / ( (4/3) * pi * R^3 )
   sim_volMassRate = (3.0 - sim_Delta) / 3.0 * sim_Mdot_T / ((4.0/3.0) * 3.141592653589793 * sim_Rsph**3)
   sim_volEnergyRate = (3.0 - sim_Delta) / 3.0 * sim_Edot_T / ((4.0/3.0) * 3.141592653589793 * sim_Rsph**3)

   !! diagnostics (printed by rank 0 only)
   if (sim_meshMe == MASTER_PE) then
      write(*,*)
      write(*,'(A)') "========================================================"
      write(*,'(A)') "  M82 Biconical Wind -- Simulation_init summary"
      write(*,'(A)') "  Ref: Boettcher et al. 2026, Nature 651, 909"
      write(*,'(A)') "========================================================"
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
      write(*,'(A)') "  --- Injection source ---"
      write(*,'(A,ES12.4,A)') "  Mdot_T            = ", sim_Mdot_T, " g/s"
      write(*,'(A,ES12.4,A)') "  Edot_T            = ", sim_Edot_T, " erg/s"
      write(*,'(A,F8.2)')     "  Delta             = ", sim_Delta
      write(*,'(A,ES12.4,A)') "  Q0 (mass/vol/s)   = ", sim_volMassRate, " g/cm^3/s"
      write(*,'(A,ES12.4,A)') "  q0 (energy/vol/s) = ", sim_volEnergyRate, " erg/cm^3/s"
      write(*,'(A)') "========================================================"
      write(*,*)
   end if

end subroutine Simulation_init
