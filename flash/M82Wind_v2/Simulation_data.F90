!! Simulation_data.F90
!!
!! Module holding all runtime parameters and derived constants for the
!! "M82 Biconical Wind with Two-Temperature Plasma" FLASH simulation.

module Simulation_data

   implicit none

#include "constants.h"
#include "Flash.h"

   !! Domain and geometry parameters

   !! Half-side of the cubic box [cm].
   real, save :: sim_Lbox

   !! Radius of the starburst source
   real, save :: sim_Rsph

   !! Half-opening angle of the biconical funnel [degrees].
   real, save :: sim_coneHalfAngle

   !! Precomputed tan(theta_half)
   real, save :: sim_tanHalfAngle

   !! Ambient (background / inter-cone) initial conditions

   !! Ambient number density [cm^-3]
   real, save :: sim_nAmb

   !! Ambient temperature for BOTH ion and electron fluids [K].
   real, save :: sim_Tamb

   !! Ambient mass density [g/cm^3] = mu * m_p * n_amb
   real, save :: sim_rhoAmb

   !! Ambient specific internal energy for EACH fluid [erg/g].
   real, save :: sim_eAmbIon
   real, save :: sim_eAmbEle

   !! 3. Sphere (starburst source) initial conditions

   !! Number density inside the sphere
   real, save :: sim_nSph

   !! Ion (proton) temperature inside the sphere
   real, save :: sim_TiSph

   !! Electron temperature inside the sphere
   real, save :: sim_TeSph

   !! Derived: sphere mass density
   real, save :: sim_rhoSph

   !! Derived: sphere specific internal energy for ions
   real, save :: sim_eSphIon

   !! Derived: sphere specific internal energy for electrons
   real, save :: sim_eSphEle

   !! Physical constants and gas parameters

   !! Boltzmann constant
   real, save :: sim_kB

   !! Proton mass
   real, save :: sim_mp

   !! Mean molecular weight
   real, save :: sim_mu

   !! Adiabatic index (gamma = 5/3 for monatomic ideal gas)
   real, save :: sim_gamma

   !! Coulomb logarithm for ion-electron energy exchange.
   real, save :: sim_coulombLog

   !! MPI rank of the master process (set in Simulation_init)
   integer, save :: sim_meshMe

   !! Conversion
   real, parameter :: sim_deg2rad = 3.14159265358979323846e0 / 180.0e0

   !! Injection source
   
   !! Total mass injection rate
   real, save :: sim_Mdot_T

   !! Total energy injection rate [erg/s]
   real, save :: sim_Edot_T

   !! Power-law index (0 = uniform CC85)
   real, save :: sim_Delta

   !! Derived volumetric injection rates [g/cm^3/s] and [erg/cm^3/s]
   real, save :: sim_volMassRate
   real, save :: sim_volEnergyRate

end module Simulation_data
