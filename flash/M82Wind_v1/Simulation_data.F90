!! Simulation_data.F90
!!
!! Module holding all runtime parameters 

module Simulation_data

   implicit none

#include "constants.h"
#include "Flash.h"


   real, save :: sim_Lbox
   real, save :: sim_Rsph

   !! Half-opening angle of the biconical funnel [degrees].
   real, save :: sim_coneHalfAngle

   !! Precomputed tan(theta_half) [dimensionless]it.
   real, save :: sim_tanHalfAngle

   !! Ambient number density [cm^-3]
   real, save :: sim_nAmb

   !! Ambient temperature for BOTH ion and electron fluids [K].
   real, save :: sim_Tamb

   !! Derived: ambient mass density [g/cm^3] = mu * m_p * n_amb
   real, save :: sim_rhoAmb

   !! Derived: ambient specific internal energy for EACH fluid [erg/g].
   !! Since T_ion = T_ele = T_amb, e_ion = e_ele = (kB*T_amb)/(gamma-1)/mu/m_p
   real, save :: sim_eAmbIon
   real, save :: sim_eAmbEle

   !! Number density inside the sphere [cm^-3]
   real, save :: sim_nSph

   !! Ion (proton) temperature inside the sphere [K].
   real, save :: sim_TiSph

   !! Electron temperature inside the sphere [K].
   real, save :: sim_TeSph

   !! Derived: sphere mass density [g/cm^3]
   real, save :: sim_rhoSph

   !! Derived: sphere specific internal energy for ions [erg/g]
   real, save :: sim_eSphIon

   !! Derived: sphere specific internal energy for electrons [erg/g]
   real, save :: sim_eSphEle--

   !! Boltzmann constant [erg/K]
   real, save :: sim_kB

   !! Proton mass [g]
   real, save :: sim_mp

   !! Mean molecular weight (mu = 0.5 for fully ionised pure hydrogen)
   real, save :: sim_mu

   !! Adiabatic index (gamma = 5/3 for monatomic ideal gas)
   real, save :: sim_gamma

   !! Coulomb logarithm for ion-electron energy exchange.
   real, save :: sim_coulombLog

   !! MPI rank of the master process (set in Simulation_init)
   integer, save :: sim_meshMe

   !! Conversion: degrees -> radians (set once in Simulation_init)
   real, parameter :: sim_deg2rad = 3.14159265358979323846e0 / 180.0e0

end module Simulation_data
