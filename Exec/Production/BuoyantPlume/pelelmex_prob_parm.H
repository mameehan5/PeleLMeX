#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
   // Shared
   int probType = 0;
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real Vin      = -1.0;
   amrex::Real Tin      = 298.0;
   amrex::Real Vcoflow  = 0.0;
   amrex::Real Tcoflow  = 298.0;
   amrex::Real Lscale = 0.0;
   amrex::Real smoothing_factor = 0.0;
   int nside = 100;
   int njet = 1;
   amrex::Real spacing = 0.0;
   amrex::Real aspect_ratio = 0.0;
   amrex::Real width = 0.0;
   amrex::Real length = 0.0;
   amrex::Real equivalence_ratio = 0.0;
   amrex::Real ce = 0.0;

   amrex::Real wall_Lscale = 0.0;
   amrex::Real Twall = Tcoflow;
   int wall_nside = 100;
   amrex::Real wall_aspect_ratio = aspect_ratio;
   amrex::Real wall_smoothing_factor = smoothing_factor;

   amrex::Real heating_max_temp = Tin;
   amrex::Real heating_time_low = -1.0;
   amrex::Real heating_time_high = -1.0;
   amrex::Real heating_smoothing_factor = 0.0;
};
#endif
