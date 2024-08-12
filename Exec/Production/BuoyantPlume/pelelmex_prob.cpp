#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>
#include <PMFData.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");

   std::string type;
   pp.query("type", type);

   pp.query("P_mean", prob_parm->P_mean);
   pp.query("Vin", PeleLM::prob_parm->Vin);
   pp.query("Tin", PeleLM::prob_parm->Tin);
   pp.query("Vcoflow", PeleLM::prob_parm->Vcoflow);
   pp.query("Tcoflow", PeleLM::prob_parm->Tcoflow);
   pp.query("Lscale", PeleLM::prob_parm->Lscale);
   pp.query("smoothing_factor", PeleLM::prob_parm->smoothing_factor);
   pp.query("nside", PeleLM::prob_parm->nside);
   pp.query("njet", PeleLM::prob_parm->njet);
   pp.query("spacing", PeleLM::prob_parm->spacing);
   pp.query("aspect_ratio", PeleLM::prob_parm->aspect_ratio);
   pp.query("width", PeleLM::prob_parm->width);
   pp.query("length", PeleLM::prob_parm->length);
   pp.query("equivalence_ratio", PeleLM::prob_parm->equivalence_ratio);
   pp.query("ce", PeleLM::prob_parm->ce);

   pp.query("wall_Lscale", PeleLM::prob_parm->wall_Lscale);
   pp.query("Twall", PeleLM::prob_parm->Twall);
   pp.query("wall_nside", PeleLM::prob_parm->wall_nside);
   pp.query("wall_aspect_ratio", PeleLM::prob_parm->wall_aspect_ratio);
   pp.query("wall_smoothing_factor", PeleLM::prob_parm->wall_smoothing_factor);

   pp.query("heating_max_temp",PeleLM::prob_parm->heating_max_temp);
   pp.query("heating_time_low",PeleLM::prob_parm->heating_time_low);
   pp.query("heating_time_high",PeleLM::prob_parm->heating_time_high);
   pp.query("heating_smoothing_factor",PeleLM::prob_parm->heating_smoothing_factor);


   if ( type == "NonReacting" ) {
      PeleLM::prob_parm->probType = 0;
   } else if ( type == "Reacting" ) {
      PeleLM::prob_parm->probType = 1;
   } 

   //PeleLM::pmf_data.initialize(); 
}
