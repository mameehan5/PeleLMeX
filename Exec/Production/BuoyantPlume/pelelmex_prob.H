#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelmex_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real u[3] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};  
    amrex::Real y1, y2;
    amrex::Real splitx = prob_lo[0] + 0.5*Lx;
    amrex::Real splity = prob_lo[1] + 0.5*Ly;
#if ( AMREX_SPACEDIM == 3 )
    amrex::Real splitz = prob_lo[2] + 0.5*Lz;
#endif
    amrex::Real phi = prob_parm.equivalence_ratio;
    massfrac_air[N2_ID] = 0.767;
    massfrac_air[O2_ID] = 0.233;
    if ( phi < 0.01) {
        massfrac_jet[N2_ID] = 0.767;
        massfrac_jet[O2_ID] = 0.233;    
#ifdef CH4_ID
    } else if ( phi > 10.0 ) {
        massfrac_jet[CH4_ID] = 1.0;
    } else {
        amrex::Real mf_ch4 = 0.0584*phi/(1.0 + 0.0584*phi);
        massfrac_jet[CH4_ID] = mf_ch4;
        massfrac_jet[O2_ID] = 0.233*(1.0-mf_ch4);
        massfrac_jet[N2_ID] = 1.0 - massfrac_jet[O2_ID] - massfrac_jet[CH4_ID];
#endif
#ifdef HE_ID
   } else if ( phi > 1.0 ) {
      massfrac_jet[HE_ID] = 1.0;
   } else {
      massfrac_jet[HE_ID] = phi;
      massfrac_jet[O2_ID] = 0.233*(1.0-phi);
      massfrac_jet[N2_ID] = 0.767*(1.0-phi);
#endif
    }

    //amrex::Real Rad = prob_parm.Lscale;
    //amrex::Real rad = std::abs(x - splitx);
    //amrex::Real sf = prob_parm.smoothing_factor;
    //amrex::Real eta = 0.5*(1.0 - tanh(2.0*((rad-Rad)/sf)));
    for (int n = 0; n < NUM_SPECIES; n++){
        //massfrac[n] = eta*massfrac_jet[n] + (1.0-eta)*massfrac_air[n];
        massfrac[n] = massfrac_air[n];
    }
    //state(i,j,k,TEMP) = eta*prob_parm.Tin + (1.0-eta)*prob_parm.Tcoflow;
    state(i,j,k,TEMP) = prob_parm.Tcoflow;
    state(i,j,k,VELX) = 0.0;
    //state(i,j,k,VELY) = eta*prob_parm.Vin + (1.0-eta)*prob_parm.Vcoflow;
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = 0.0;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
    state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
       state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();
   const amrex::Real* prob_hi = geomdata.ProbHi();

   AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                const amrex::Real Lz = prob_hi[2] - prob_lo[2]);
   AMREX_D_TERM(amrex::Real splitx = prob_lo[0] + 0.5 * Lx;,
                amrex::Real splity = prob_lo[1] + 0.5 * Ly;,
                amrex::Real splitz = prob_lo[2] + 0.5 * Lz;)
   amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
   amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
   amrex::Real massfrac[NUM_SPECIES]     = {0.0};
   amrex::Real phi = prob_parm.equivalence_ratio; 

   massfrac_air[N2_ID] = 0.767;
   massfrac_air[O2_ID] = 0.233;
   
   if ( phi < 0.01) {
      massfrac_jet[N2_ID] = 0.767;
      massfrac_jet[O2_ID] = 0.233;    
#ifdef CH4_ID
   } else if ( phi > 10.0 ) {
      massfrac_jet[CH4_ID] = 1.0;
   } else {
      amrex::Real mf_ch4 = 0.0584*phi/(1.0 + 0.0584*phi);
      massfrac_jet[CH4_ID] = mf_ch4;
      massfrac_jet[O2_ID] = 0.233*(1.0-mf_ch4);
      massfrac_jet[N2_ID] = 1.0 - massfrac_jet[O2_ID] - massfrac_jet[CH4_ID];
#endif
#ifdef HE_ID
   } else if ( phi > 1.0 ) {
      massfrac_jet[HE_ID] = 1.0;
   } else {     
      massfrac_jet[HE_ID] = phi;
      massfrac_jet[O2_ID] = 0.233*(1.0-phi);
      massfrac_jet[N2_ID] = 0.767*(1.0-phi);
#endif
   }    

   auto eos = pele::physics::PhysicsType::eos();

   for (int n = 0; n < NVAR; n++){
      s_ext[n] = 0.0;
   }

   if ( ( AMREX_SPACEDIM == 3 && prob_parm.nside != 4 && idir == 2 && sgn == 1 ) || 
        ( AMREX_SPACEDIM == 3 && prob_parm.nside == 4 && idir == 1 && sgn == 1 ) ||
        ( AMREX_SPACEDIM == 2 && idir == 1                         && sgn == 1 ) ) {
      
      amrex::Real Rad = prob_parm.Lscale;
      amrex::Real sf = prob_parm.smoothing_factor;
      amrex::Real ht0 = prob_parm.heating_time_low;
      amrex::Real htf = prob_parm.heating_time_high;
      amrex::Real hsf = prob_parm.heating_smoothing_factor;

      amrex::Real eta = 0.0;
      amrex::Real rad = 0.0;

      amrex::Real etaW = 0.0;
      amrex::Real wRad = prob_parm.wall_Lscale;
      amrex::Real wsf  = prob_parm.wall_smoothing_factor;

#if ( AMREX_SPACEDIM == 2 )
      for (int n = 0; n < prob_parm.njet; n++) {
          rad = std::abs(x[0]-splitx-prob_parm.spacing*(n-((prob_parm.njet-1)/2.0)));
          eta = eta + 0.5*(1.0 - tanh(2.0*((rad-Rad)/sf)));
          if ( wRad > Rad ) {
             etaW = etaW + 0.5*(1.0 - tanh(2.0*((rad-wRad)/wsf)));
          }
      }
#elif ( AMREX_SPACEDIM == 3 )
      if ( prob_parm.nside == 4 ) {
          amrex::Real slength, llength; 
          if ( prob_parm.width > 0.0 and prob_parm.length > 0.0 ) {
              slength = prob_parm.width; 
              llength = prob_parm.length;
          } else {
              slength = prob_parm.Lscale/(prob_parm.aspect_ratio/(2.0 + 2.0*prob_parm.aspect_ratio));
              llength = prob_parm.aspect_ratio*slength;
          }
          amrex::Real BCx1 = -slength/2.0;
          amrex::Real BCx2 =  slength/2.0;
          amrex::Real BCz1 = -llength/2.0;
          amrex::Real BCz2 =  llength/2.0;

          amrex::Real eta_xmin = 0.5*(1.0 - tanh(2.0*(BCx1-x[0])/sf));
          amrex::Real eta_xmax = 0.5*(1.0 - tanh(2.0*(x[0]-BCx2)/sf));
          amrex::Real eta_zmin = 0.5*(1.0 - tanh(2.0*(BCz1-x[2])/sf));
          amrex::Real eta_zmax = 0.5*(1.0 - tanh(2.0*(x[2]-BCz2)/sf));
          eta = eta_xmin*eta_xmax*eta_zmin*eta_zmax;
          if ( wRad > Rad ) {
              amrex::Real slengthW = prob_parm.wall_Lscale/(prob_parm.wall_aspect_ratio/(2.0 + 2.0*prob_parm.wall_aspect_ratio));
              amrex::Real llengthW = prob_parm.wall_aspect_ratio*slengthW;
              amrex::Real wBCx1 = -slengthW/2.0;
              amrex::Real wBCx2 =  slengthW/2.0;
              amrex::Real wBCz1 = -llengthW/2.0;
              amrex::Real wBCz2 =  llengthW/2.0;

              amrex::Real etaW_xmin = 0.5*(1.0 - tanh(2.0*(wBCx1-x[0])/wsf));
              amrex::Real etaW_xmax = 0.5*(1.0 - tanh(2.0*(x[0]-wBCx2)/wsf));
              amrex::Real etaW_zmin = 0.5*(1.0 - tanh(2.0*(wBCz1-x[2])/wsf));
              amrex::Real etaW_zmax = 0.5*(1.0 - tanh(2.0*(x[2]-wBCz2)/wsf));
              etaW = etaW_xmin*etaW_xmax*etaW_zmin*etaW_zmax;
          }
      } else if ( prob_parm.nside > 20 ) {
          Rad = 2.0*prob_parm.Lscale;
          rad = std::sqrt((x[0]-splitx)*(x[0]-splitx) + (x[1]-splity)*(x[1]-splity));
          eta = 0.5*(1.0 + tanh(0.25*(Rad/sf)*(Rad/rad - rad/Rad)));
          if ( wRad > Rad ) {
              wRad = 2.0*prob_parm.wall_Lscale;
              etaW = 0.5*(1.0 + tanh(0.25*(wRad/wsf)*(wRad/rad - rad/wRad)));
          }
      }
#endif
      if ( wRad > Rad ) {
          etaW = etaW - eta;
      }
      amrex::Real etaT0 = 0.5*(1.0 - tanh(2.0*(ht0-time)/hsf));
      amrex::Real etaTf = 0.5*(1.0 - tanh(2.0*(time-htf)/hsf));
      amrex::Real etaT = etaT0*etaTf;

#if ( AMREX_SPACEDIM == 2 )
      s_ext[VELX] = -(1.0-eta)*prob_parm.ce;
      s_ext[VELY] = eta*prob_parm.Vin + (1.0-eta-etaW)*prob_parm.Vcoflow;
#elif ( AMREX_SPACEDIM == 3 )
      if ( prob_parm.nside == 4 ) { 
          s_ext[VELX] = -(x[0]/rad)*(1.0-eta)*prob_parm.ce/rad;
          s_ext[VELY] = eta*prob_parm.Vin + (1.0-eta-etaW)*prob_parm.Vcoflow;
          s_ext[VELZ] = -(x[2]/rad)*(1.0-eta)*prob_parm.ce/rad;;
      } else if ( prob_parm.nside > 20 ) {
          s_ext[VELX] = -(x[0]/rad)*(1.0-eta)*prob_parm.ce/rad;
          s_ext[VELY] = -(x[1]/rad)*(1.0-eta)*prob_parm.ce/rad;
          s_ext[VELZ] = eta*prob_parm.Vin + (1.0-eta-etaW)*prob_parm.Vcoflow;
      }
#endif
      amrex::Real Tinlet = etaT*prob_parm.heating_max_temp + (1.0-etaT)*prob_parm.Tin;
      s_ext[TEMP] = eta*Tinlet + etaW*prob_parm.Twall + (1.0-eta-etaW)*prob_parm.Tcoflow;

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
         massfrac[n] = eta*massfrac_jet[n] + (1.0-eta)*massfrac_air[n];
      }     

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
         s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
