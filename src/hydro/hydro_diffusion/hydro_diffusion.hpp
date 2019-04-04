#ifndef HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
#define HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_diffusion.hpp
//  \brief defines class HydroDiffusion
//  Contains data and functions that implement the diffusion processes

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;
class Coordinates;
class HydroDiffusion;


void ConstViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bc, int is, int ie, int js, int je, int ks, int ke);

void  ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bc, int is, int ie, int js, int je, int ks, int ke);


enum {ISO=0, ANI=1};
enum {FLXEN=0, FLXMU=1};

//! \class HydroDiffusion
//  \brief data and functions for physical diffusion processes in the hydro

class HydroDiffusion {
public:
  HydroDiffusion(Hydro *phyd, ParameterInput *pin);
  ~HydroDiffusion();

  // data
  bool hydro_diffusion_defined;
  Real nu_iso, nu_aniso; // kinetmatic viscosity coeffs
  AthenaArray<Real> visflx[3]; // viscous stress tensor
  AthenaArray<Real> nu; // viscosity array
  bool mirror_limit, firehose_limit; // Only used if nu_aniso defined, limits fluxes
  // As defined Dp = nu_aniso*rho*(BBdV - 1/3*div(V))

  Real kappa_iso, kappa_aniso; // thermal conduction coeff
  AthenaArray<Real> cndflx[3]; // thermal stress tensor
  AthenaArray<Real> kappa; // conduction array
  
  // Storage for parallel gradients in LF calculation
  AthenaArray<Real> dprl_cond;
  
  // Landau-fluid |k| for evaluating heat fluxes
  // If param fft_conduct=1, enables Fourier calculation of 1/|k_prl| ~ 1/|B0.k| in
  //    CalcParallelHeatFluxesFFT. Otherwise, uses 1/|k_prl|=1/kl_lf, with kl_lf set with
  //    kl_landau, as a standard diffusion operator
  Real kl_lf; // See Sharma et al. 2006, equations (13) and (14)
  bool using_fft_for_conduction;

  // functions
  void CalcHydroDiffusionFlux(const AthenaArray<Real> &p, const AthenaArray<Real> &c,
                                    AthenaArray<Real> *flx,
                              const FaceField &b, const AthenaArray<Real> &bcc);
  void AddHydroDiffusionFlux(AthenaArray<Real> *flx_src, AthenaArray<Real> *flx_des);
  void AddHydroDiffusionEnergyFlux(AthenaArray<Real> *flux_src,
                                   AthenaArray<Real> *flux_des);
  void AddHydroDiffusionCGLHeatFlux(AthenaArray<Real> *flux_src,
                                   AthenaArray<Real> *flux_des);
  void ClearHydroFlux(AthenaArray<Real> *flx);
  void SetHydroDiffusivity(AthenaArray<Real> &w, AthenaArray<Real> &bc);
  void NewHydroDiffusionDt(Real &dt_vis, Real &dt_cnd);

  // viscosity
  void ViscousFlux_iso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                             AthenaArray<Real> *flx);
  void ViscousFlux_aniso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                               AthenaArray<Real> *flx,
                         const FaceField &b, const AthenaArray<Real> &bcc);

  // thermal conduction
  void ThermalFlux_iso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                             AthenaArray<Real> *flx);
  void ThermalFlux_aniso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                               AthenaArray<Real> *flx);
  
  // Thermal conduction (heat fluxes) in CGL. This version (without fft) assumes
  // 1/|k_prl|~~kl_lf (kL in Sharma et al. 2006)
  void ThermalFlux_anisoCGL(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                         AthenaArray<Real> *flx,
                    const FaceField &b, const AthenaArray<Real> &bcc);
  
  // Heat fluxes in CGL using the "Landau Fluid" form (Synder et al. 1997)
  // If enabled, this function is called before ThermalFlux_anisoCGL, to allow
  //   boundary conditions between blocks to be sorted in between (see task list).
  void CalcParallelGradientsFFT(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                                 const FaceField &b, const AthenaArray<Real> &bcc);
  // Compute heat fluxes from parallel gradients of T_prl, T_prp, and |B|
  void ThermalFlux_anisoCGLFFT(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                               AthenaArray<Real> *flx,
                                 const FaceField &b, const AthenaArray<Real> &bcc);

private:
  MeshBlock *pmb_;    // ptr to meshblock containing this HydroDiffusion
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroDiffusion
  Coordinates *pco_;  // ptr to coordinates class
  AthenaArray<Real> divv_; // divergence of velocity
  AthenaArray<Real> x1area_,x2area_,x2area_p1_,x3area_,x3area_p1_;
  AthenaArray<Real> vol_;
  AthenaArray<Real> fx_,fy_,fz_;
  AthenaArray<Real> dx1_,dx2_,dx3_;
  AthenaArray<Real> nu_tot_,kappa_tot_;
  AthenaArray<Real> bmagcc_; // Cell centered |B| to speed up LF calculation
  
  // For use in Landau fluid CFL limit
  Real csprl_, rhomean_, nu_c_;
  Real *bhat_mean_;
  AthenaArray<Real> csprl_iloop_, rhomean_iloop_;
  Real sqrt_twopi_, threepi_m_eight_;

  // functions pointer to calculate spatial dependent coefficients
  ViscosityCoeff_t CalcViscCoeff_;
  ConductionCoeff_t CalcCondCoeff_;

  // auxiliary functions to calculate viscous flux
  void Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv);
  void FaceXdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceXdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceXdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
};
#endif // HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
