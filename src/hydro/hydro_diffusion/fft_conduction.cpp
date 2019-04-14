//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// Written by J. Squire
//========================================================================================
//! \file fft_conduction.cpp
//  \brief implementation of functions in class FFTConduction for applying 1/|k| operator
//         in the Landau Fluid, CGL model

// C/C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>

// Athena++ headers
#include "fft_conduction.hpp"
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"

#define FFTW_PLANNING FFTW_MEASURE

//----------------------------------------------------------------------------------------
//! \fn FFTConductionDriver::FFTConductionDriver(Mesh *pm, ParameterInput *pin)
//  \brief FFTConduction constructor

FFTConduction::FFTConduction(HydroDiffusion *phdif, ParameterInput *pin) {
  
  pmy_hdif_ = phdif;
  pmy_hydro_ = pmy_hdif_->pmy_hydro_;
  pmb_ = pmy_hydro_->pmy_block;
  Mesh *pm = pmb_->pmy_mesh;
  
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTConduction::FFTConduction" << std::endl
    << "FFT not defined." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
#endif // FFT
  // Chosen direction for B0
  xdir = pin->GetInteger("problem","fft_conduct");
  
  // Start by checking decomposition, make sure it is pencil-like along
  // conduction direction. Do this for meshblocks directly (does not allow
  // multiple meshblocks on one process).
  if (pm->use_uniform_meshgen_fn_[X1DIR]==false
      || pm->use_uniform_meshgen_fn_[X2DIR]==false
      || pm->use_uniform_meshgen_fn_[X3DIR]==false) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTConduction::FFTConduction" << std::endl
    << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  
  int64_t npx1 = pm->nrbx1;
  int64_t npx2 = pm->nrbx2;
  int64_t npx3 = pm->nrbx3;
  
  bool incompat = false;
  if ( xdir==1 && npx1 > 1) incompat=true;
  if ( xdir==2 && (npx2 > 1 || pmb_->block_size.nx2==1)) incompat=true;
  if ( xdir==3 && (npx3 > 1 || pmb_->block_size.nx3==1)) incompat=true;
  if ( incompat==true ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTConduction::FFTConduction" << std::endl
    << "Chosen conduction direction x" << xdir << " is not compatible "
    << "with MeshBlock decomposition or grid dimensions." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  
#ifdef FFT
  // Setup fftw plan. Done in one go for each meshblock using "stride" parameter.
  // Quantities to be transformed are copied into in_. This is much more straightforward
  // due to the ghost zones, and annoyance of transforms in the "2" direction of 3D array.
  
  // Use the "r2r" format, which allows real and complex data to be stored in an
  // array the same size as the real array, and is slightly easier to deal with.
  // Fourier space format is real in h[k], and imag in h[n-k], with h[0] and h[n/2] real.
  nx1_ = pmb_->block_size.nx1;
  nx2_ = pmb_->block_size.nx2;
  nx3_ = pmb_->block_size.nx3;
  ntot_ = nx1_*nx2_*nx3_;
  if ( xdir==1 ){
    ibx = 0; iby = 1; ibz = 2;
    xsize_ = pmb_->block_size.x1max - pmb_->block_size.x1min;
    nft_ = nx1_;
    howmany_ = nx2_*nx3_;
    stride_ = 1;
    dist_ = nx1_;
    in_ = (Real *) fftw_malloc(sizeof(Real) * ntot_);
  }
  if ( xdir==2 ){
    ibx = 1; iby =  0; ibz = 2; // Handles both 2D and 3D (unlike 1,2,0)
    xsize_ = pmb_->block_size.x2max - pmb_->block_size.x2min;
    nft_ = nx2_;
    howmany_ = nx1_;
    stride_ = nx1_;
    dist_ = 1;
    in_ = (Real *) fftw_malloc(sizeof(Real) * nx1_*nx2_);
    // It does not seem to be possible to run 2-direction ffts with a
    // single fftw call, so copy into a series of slices to run fft in this case.
  }
  if ( xdir==3 ){
    ibx = 2; iby = 0; ibz = 1;
    xsize_ = pmb_->block_size.x3max - pmb_->block_size.x3min;
    nft_ = nx3_;
    howmany_ = nx1_*nx2_;
    stride_ = nx1_*nx2_;
    dist_ = 1;
    in_ = (Real *) fftw_malloc(sizeof(Real) * ntot_);
  }
  kprl_.NewAthenaArray(nft_);
  ConstructKprl();
  kl_lf_ = phdif->kl_lf;

  // Plans
  fftw_r2r_kind ftype=FFTW_R2HC;
  fftw_r2r_kind btype=FFTW_HC2R;
  forward_ = fftw_plan_many_r2r(1, &nft_, howmany_,
                              in_, NULL, stride_, dist_,
                              in_, NULL, stride_, dist_,
                              &ftype, FFTW_PLANNING);
  backward_ = fftw_plan_many_r2r(1, &nft_, howmany_,
                                 in_, NULL, stride_, dist_,
                                 in_, NULL, stride_, dist_,
                                 &btype, FFTW_PLANNING);

  std::cout << "Finished creating plan for conduction driver\n ";
  
  
#endif // FFT

  
}

FFTConduction::~FFTConduction() {
#ifdef FFT
  fftw_destroy_plan(forward_);
  fftw_destroy_plan(backward_);
  fftw_free(in_);
  kprl_.DeleteAthenaArray();
#endif //FFT
}

//----------------------------------------------------------------------------------------
//! \fn void FFTConduction::Solve(Real csprl, Real rho, Real nu_c)
//  \brief Apply 1/|k| operator to each variable

void FFTConduction::Solve(Real csprl, Real rho, Real nu_c) {
  
  AthenaArray<Real> data, dslice;
  int n;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  
  // Parallel operator
  Real numerprl = 8.*rho*SQR(csprl)/ static_cast<Real>(nft_); //
  Real denomprl = std::sqrt(8.*PI)*csprl;
  Real nuprl = (3.*PI-8.)*nu_c;
  // Perpendicular operator
  Real numerprp = 2.*rho*SQR(csprl)/ static_cast<Real>(nft_);
  Real denomprp = std::sqrt(2.*PI)*csprl;
  Real nuprp = nu_c;
  if ( (xdir==1) || (xdir==3) ){ // FT along x1 or x3
    // ---------------
    // Tprl
    data.InitWithShallowSlice(pmy_hdif_->dprl_cond,4,ICPR,1);
    // Copy data, take fft, apply kernal, take backwards fft, then copy back.
    ExecuteFFT13( 1, data, is, ie, js, je, ks, ke );
    ApplyKernel13(numerprl, denomprl, nuprl);
    ExecuteFFT13(-1, data, is, ie, js, je, ks, ke );
    // ---------------
    // Tprp
    data.InitWithShallowSlice(pmy_hdif_->dprl_cond,4,ICPP,1);
    ExecuteFFT13( 1, data, is, ie, js, je, ks, ke );
    ApplyKernel13(numerprp, denomprp, nuprp);
    ExecuteFFT13(-1, data, is, ie, js, je, ks, ke );
    // ---------------
    // Bmag
    data.InitWithShallowSlice(pmy_hdif_->dprl_cond,4,ICBM,1);
    ExecuteFFT13( 1, data, is, ie, js, je, ks, ke );
    ApplyKernel13(numerprp, denomprp, nuprp);
    ExecuteFFT13(-1, data, is, ie, js, je, ks, ke );
  } else if (xdir==2) { // FT along x2
    // Have to operator separately for each k slice
    // ---------------
    // Tprl
    data.InitWithShallowSlice(pmy_hdif_->dprl_cond,4,ICPR,1);
    for (int k=ks; k<=ke; ++k){
      dslice.InitWithShallowSlice(data,3,k,1);
      ExecuteFFT2( 1, dslice, is, ie, js, je);
      ApplyKernel2(numerprl, denomprl, nuprl);
      ExecuteFFT2(-1, dslice, is, ie, js, je);
    }
    // ---------------
    // Tprp
    data.InitWithShallowSlice(pmy_hdif_->dprl_cond,4,ICPP,1);
    for (int k=ks; k<=ke; ++k){
      dslice.InitWithShallowSlice(data,3,k,1);
      ExecuteFFT2( 1, dslice, is, ie, js, je);
      ApplyKernel2(numerprp, denomprp, nuprp);
      ExecuteFFT2(-1, dslice, is, ie, js, je);
    }
    // ---------------
    // Bmag
    data.InitWithShallowSlice(pmy_hdif_->dprl_cond,4,ICBM,1);
    for (int k=ks; k<=ke; ++k){
      dslice.InitWithShallowSlice(data,3,k,1);
      ExecuteFFT2( 1, dslice, is, ie, js, je);
      ApplyKernel2(numerprp, denomprp, nuprp);
      ExecuteFFT2(-1, dslice, is, ie, js, je);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTConduction::ApplyKernel13()
//  \brief Apply 1/|B0.k| operator to out_ in k space for 1 or 3 direction
void FFTConduction::ApplyKernel13(Real numer, Real denom, Real nu_c) {
  if ( xdir==1 ){
    for (int k=0; k<nx3_; k++) {
      for (int j=0; j<nx2_; j++) {
#pragma omp simd
        for (int i=0; i<nx1_; i++) {
          int64_t id = (nx2_*k + j)*nx1_ + i;
          Real kp = kprl_(i);
          kp = kp > kl_lf_ ? kp : kl_lf_;
          // kl_lf is used as the minimum kprl
          in_[id] *= - numer / ( denom*kp + nu_c );
    }}}
  } else if ( xdir==3 ) {
    for (int k=0; k<nx3_; k++) {
      for (int j=0; j<nx2_; j++) {
#pragma omp simd
        for (int i=0; i<nx1_; i++) {
          int64_t id = (nx2_*k + j)*nx1_ + i;
          Real kp = kprl_(k);
          kp = kp > kl_lf_ ? kp : kl_lf_;
          // kl_lf is used as the minimum kprl
          in_[id] *= - numer / ( denom*kp + nu_c );
    }}}
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FFTConduction::ApplyKernel2(Real numer, Real denom, Real nu_c)
//  \brief Apply 1/|B0.k| operator to out_ in k space for 2 direction
void FFTConduction::ApplyKernel2(Real numer, Real denom, Real nu_c) {
  for (int j=0; j<nx2_; j++) {
#pragma omp simd
    for (int i=0; i<nx1_; i++) {
      int64_t id = j*nx1_ + i;
      Real kp = kprl_(j);
      kp = kp > kl_lf_ ? kp : kl_lf_;
      // kl_lf is used as the minimum kprl
      in_[id] *= - numer / ( denom*kp + nu_c );
  }}
  return;
}


//----------------------------------------------------------------------------------------
//! \fn FFTConductionDriver::ExecuteFFT13(int, AthenaArray<Real>, int,int,int,int,int,int)
//  \brief Copy data and apply transform in x1 or x3 direction
void FFTConduction::ExecuteFFT13(int ftdir, AthenaArray<Real> &data,
                          int is, int ie, int js, int je, int ks, int ke){
  // Runs fft in either forward (ftdir=1) or backward (ftdir=-1) for xdir==1
  // x1 and x3 directions simple, x2 needs a copy of each plane
  int n;
  
  // Apply the transform
  if ( ftdir == 1 ){ // Forward
    // Copy
    n=0;
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
#pragma omp simd
        for (int i=is; i<=ie; i++){
          in_[n++] = data(k,j,i);
    }}}
    // Execute
    fftw_execute(forward_);
  } else if ( ftdir == -1 ) { // Backward
    // Execute
    fftw_execute(backward_);
    // Copy
    n=0;
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
#pragma omp simd
        for (int i=is; i<=ie; i++){
          data(k,j,i) = in_[n++];
    }}}
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn FFTConductionDriver::ExecuteFFT2(int, AthenaArray<Real>, int,int,int,int)
//  \brief Copy data and apply transform in x2 direction
void FFTConduction::ExecuteFFT2(int ftdir, AthenaArray<Real> &dslice,
                                 int is, int ie, int js, int je){
  int n;
  if ( ftdir == 1 ){ // Forward
    // FFT for each k slice
    n=0;
    for (int j=js; j<=je; j++){
#pragma omp simd
      for (int i=is; i<=ie; i++){
        in_[n++] = dslice(j,i);
    }}
    fftw_execute(forward_);
  } else if ( ftdir == -1 ) {
    n=0;
    fftw_execute(backward_);
    for (int j=js; j<=je; j++){
#pragma omp simd
      for (int i=is; i<=ie; i++){
        dslice(j,i) = in_[n++];
    }}
  }
}

//----------------------------------------------------------------------------------------
//! \fn FFTConductionDriver::ConstuctKprl(void)
//  \brief Constructs Kprl along chosen direction
void FFTConduction::ConstructKprl(void) {
  Real k0 = 2.*PI/xsize_;
  // See fftw r2r transform format. Works with both even and odd grid.
  for (int i=0; i<nft_/2+1; i++){
    kprl_(i) = i*k0;
  }
  for (int i=nft_/2+1; i<nft_; i++){
    kprl_(i) = (nft_-i)*k0;
  }
};
