#ifndef HYDRO_HYDRO_DIFFUSION_FFT_CONDUCTION_HPP_
#define HYDRO_HYDRO_DIFFUSION_FFT_CONDUCTION_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// Written by J. Squire
//========================================================================================
//! \file fftgravity.hpp
//  \brief defines FFTConduction class

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

#ifdef FFT
#include <fftw3.h>
#endif

class Mesh;
class MeshBlock;
class ParameterInput;
class Coordinates;
class Hydro;
class HydroDiffusion;


//! \class FFTConduction
//  \brief FFT solver for 1/|B0.k| operator on each block
// Initialized from meshblock constructor.

class FFTConduction {
public:
  FFTConduction(HydroDiffusion *phdif, ParameterInput *pin);
  ~FFTConduction();

  void Solve(Real csprl_, Real rho_, Real nu_c);
  
  // Multiply by 1/kprl operator in Fourier space
  void ApplyKernel13(Real numer, Real denom, Real nu_c);
  void ApplyKernel2(Real numer, Real denom, Real nu_c);
  
  int xdir;// Direction of the fft
  int ibx, iby, ibz; // Parallel, perp directions
  void ConstructKprl(void);
  
  // Apply Fourier transforms for different xdir choices
  void ExecuteFFT13(int ftdir, AthenaArray<Real> &data,
             int is, int ie, int js, int je, int ks, int ke);
  // Apply Fourier transforms for different xdir choices
  void ExecuteFFT2(int ftdir, AthenaArray<Real> &dslice,
                    int is, int ie, int js, int je);
  
private:
  MeshBlock *pmb_;    // ptr to meshblock containing this
  Hydro *pmy_hydro_;  // ptr to Hydro containing this
  HydroDiffusion *pmy_hdif_;

  Real xsize_;
  int nx1_, nx2_,nx3_, ntot_;
  int nft_, howmany_; // Size of Fourier transform. How many of them
  int stride_, dist_; // Stride and distance between ffts (see fftw)
  // Plans
#ifdef FFT
  fftw_plan forward_, backward_;
  
#endif
  
  Real *in_; // Copy data into *in to take FFT
  AthenaArray<Real> kprl_; // |k_0prl|
  Real kl_lf_; // kl_lf is used as the minimum k for which Fourier operator applies
  
};



#endif // HYDRO_HYDRO_DIFFUSION_FFT_CONDUCTION_HPP_


