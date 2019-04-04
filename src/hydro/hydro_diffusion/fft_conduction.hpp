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
#include "../../fft/athena_fft.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class FFTBlock;
class FFTDriver;
class ConductionSolverTaskList;

//! \class FFTConduction
//  \brief FFT solver for 1/|B0.k| operator on each block

class FFTConduction : public FFTBlock {
public:
  FFTConduction(FFTDriver *pfd, LogicalLocation iloc, int igid,
             RegionSize msize, RegionSize bsize)
  : FFTBlock(pfd, iloc, igid, msize, bsize) {}
  ~FFTConduction() {};
  void ApplyKernel(int mode, Real *params);
};


//! \class FFTDriver
//  \brief FFT solver for 1/|B0.k| operator

class FFTConductionDriver : public FFTDriver{
public:
  FFTConductionDriver(Mesh *pm, ParameterInput *pin);
  ~FFTConductionDriver();
  void Solve( int mode, Real *bmean);
    
  // Task list for sending/receiving boundaries.
//private:
//  ConductionSolverTaskList *cond_tlist_;
};

#endif // HYDRO_HYDRO_DIFFUSION_FFT_CONDUCTION_HPP_