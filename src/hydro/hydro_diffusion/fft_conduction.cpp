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
#include "../../fft/athena_fft.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn FFTConductionDriver::FFTConductionDriver(Mesh *pm, ParameterInput *pin)
//  \brief FFTConductionDriver constructor

FFTConductionDriver::FFTConductionDriver(Mesh *pm, ParameterInput *pin)
: FFTDriver(pm, pin) {

  // initialize using FFTConduction
  
  int igid=Globals::my_rank;
  pmy_fb=new FFTConduction(this, fft_loclist_[igid], igid, fft_mesh_size_, fft_block_size_);
  pmy_fb->SetNormFactor(1./gcnt_);
  
  QuickCreatePlan();
  
  std::cout << "Finished creating plan for conduction driver\n ";
  
}

//----------------------------------------------------------------------------------------
//! \fn void FFTConductionDriver::Solve(int stage)
//  \brief load the data and apply 1/|k| operator to each variable

void FFTConductionDriver::Solve(void) {
  FFTBlock *pfb=pmy_fb;
  AthenaArray<Real> in;
  // Load the source
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
//  for (int igid=nbs;igid<=nbe;igid++) {
//    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
//    if (pmb!=NULL) {
//      in.InitWithShallowSlice(pmb->phydro->u,4,IDN,1);
//      pfb->LoadSource(in, 1, NGHOST, pmb->loc, pmb->block_size);
//    }
//    //    else { // on another process
//    //    }
//  }
  
//  pfb->ExecuteForward();
//  pfb->ApplyKernel(mode);
//  pfb->ExecuteBackward();
//
//  // Return the result
//  for (int igid=nbs;igid<=nbe;igid++) {
//    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
//    if (pmb!=NULL) {
//      pfb->RetrieveResult(pmb->pgrav->phi, 1, NGHOST,
//                          pmb->loc, pmb->block_size);
//    }
//    //    else { // on another process
//    //    }
//  }
  
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTConduction::ApplyKOperator()
//  \brief Apply 1/|k| operator to out_ in k space
void FFTConduction::ApplyKOperator(void) {
  Real pcoeff;
  Real dx1sq=SQR(2*PI/(kNx[0]*dkx[0]));
  Real dx2sq=SQR(2*PI/(kNx[1]*dkx[1]));
  Real dx3sq=SQR(2*PI/(kNx[2]*dkx[2]));
  for (int k=0; k<knx[2]; k++) {
    for (int j=0; j<knx[1]; j++) {
      for (int i=0; i<knx[0]; i++) {
        int64_t gidx = GetGlobalIndex(i,j,k);
        if (gidx == 0) {
          pcoeff = 0.0;
        } else {
          Real kx=(i+kdisp[0]);
          Real ky=(j+kdisp[1]);
          Real kz=(k+kdisp[2]);
          if (kx > 0.5*kNx[0]) kx -= kNx[0];
          if (ky > 0.5*kNx[1]) ky -= kNx[1];
          if (kz > 0.5*kNx[2]) kz -= kNx[2];
          
          kx *= dkx[0];
          ky *= dkx[1];
          kz *= dkx[2];
          pcoeff = -kx*kx;
          if (dim_ > 1) pcoeff -= ky*ky;
          if (dim_ > 2) pcoeff -= kz*kz;
        
          pcoeff = 1.0/pcoeff;
        }
        
        int64_t idx_in=GetIndex(i,j,k,b_in_);
        int64_t idx_out=GetIndex(i,j,k,f_out_);
        in_[idx_in][0] = out_[idx_out][0]*pcoeff;
        in_[idx_in][1] = out_[idx_out][1]*pcoeff;
      }
    }
  }
  return;
}
