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
  // Does not work in 1D with MPI
#ifdef MPI_PARALLEL
  if ( pm->mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTConductionDriver::FFTConductionDriver" << std::endl
    << "Fourier transform in 1D with MPI are not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
#endif
  // initialize using FFTConduction
  int igid=Globals::my_rank;
  pmy_fb=new FFTConduction(this, fft_loclist_[igid], igid, fft_mesh_size_, fft_block_size_);
  pmy_fb->SetNormFactor(1./gcnt_);
  
  QuickCreatePlan();
  
  std::cout << "Finished creating plan for conduction driver\n ";
  
}

FFTConductionDriver::~FFTConductionDriver() {
}

//----------------------------------------------------------------------------------------
//! \fn void FFTConductionDriver::Solve(int stage)
//  \brief load the data and apply 1/|k| operator to each variable

void FFTConductionDriver::Solve(int mode, Real *meandata) {
  FFTBlock *pfb=pmy_fb;
  AthenaArray<Real> in, out;
  
  // Load the source
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
  
  // Apply 1/|B0.k| operator to tprp, tprl, and bmag
  
  //  meandata={|Bx1|,|Bx2|,|Bx3|,rho,cs_prl,nu_c}
  // Mean B, numerator, and denomenators for operator
  Real bm_num_denom[6];
  bm_num_denom[0] = meandata[0];
  bm_num_denom[1] = meandata[1];
  bm_num_denom[2] = meandata[2];
  // Parallel operator
  bm_num_denom[3] = 8.*meandata[3]*SQR(meandata[4]);
  bm_num_denom[4] = std::sqrt(8.*PI)*meandata[4];
  bm_num_denom[5] = (3.*PI-8.)*meandata[5];
  // Tprl
  for (int igid=nbs;igid<=nbe;igid++) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if (pmb!=NULL) {
      in.InitWithShallowSlice(pmb->phydro->phdif->dprl_cond,4,ICPR,1);
      pfb->LoadSource(in, 1, NGHOST, pmb->loc, pmb->block_size);
    }
  }
  pfb->ExecuteForward();
  pfb->ApplyKernel(mode,bm_num_denom);
  pfb->ExecuteBackward();
  // Return the result
  for (int igid=nbs;igid<=nbe;igid++) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if (pmb!=NULL) {
      out.InitWithShallowSlice(pmb->phydro->phdif->dprl_cond,4,ICPR,1);
      pfb->RetrieveResult(out, 1, NGHOST,
                          pmb->loc, pmb->block_size);
    }
  }
  // Perpendicular operator
  bm_num_denom[3] = 2.*meandata[3]*SQR(meandata[4]);
  bm_num_denom[4] = std::sqrt(2.*PI)*meandata[4];
  bm_num_denom[5] = meandata[5];
  // Tprp
  for (int igid=nbs;igid<=nbe;igid++) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if (pmb!=NULL) {
      in.InitWithShallowSlice(pmb->phydro->phdif->dprl_cond,4,ICPP,1);
      pfb->LoadSource(in, 1, NGHOST, pmb->loc, pmb->block_size);
    }
  }
  pfb->ExecuteForward();
  pfb->ApplyKernel(mode,bm_num_denom);
  pfb->ExecuteBackward();
  // Return the result
  for (int igid=nbs;igid<=nbe;igid++) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if (pmb!=NULL) {
      out.InitWithShallowSlice(pmb->phydro->phdif->dprl_cond,4,ICPP,1);
      pfb->RetrieveResult(out, 1, NGHOST,
                          pmb->loc, pmb->block_size);
    }
  }
  // |B|
  for (int igid=nbs;igid<=nbe;igid++) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if (pmb!=NULL) {
      in.InitWithShallowSlice(pmb->phydro->phdif->dprl_cond,4,ICBM,1);
      pfb->LoadSource(in, 1, NGHOST, pmb->loc, pmb->block_size);
    }
  }
  pfb->ExecuteForward();
  pfb->ApplyKernel(mode,bm_num_denom);
  pfb->ExecuteBackward();
  // Return the result
  for (int igid=nbs;igid<=nbe;igid++) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if (pmb!=NULL) {
      out.InitWithShallowSlice(pmb->phydro->phdif->dprl_cond,4,ICBM,1);
      pfb->RetrieveResult(out, 1, NGHOST,
                          pmb->loc, pmb->block_size);
    }
  }
    
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTConduction::ApplyKernel()
//  \brief Apply 1/|B0.k| operator to out_ in k space
void FFTConduction::ApplyKernel(int mode, Real *md) {
  // Operator is -md[3]/(md[4]*kprl + md[5]) where kprl=kx*md[0] + ky*md[1] + ky*md[2];
  //  Or, if |kprl| is larger than the box, just -1
  Real pcoeff;
  Real boxk = std::min(dkx[0],dkx[1]);
  boxk = 0.5*std::min(boxk,dkx[2]);
//  std::cout << "I'm in ApplyKernal" << std::endl;
  for (int k=0; k<knx[2]; k++) {
    for (int j=0; j<knx[1]; j++) {
      for (int i=0; i<knx[0]; i++) {
        int64_t gidx = GetGlobalIndex(i,j,k);
        
        Real kx=(i+kdisp[0]);
        Real ky=(j+kdisp[1]);
        Real kz=(k+kdisp[2]);
        if (kx > 0.5*kNx[0]) kx -= kNx[0];
        if (ky > 0.5*kNx[1]) ky -= kNx[1];
        if (kz > 0.5*kNx[2]) kz -= kNx[2];
        
        kx *= dkx[0];
        ky *= dkx[1];
        kz *= dkx[2];
        pcoeff = fabs(kx*md[0] + ky*md[1] + kz*md[2]);
      
        // No effect on quantities that are constant along the field.
        if (pcoeff >= boxk)
          pcoeff = - md[3]/(md[4]*fabs(pcoeff) + md[5]);
        else
          pcoeff = -1.0;
        
        int64_t idx_in=GetIndex(i,j,k,b_in_);
        int64_t idx_out=GetIndex(i,j,k,f_out_);
        
        in_[idx_in][0] = out_[idx_out][0]*pcoeff;
        in_[idx_in][1] = out_[idx_out][1]*pcoeff;
        
      }
    }
  }
//  std::cout << "\n\n";
  return;
}


