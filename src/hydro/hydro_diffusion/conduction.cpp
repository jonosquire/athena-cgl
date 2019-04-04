//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C++ headers
#include <cmath>
#include <iostream>   // endl
#include <sstream>    // sstream
#include <string>     // c_str()

// Athena++ headers
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"
#include "fft_conduction.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

// Declarations
static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);

//---------------------------------------------------------------------------------------
// Calculate isotropic thermal conduction

void HydroDiffusion::ThermalFlux_iso(const AthenaArray<Real> &prim,
              const AthenaArray<Real> &cons, AthenaArray<Real> *cndflx) {
  
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real kappaf, denf, dTdx, dTdy, dTdz;

  // i-direction
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) {
      if(pmb_->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        dTdx = (prim(IPR,k,j,i)/prim(IDN,k,j,i) - prim(IPR,k,j,i-1)/
                prim(IDN,k,j,i-1))/pco_->dx1v(i-1);
        x1flux(k,j,i) -= kappaf*denf*dTdx;
      }
    }
  }

  // j-direction
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx3 == 1) // 2D
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j-1,i));
          dTdy = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k,j-1,i)/
                    prim(IDN,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1);
          x2flux(k,j,i) -= kappaf*denf*dTdy;
        }
      }
    }
  } // zero flux for 1D

  // k-direction
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) // 2D or 3D
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D
      il=is-1, iu=ie+1;
  }
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k-1,j,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k-1,j,i));
          dTdz = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k-1,j,i)/
                   prim(IDN,k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
          x3flux(k,j,i) -= kappaf*denf*dTdz;
        }
      }
    }
  } // zero flux for 1D/2D

  return;
}


//---------------------------------------------------------------------------------------
// Calculate anisotropic thermal conduction

void HydroDiffusion::ThermalFlux_aniso(const AthenaArray<Real> &p,
                 const AthenaArray<Real> &c, AthenaArray<Real> *flx) {
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::CalcParallelGradientsFFT
//  \brief Calculate heat fluxes using FFT method. Called before ThermalFlux_anisoCGL
//   Computes parallel gradients at cell centeres, then applies fft
void HydroDiffusion::CalcParallelGradientsFFT(const AthenaArray<Real> &prim,
                                               const AthenaArray<Real> &cons,
                                              const FaceField &b, const AthenaArray<Real> &bcc){
  

  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real bfloor = pmb_->peos->GetBFieldFloor();
  Real nu_c = pmb_->peos->GetCollisionFreq();
  
  // Compute cell centered |B|, from which we will take derivatives
  // Also compute average of |B| in each direction, to use for computing
  il=is-1, iu=ie+1, jl=js, ju=je, kl=ks, ku=ke;
  if(pmb_->block_size.nx2 > 1) {
    if(pmb_->block_size.nx3 == 1) // 2D
      jl=js-1, ju=je+1;
    else // 3D
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        bmagcc_(k,j,i) = std::sqrt( bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                                   bcc(IB2,k,j,i)*bcc(IB2,k,j,i) +
                                   bcc(IB3,k,j,i)*bcc(IB3,k,j,i) );
        bmagcc_(k,j,i) = (bmagcc_(k,j,i) > bfloor) ?  bmagcc_(k,j,i) : bfloor;
      }
    }
  }
  // Compute mean quantities to use in Fourier operator.
  // meandata={|bhx1|,|bhx2|,|bhx3|,rho,cs_prl,nu_c}
  Real meandata[6] = {0}, gmean[6];
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        meandata[0] += bcc(IB1,k,j,i);
        meandata[1] += bcc(IB2,k,j,i);
        meandata[2] += bcc(IB3,k,j,i);
        meandata[3] += prim(IDN,k,j,i);
        meandata[4] += prim(IPR,k,j,i)/prim(IDN,k,j,i);
      }
    }
  }
  
#ifdef MPI_PARALLEL
  // Sum the mean field, cs_prl, rho over all processors
  int mpierr = MPI_Allreduce(meandata, gmean, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    std::stringstream msg;
    msg << "[normalize]: MPI_Allreduce error = "
    << mpierr << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int n=0; n<5; n++) meandata[n]=gmean[n];
#endif // MPI_PARALLEL
  for (int n=0; n<5; n++) meandata[n] /= (pmb_->pmy_mesh->mesh_size.nx1 *
                                          pmb_->pmy_mesh->mesh_size.nx2 *
                                          pmb_->pmy_mesh->mesh_size.nx3);
  meandata[4] = std::sqrt(meandata[4]);
  meandata[5] = nu_c;
  // meandata[0-2] is a unit vector direction for 1/|b0.k|
  Real meanb_mag = std::sqrt( SQR(meandata[0]) + SQR(meandata[1]) + SQR(meandata[2]));
  for (int n=0;n<3;n++) meandata[n] /= meanb_mag;
  meandata[0]=0;meandata[1]=1.;meandata[2]=0.;
//  std::cout << "bhat=" <<meandata[0]<<","<<meandata[1]<<","<<meandata[2]<<"\n";
  
  // Need details of 1/|b0.k| operator for CFL limit
  for (int n=0;n<3;n++) bhat_mean_[n] = meandata[n];
  csprl_ = meandata[4];
  rhomean_ = meandata[3];
  nu_c_ = nu_c;
  
  // Compute cell centered derivatives along the field lines
  if (pmb_->block_size.nx3 > 1) {// 3D
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          
          Real bhx = bcc(IB1,k,j,i)/bmagcc_(k,j,i);
          Real bhy = bcc(IB2,k,j,i)/bmagcc_(k,j,i);
          Real bhz = bcc(IB3,k,j,i)/bmagcc_(k,j,i);
          
          Real tprlcc = prim(IPR,k,j,i) / prim(IDN,k,j,i);
          Real tprpcc = prim(IPP,k,j,i) / prim(IDN,k,j,i);
          
          Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
          Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
          Real dz = 0.5*(pco_->dx1v(k-1) + pco_->dx1v(k));
          
          Real dtprl_dx = vanleer( prim(IPR,k,j,i+1)/prim(IDN,k,j,i+1) - tprlcc,
                                  -prim(IPR,k,j,i-1)/prim(IDN,k,j,i-1) + tprlcc) / dx;
          Real dtprp_dx = vanleer( prim(IPP,k,j,i+1)/prim(IDN,k,j,i+1) - tprpcc,
                                  -prim(IPP,k,j,i-1)/prim(IDN,k,j,i-1) + tprpcc) / dx;
          Real dbmag_dx = vanleer( bmagcc_(k,j,i+1) - bmagcc_(k,j,i),
                                  -bmagcc_(k,j,i-1) + bmagcc_(k,j,i)) / dx;
          
          Real dtprl_dy = vanleer( prim(IPR,k,j+1,i)/prim(IDN,k,j+1,i) - tprlcc,
                                  -prim(IPR,k,j-1,i)/prim(IDN,k,j-1,i) + tprlcc) / dy;
          Real dtprp_dy = vanleer( prim(IPP,k,j+1,i)/prim(IDN,k,j+1,i) - tprpcc,
                                  -prim(IPP,k,j-1,i)/prim(IDN,k,j-1,i) + tprpcc) / dy;
          Real dbmag_dy = vanleer( bmagcc_(k,j+1,i) - bmagcc_(k,j,i),
                                  -bmagcc_(k,j-1,i) + bmagcc_(k,j,i)) / dy;
          
          Real dtprl_dz = vanleer( prim(IPR,k+1,j,i)/prim(IDN,k+1,j,i) - tprlcc,
                                  -prim(IPR,k-1,j,i)/prim(IDN,k-1,j,i) + tprlcc) / dz;
          Real dtprp_dz = vanleer( prim(IPP,k+1,j,i)/prim(IDN,k+1,j,i) - tprpcc,
                                  -prim(IPP,k-1,j,i)/prim(IDN,k-1,j,i) + tprpcc) / dz;
          Real dbmag_dz = vanleer( bmagcc_(k+1,j,i) - bmagcc_(k,j,i),
                                  -bmagcc_(k-1,j,i) + bmagcc_(k,j,i)) / dz;
          
          dprl_cond(ICPR,k,j,i) = bhx*dtprl_dx + bhy*dtprl_dy + bhz*dtprl_dz;
          dprl_cond(ICPP,k,j,i) = bhx*dtprp_dx + bhy*dtprp_dy + bhz*dtprp_dz;
          dprl_cond(ICBM,k,j,i) = bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz;
          
    }}}
  } else if (pmb_->block_size.nx2 > 1) {// 2D
    int k=0;
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        
        Real bhx = bcc(IB1,k,j,i)/bmagcc_(k,j,i);
        Real bhy = bcc(IB2,k,j,i)/bmagcc_(k,j,i);
        Real bhz = bcc(IB3,k,j,i)/bmagcc_(k,j,i);
        
        Real tprlcc = prim(IPR,k,j,i) / prim(IDN,k,j,i);
        Real tprpcc = prim(IPP,k,j,i) / prim(IDN,k,j,i);
        
        Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
        Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
        
        Real dtprl_dx = vanleer( prim(IPR,k,j,i+1)/prim(IDN,k,j,i+1) - tprlcc,
                                -prim(IPR,k,j,i-1)/prim(IDN,k,j,i-1) + tprlcc) / dx;
        Real dtprp_dx = vanleer( prim(IPP,k,j,i+1)/prim(IDN,k,j,i+1) - tprpcc,
                                -prim(IPP,k,j,i-1)/prim(IDN,k,j,i-1) + tprpcc) / dx;
        Real dbmag_dx = vanleer( bmagcc_(k,j,i+1) - bmagcc_(k,j,i),
                                -bmagcc_(k,j,i-1) + bmagcc_(k,j,i)) / dx;
        
        Real dtprl_dy = vanleer( prim(IPR,k,j+1,i)/prim(IDN,k,j+1,i) - tprlcc,
                                -prim(IPR,k,j-1,i)/prim(IDN,k,j-1,i) + tprlcc) / dy;
        Real dtprp_dy = vanleer( prim(IPP,k,j+1,i)/prim(IDN,k,j+1,i) - tprpcc,
                                -prim(IPP,k,j-1,i)/prim(IDN,k,j-1,i) + tprpcc) / dy;
        Real dbmag_dy = vanleer( bmagcc_(k,j+1,i) - bmagcc_(k,j,i),
                                -bmagcc_(k,j-1,i) + bmagcc_(k,j,i)) / dy;
        
        dprl_cond(ICPR,k,j,i) = bhx*dtprl_dx + bhy*dtprl_dy ;
        dprl_cond(ICPP,k,j,i) = bhx*dtprp_dx + bhy*dtprp_dy ;
        dprl_cond(ICBM,k,j,i) = bhx*dbmag_dx + bhy*dbmag_dy ;
        
      }}
  } else { // 1D
    int k=0, j=0;
#pragma omp simd
    for (int i=is; i<=ie; ++i) {
      
      Real bhx = bcc(IB1,k,j,i)/bmagcc_(k,j,i);
      Real bhy = bcc(IB2,k,j,i)/bmagcc_(k,j,i);
      Real bhz = bcc(IB3,k,j,i)/bmagcc_(k,j,i);
      
      Real tprlcc = prim(IPR,k,j,i) / prim(IDN,k,j,i);
      Real tprpcc = prim(IPP,k,j,i) / prim(IDN,k,j,i);
      
      Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
      
      Real dtprl_dx = vanleer( prim(IPR,k,j,i+1)/prim(IDN,k,j,i+1) - tprlcc,
                              -prim(IPR,k,j,i-1)/prim(IDN,k,j,i-1) + tprlcc) / dx;
      Real dtprp_dx = vanleer( prim(IPP,k,j,i+1)/prim(IDN,k,j,i+1) - tprpcc,
                              -prim(IPP,k,j,i-1)/prim(IDN,k,j,i-1) + tprpcc) / dx;
      Real dbmag_dx = vanleer( bmagcc_(k,j,i+1) - bmagcc_(k,j,i),
                              -bmagcc_(k,j,i-1) + bmagcc_(k,j,i)) / dx;

      dprl_cond(ICPR,k,j,i) = bhx*dtprl_dx ;
      dprl_cond(ICPP,k,j,i) = bhx*dtprp_dx ;
      dprl_cond(ICBM,k,j,i) = bhx*dbmag_dx ;
      
    }
  }
  
  //  std::cout << "average quantites " << meandata[0] << " "<< meandata[1] << " "<< meandata[2] << " "<< meandata[3] << " "<< meandata[4] << " " << meandata[5] << "\n";
  
  // Apply -2cs_prl^2/(sqrt(2pi)*cs_prl*|B0.k| + nu_c) operator for Tprp and B
  // -8cs_prl^2/(sqrt(8pi)*cs_prl*|B0.k| + (3pi-8)nu_c) operator for Tprl
//  bool print = true;
//  if (print){
//    int k=0, j=0;
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (Globals::my_rank == 0) {
//      std::cout << "Before dprl_tprl p1 " <<  meandata[0] << " " <<  meandata[1] << " " <<  meandata[2] <<" " <<  meandata[3] << " " <<  meandata[4] << " " <<  meandata[5] <<  "\n";
//      for (int i=is-1; i<=ie+1; ++i) {
//        std::cout << dprl_cond(ICPR,k,j,i) << ", ";
//      }
//      std::cout << "\n";
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (Globals::my_rank == 1) {
//      std::cout << "Before dprl_tprl p2 " <<  meandata[0] << " " <<  meandata[1] << " " <<  meandata[2] << "\n";
//      for (int i=is-1; i<=ie+1; ++i) {
//        std::cout << dprl_cond(ICPR,k,j,i) << ", ";
//      }
//      std::cout << "\n";
//    }
//  }
  
  int mode = 0;
  pmb_->pmy_mesh->pfcondd->Solve( mode, meandata);
  
  // Code to test just using kL -- to delete (or could delete kL code...)
//  Real prlmult[2], prpmult[2];
//  prpmult[0] = 2.*meandata[3]*SQR(meandata[4]);
//  prpmult[1] = std::sqrt(2.*PI)*meandata[4];
//  prlmult[0] = 8.*meandata[3]*SQR(meandata[4]);
//  prlmult[1] = std::sqrt(8.*PI)*meandata[4];
//  Real kL = 100.;
//
//
//  il=is-1, iu=ie+1, jl=js, ju=je, kl=ks, ku=ke;
//  for (int k=kl; k<=ku; ++k) {
//    for (int j=jl; j<=ju; ++j) {
//#pragma omp simd
//      for (int i=il; i<=iu; ++i) {
//        dprl_cond(ICPR,k,j,i) *= - prlmult[0]/(prlmult[1]*kL);
//        dprl_cond(ICPP,k,j,i) *= - prpmult[0]/(prpmult[1]*kL);
//        dprl_cond(ICBM,k,j,i) *= - prpmult[0]/(prpmult[1]*kL);
//      }
//    }
//  }
  // Compute
  
//  if (print){
//    int k=0, j=0;
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (Globals::my_rank == 0) {
//      std::cout << "After dprl_tprl p1 ";
//      for (int i=is-1; i<=ie+1; ++i) {
//        std::cout << dprl_cond(ICPR,k,j,i) << ", ";
//      }
//      std::cout << "\n";
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (Globals::my_rank == 1) {
//      std::cout << "After dprl_tprl p2 ";
//      for (int i=is-1; i<=ie+1; ++i) {
//        std::cout << dprl_cond(ICPR,k,j,i) << ", ";
//      }
//      std::cout << "\n";
//    }
//  }

  // After solve computes operators using FFT, boundary calls are handled through the
  // ConductionBoundaryValues class

  return;
  
}

//----------------------------------------------------------------------------------------
// Parallel heat fluxes (of perpendicular and parallel heat) in CGL. This version is
// used if the paralle gradients are already known from CalcParallelGradientsFFT
void HydroDiffusion::ThermalFlux_anisoCGLFFT(const AthenaArray<Real> &prim,const AthenaArray<Real> &cons,
                                             AthenaArray<Real> *flx,
                                             const FaceField &b, const AthenaArray<Real> &bcc){
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  

//  bool print = true;
//  if (print){
//    int k=0, j=0;
////    MPI_Barrier(MPI_COMM_WORLD);
//    if (Globals::my_rank == 0) {
//      std::cout << "After bvals dprl_tprl p1 ";
//      for (int i=is-1; i<=ie+1; ++i) {
//        std::cout << dprl_cond(ICPR,k,j,i) << ", ";
//      }
//      std::cout << "\n";
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (Globals::my_rank == 1) {
//      std::cout << "After bvals dprl_tprl p2 ";
//      for (int i=is-1; i<=ie+1; ++i) {
//        std::cout << dprl_cond(ICPR,k,j,i) << ", ";
//      }
//      std::cout << "\n";
//    }
//  }
  
  // Form heat fluxes in each direction
  // i-direction
  jl=js, ju=je, kl=ks, ku=ke;
  if(pmb_->block_size.nx2 > 1) {
    if(pmb_->block_size.nx3 == 1) // 2D
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    else // 3D
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        Real qprl = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j,i-1));
        Real qprp = 0.5*(dprl_cond(ICPP,k,j,i) + dprl_cond(ICPP,k,j,i-1));
//        Real qprl = vanleer(dprl_cond(ICPR,k,j,i), dprl_cond(ICPR,k,j,i-1));
//        Real qprp = vanleer(dprl_cond(ICPP,k,j,i), dprl_cond(ICPP,k,j,i-1));
//        qprp -= 0.5*prim(IPP,k,j,i)/prim(IDN,k,j,i) * (1 - prim(IPP,k,j,i)/prim(IPR,k,j,i) )/
//                bmagcc_(k,j,i) * dprl_cond(ICBM,k,j,i);
//        qprp -= 0.5*prim(IPP,k,j,i-1)/prim(IDN,k,j,i-1)*(1-prim(IPP,k,j,i-1)/prim(IPR,k,j,i-1) )/
//                bmagcc_(k,j,i-1) * dprl_cond(ICBM,k,j,i-1);
        
        Real bmag = 0.5*(bmagcc_(k,j,i) + bmagcc_(k,j,i-1));
        
        x1flux(FLXEN,k,j,i) += b.x1f(k,j,i) / bmag * (qprp + qprl/2.);
        x1flux(FLXMU,k,j,i) += b.x1f(k,j,i) * qprp / (bmag * bmag);
        
      }
    }
  }
  
  // j-direction
  il=is, iu=ie, kl=ks, ku=ke;
  if(pmb_->block_size.nx3 == 1) // 2D
    il=is-1, iu=ie+1, kl=ks, ku=ke;
  else // 3D
    il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
//          Real qprl = vanleer(dprl_cond(ICPR,k,j,i), dprl_cond(ICPR,k,j-1,i));
//          Real qprp = vanleer(dprl_cond(ICPP,k,j,i), dprl_cond(ICPP,k,j-1,i));
          Real qprl = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j-1,i));
          Real qprp = 0.5*(dprl_cond(ICPP,k,j,i) + dprl_cond(ICPP,k,j-1,i));
//          qprp -= 0.5*prim(IPP,k,j,i)/prim(IDN,k,j,i) * (1 - prim(IPP,k,j,i)/prim(IPR,k,j,i) )/
//                  bmagcc_(k,j,i) * dprl_cond(ICBM,k,j,i);
//          qprp -= 0.5*prim(IPP,k,j-1,i)/prim(IDN,k,j-1,i)*(1-prim(IPP,k,j-1,i)/prim(IPR,k,j-1,i) )/
//                  bmagcc_(k,j-1,i) * dprl_cond(ICBM,k,j-1,i);
          
          Real bmag = 0.5*(bmagcc_(k,j,i) + bmagcc_(k,j-1,i));
          
          x2flux(FLXEN,k,j,i) += b.x2f(k,j,i) / bmag * (qprp + qprl/2.);
          x2flux(FLXMU,k,j,i) += b.x2f(k,j,i) * qprp / (bmag * bmag);
          
        }
      }
    }
  } // zero flux for 1D
  
  // k-direction
  il=is-1, iu=ie+1, jl=js-1, ju=je+1;

  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real qprl = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k-1,j,i));
          Real qprp = 0.5*(dprl_cond(ICPP,k,j,i) + dprl_cond(ICPP,k-1,j,i));
          qprp -= 0.5*prim(IPP,k,j,i)/prim(IDN,k,j,i) * (1 - prim(IPP,k,j,i)/prim(IPR,k,j,i) )/
                  bmagcc_(k,j,i) * dprl_cond(ICBM,k,j,i);
          qprp -= 0.5*prim(IPP,k-1,j,i)/prim(IDN,k-1,j,i)*(1-prim(IPP,k-1,j,i)/prim(IPR,k-1,j,i) )/
                  bmagcc_(k-1,j,i) * dprl_cond(ICBM,k-1,j,i);
          
          Real bmag = 0.5*(bmagcc_(k,j,i) + bmagcc_(k-1,j,i));
          
          x3flux(FLXEN,k,j,i) += b.x3f(k,j,i) / bmag * (qprp + qprl/2.);
          x3flux(FLXMU,k,j,i) += b.x3f(k,j,i) * qprp / (bmag * bmag);
          
        }
      }
    }
  } // zero flux for 1D/2D
    
  return;
}


//----------------------------------------------------------------------------------------
// Parallel heat fluxes (of perpendicular and parallel heat) in CGL. If kl_lf>0, this
// function is all that is evaluated.
void HydroDiffusion::ThermalFlux_anisoCGL(const AthenaArray<Real> &prim,const AthenaArray<Real> &cons,
                          AthenaArray<Real> *flx,
                          const FaceField &b, const AthenaArray<Real> &bcc){
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real bfloor = pmb_->peos->GetBFieldFloor();
  Real nu_c = pmb_->peos->GetCollisionFreq();
  
  // To compute CFL conditions, csprl_iloop will be block maximum of
  int ncells1 = pmb_->block_size.nx1;
  nu_c_ = nu_c;
  
  // Store pprp/rho and pprl/rho to save lots of computation
  il=is-1, iu=ie+1, jl=js, ju=je, kl=ks, ku=ke;
  if(pmb_->block_size.nx2 > 1) {
    if(pmb_->block_size.nx3 == 1) // 2D
      jl=js-1, ju=je+1;
    else // 3D
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        // dprl_cond is useful storage for Tprl and Tprp here, saves computation
        dprl_cond(ICPR,k,j,i) = prim(IPR,k,j,i) / prim(IDN,k,j,i);
        dprl_cond(ICPP,k,j,i) = prim(IPP,k,j,i) / prim(IDN,k,j,i);
        bmagcc_(k,j,i) = std::sqrt( bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                                   bcc(IB2,k,j,i)*bcc(IB2,k,j,i) +
                                   bcc(IB3,k,j,i)*bcc(IB3,k,j,i) );
        bmagcc_(k,j,i) = (bmagcc_(k,j,i) > bfloor) ?  bmagcc_(k,j,i) : bfloor;
      }
    }
  }
  
  
  // ----------------------------------------------------------------- //
  // i direction
  jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  
  if (pmb_->block_size.nx3 > 1) {// 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          // Compute bhat at the i-1/2 interface
          Real bx = b.x1f(k,j,i);
          Real by = 0.5*(bcc(IB2,k,j,i) + bcc(IB2,k,j,i-1));
          Real bz = 0.5*(bcc(IB3,k,j,i) + bcc(IB3,k,j,i-1));
          Real bmag = std::sqrt(bx*bx + by*by + bz*bz);
          bmag = (bmag > bfloor) ?  bmag : bfloor;
          Real bhx = bx/bmag;
          Real bhy = by/bmag;
          Real bhz = bz/bmag;
          
          // rho, pprl and prp to multiply derivatives
          Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
          Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k,j,i-1));
          Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k,j,i-1));
          Real csprl2 = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j,i-1));
          Real csprl = std::sqrt(csprl2);
          csprl_iloop_(i) = csprl;
          rhomean_iloop_(i) = rho;
          
          Real dx = pco_->dx1v(i-1);
          Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
          Real dz = 0.5*(pco_->dx3v(k-1) + pco_->dx3v(k));
          
          // x gradients at i-1/2 j k
          Real dplrho_dx = (dprl_cond(ICPR,k,j,i) - dprl_cond(ICPR,k,j,i-1) ) / dx;
          Real dpprho_dx = (dprl_cond(ICPP,k,j,i) - dprl_cond(ICPP,k,j,i-1) ) / dx;
          Real dbmag_dx = (bmagcc_(k,j,i) - bmagcc_(k,j,i-1) ) / dx;
          
          // y gradients at i-1/2 j k
          Real dplrho_dy = limiter4(dprl_cond(ICPR,k,j+1,i  ) - dprl_cond(ICPR,k,j  ,i  ),
                                    dprl_cond(ICPR,k,j  ,i  ) - dprl_cond(ICPR,k,j-1,i  ),
                                    dprl_cond(ICPR,k,j+1,i-1) - dprl_cond(ICPR,k,j  ,i-1),
                                    dprl_cond(ICPR,k,j  ,i-1) - dprl_cond(ICPR,k,j-1,i-1));
          dplrho_dy /= dy;
          Real dpprho_dy = limiter4(dprl_cond(ICPP,k,j+1,i  ) - dprl_cond(ICPP,k,j  ,i  ),
                                    dprl_cond(ICPP,k,j  ,i  ) - dprl_cond(ICPP,k,j-1,i  ),
                                    dprl_cond(ICPP,k,j+1,i-1) - dprl_cond(ICPP,k,j  ,i-1),
                                    dprl_cond(ICPP,k,j  ,i-1) - dprl_cond(ICPP,k,j-1,i-1));
          dpprho_dy /= dy;
          Real dbmag_dy = limiter4( bmagcc_(k,j+1,i  ) - bmagcc_(k,j  ,i  ),
                                   bmagcc_(k,j  ,i  ) - bmagcc_(k,j-1,i  ),
                                   bmagcc_(k,j+1,i-1) - bmagcc_(k,j  ,i-1),
                                   bmagcc_(k,j  ,i-1) - bmagcc_(k,j-1,i-1));
          dbmag_dy /= dy;
          
          // z gradients at i-1/2 j k
          Real dplrho_dz = limiter4(dprl_cond(ICPR,k+1,j,i  ) - dprl_cond(ICPR,k  ,j,i  ),
                                    dprl_cond(ICPR,k  ,j,i  ) - dprl_cond(ICPR,k-1,j,i  ),
                                    dprl_cond(ICPR,k+1,j,i-1) - dprl_cond(ICPR,k  ,j,i-1),
                                    dprl_cond(ICPR,k  ,j,i-1) - dprl_cond(ICPR,k-1,j,i-1));
          dplrho_dz /= dz;
          Real dpprho_dz = limiter4(dprl_cond(ICPP,k+1,j,i  ) - dprl_cond(ICPP,k  ,j,i  ),
                                    dprl_cond(ICPP,k  ,j,i  ) - dprl_cond(ICPP,k-1,j,i  ),
                                    dprl_cond(ICPP,k+1,j,i-1) - dprl_cond(ICPP,k  ,j,i-1),
                                    dprl_cond(ICPP,k  ,j,i-1) - dprl_cond(ICPP,k-1,j,i-1));
          dpprho_dz /= dz;
          Real dbmag_dz = limiter4( bmagcc_(k+1,j,i  ) - bmagcc_(k  ,j,i  ),
                                   bmagcc_(k  ,j,i  ) - bmagcc_(k-1,j,i  ),
                                   bmagcc_(k+1,j,i-1) - bmagcc_(k  ,j,i-1),
                                   bmagcc_(k  ,j,i-1) - bmagcc_(k-1,j,i-1));
          dbmag_dz /= dz;
          
          // Store maximum sound speed and density to use in computing CFL condition
          Real kappa_prp = 2.*csprl2/(sqrt_twopi_*csprl*kl_lf + nu_c);
          Real kappa_prl = 8.*csprl2/(2.*sqrt_twopi_*csprl*kl_lf + threepi_m_eight_*nu_c);
          
          Real qprp = - kappa_prp * (rho*(bhx*dpprho_dx + bhy*dpprho_dy + bhz*dpprho_dz) -
                                      pprp*(1.-pprp/pprl)*
                                     (bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz)/bmag);
          Real qprl = - kappa_prl * rho * (bhx*dplrho_dx + bhy*dplrho_dy + bhz*dplrho_dz);
          
          x1flux(FLXEN,k,j,i) += bhx*(qprp + qprl/2.);
          x1flux(FLXMU,k,j,i) += bhx * qprp / bmag;
        }
        // CFL conditions using the same idea as the FFT version
        for (int i=is; i<=ie+1; ++i) {
          csprl_ = std::max(csprl_iloop_(i),csprl_);
          rhomean_ = std::max(rhomean_iloop_(i),rhomean_);
        }
      }}
  } else if (pmb_->block_size.nx2 > 1) { // 2D
    int k=0;
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        // Compute bhat at the i-1/2 interface
        Real bx = b.x1f(k,j,i);
        Real by = 0.5*(bcc(IB2,k,j,i) + bcc(IB2,k,j,i-1));
        Real bz = 0.5*(bcc(IB3,k,j,i) + bcc(IB3,k,j,i-1));
        Real bmag = std::sqrt(bx*bx + by*by + bz*bz);
        bmag = (bmag > bfloor) ?  bmag : bfloor;
        Real bhx = bx/bmag;
        Real bhy = by/bmag;
        Real bhz = bz/bmag;
        
        // rho, pprl and prp to multiply derivatives
        Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
        Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k,j,i-1));
        Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k,j,i-1));
        Real csprl2 = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j,i-1));
        Real csprl = std::sqrt(csprl2);
        csprl_iloop_(i) = csprl;
        rhomean_iloop_(i) = rho;
        
        Real dx = pco_->dx1v(i-1);
        Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
        
        // x gradients at i-1/2 j k
        Real dplrho_dx = (dprl_cond(ICPR,k,j,i) - dprl_cond(ICPR,k,j,i-1) ) / dx;
        Real dpprho_dx = (dprl_cond(ICPP,k,j,i) - dprl_cond(ICPP,k,j,i-1) ) / dx;
        Real dbmag_dx = (bmagcc_(k,j,i) - bmagcc_(k,j,i-1) ) / dx;
        
        // y gradients at i-1/2 j k
        Real dplrho_dy = limiter4(dprl_cond(ICPR,k,j+1,i  ) - dprl_cond(ICPR,k,j  ,i  ),
                                  dprl_cond(ICPR,k,j  ,i  ) - dprl_cond(ICPR,k,j-1,i  ),
                                  dprl_cond(ICPR,k,j+1,i-1) - dprl_cond(ICPR,k,j  ,i-1),
                                  dprl_cond(ICPR,k,j  ,i-1) - dprl_cond(ICPR,k,j-1,i-1));
        dplrho_dy /= dy;
        Real dpprho_dy = limiter4(dprl_cond(ICPP,k,j+1,i  ) - dprl_cond(ICPP,k,j  ,i  ),
                                  dprl_cond(ICPP,k,j  ,i  ) - dprl_cond(ICPP,k,j-1,i  ),
                                  dprl_cond(ICPP,k,j+1,i-1) - dprl_cond(ICPP,k,j  ,i-1),
                                  dprl_cond(ICPP,k,j  ,i-1) - dprl_cond(ICPP,k,j-1,i-1));
        dpprho_dy /= dy;
        Real dbmag_dy = limiter4( bmagcc_(k,j+1,i  ) - bmagcc_(k,j  ,i  ),
                                 bmagcc_(k,j  ,i  ) - bmagcc_(k,j-1,i  ),
                                 bmagcc_(k,j+1,i-1) - bmagcc_(k,j  ,i-1),
                                 bmagcc_(k,j  ,i-1) - bmagcc_(k,j-1,i-1));
        dbmag_dy /= dy;
        
        // Store maximum sound speed and density to use in computing CFL condition
        Real kappa_prp = 2.*csprl2/(sqrt_twopi_*csprl*kl_lf + nu_c);
        Real kappa_prl = 8.*csprl2/(2.*sqrt_twopi_*csprl*kl_lf + threepi_m_eight_*nu_c);
        
        Real qprp = - kappa_prp * (rho*(bhx*dpprho_dx + bhy*dpprho_dy) -
                                    pprp*(1.-pprp/pprl)*
                                   (bhx*dbmag_dx + bhy*dbmag_dy)/bmag);
        Real qprl = - kappa_prl * rho * (bhx*dplrho_dx + bhy*dplrho_dy);
        
        x1flux(FLXEN,k,j,i) += bhx*(qprp + qprl/2.);
        x1flux(FLXMU,k,j,i) += bhx * qprp / bmag;
      }
      // CFL conditions using the same idea as the FFT version
      for (int i=is; i<=ie+1; ++i) {
        csprl_ = std::max(csprl_iloop_(i),csprl_);
        rhomean_ = std::max(rhomean_iloop_(i),rhomean_);
      }
    }
  } else { // 1D
    int k=0, j=0;
#pragma omp simd
    for (int i=is; i<=ie+1; ++i) {
      // Compute bhat at the i-1/2 interface
      Real bx = b.x1f(k,j,i);
      Real by = 0.5*(bcc(IB2,k,j,i) + bcc(IB2,k,j,i-1));
      Real bz = 0.5*(bcc(IB3,k,j,i) + bcc(IB3,k,j,i-1));
      Real bmag = std::sqrt(bx*bx + by*by + bz*bz);
      bmag = (bmag > bfloor) ?  bmag : bfloor;
      Real bhx = bx/bmag;
      Real bhy = by/bmag;
      Real bhz = bz/bmag;
      
      // rho, pprl and prp to multiply derivatives
      Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
      Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k,j,i-1));
      Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k,j,i-1));
      Real csprl2 = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j,i-1));
      Real csprl = std::sqrt(csprl2);
      csprl_iloop_(i) = csprl;
      rhomean_iloop_(i) = rho;
      
      Real dx = pco_->dx1v(i-1);
      
      // x gradients at i-1/2 j k
      Real dplrho_dx = (dprl_cond(ICPR,k,j,i) - dprl_cond(ICPR,k,j,i-1) ) / dx;
      Real dpprho_dx = (dprl_cond(ICPP,k,j,i) - dprl_cond(ICPP,k,j,i-1) ) / dx;
      Real dbmag_dx = (bmagcc_(k,j,i) - bmagcc_(k,j,i-1) ) / dx;
      
      // Store maximum sound speed and density to use in computing CFL condition
      Real kappa_prp = 2.*csprl2/(sqrt_twopi_*csprl*kl_lf + nu_c);
      Real kappa_prl = 8.*csprl2/(2.*sqrt_twopi_*csprl*kl_lf + threepi_m_eight_*nu_c);
      
      Real qprp = - kappa_prp * (rho*(bhx*dpprho_dx ) -
                                pprp*(1.-pprp/pprl)*
                                    (bhx*dbmag_dx )/bmag);
      Real qprl = - kappa_prl * rho * (bhx*dplrho_dx );
      
      x1flux(FLXEN,k,j,i) += bhx*(qprp + qprl/2.);
      x1flux(FLXMU,k,j,i) += bhx * qprp / bmag;
    }
    // CFL conditions using the same idea as the FFT version
    for (int i=is; i<=ie+1; ++i) {
      csprl_ = std::max(csprl_iloop_(i),csprl_);
      rhomean_ = std::max(rhomean_iloop_(i),rhomean_);
    }
  }
  
  // ----------------------------------------------------------------- //
  // j-direction
  il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  
  if(pmb_->block_size.nx3 > 1) { // 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          // Compute bhat at the j-1/2 interface
          Real bx = 0.5*(bcc(IB1,k,j-1,i) + bcc(IB1,k,j,i));
          Real by = b.x2f(k,j,i);
          Real bz = 0.5*(bcc(IB3,k,j-1,i) + bcc(IB3,k,j,i));
          Real bmag = std::sqrt(bx*bx + by*by + bz*bz);
          bmag = (bmag > bfloor) ?  bmag : bfloor;
          Real bhx = bx/bmag;
          Real bhy = by/bmag;
          Real bhz = bz/bmag;
          
          // rho, pprl and prp to multiply derivatives
          Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j-1,i));
          Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k,j-1,i));
          Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k,j-1,i));
          Real csprl2 = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j-1,i));
          Real csprl = std::sqrt(csprl2);
          
          Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
          Real dy = pco_->dx2v(j-1);
          Real dz = 0.5*(pco_->dx3v(k-1) + pco_->dx3v(k));
          
          // x gradients at i j-1/2 k
          Real dplrho_dx = limiter4(dprl_cond(ICPR,k,j  ,i+1) - dprl_cond(ICPR,k,j  ,i  ),
                                    dprl_cond(ICPR,k,j  ,i  ) - dprl_cond(ICPR,k,j  ,i-1),
                                    dprl_cond(ICPR,k,j-1,i+1) - dprl_cond(ICPR,k,j-1,i  ),
                                    dprl_cond(ICPR,k,j-1,i  ) - dprl_cond(ICPR,k,j-1,i-1));
          dplrho_dx /= dx;
          Real dpprho_dx = limiter4(dprl_cond(ICPP,k,j  ,i+1) - dprl_cond(ICPP,k,j  ,i),
                                    dprl_cond(ICPP,k,j  ,i  ) - dprl_cond(ICPP,k,j  ,i-1),
                                    dprl_cond(ICPP,k,j-1,i+1) - dprl_cond(ICPP,k,j-1,i),
                                    dprl_cond(ICPP,k,j-1,i  ) - dprl_cond(ICPP,k,j-1,i-1));
          dpprho_dx /= dx;
          Real dbmag_dx = limiter4( bmagcc_(k,j  ,i+1) - bmagcc_(k,j  ,i),
                                   bmagcc_(k,j  ,i  ) - bmagcc_(k,j  ,i-1),
                                   bmagcc_(k,j-1,i+1) - bmagcc_(k,j-1,i),
                                   bmagcc_(k,j-1,i  ) - bmagcc_(k,j-1,i-1));
          dbmag_dx /= dx;
          
          // x gradients at i-1/2 j k
          Real dplrho_dy = (dprl_cond(ICPR,k,j,i) - dprl_cond(ICPR,k,j-1,i) ) / dy;
          Real dpprho_dy = (dprl_cond(ICPP,k,j,i) - dprl_cond(ICPP,k,j-1,i) ) / dy;
          Real dbmag_dy = (bmagcc_(k,j,i) - bmagcc_(k,j-1,i) ) / dy;
          
          // z gradients at i j-1/2 k
          Real dplrho_dz = limiter4(dprl_cond(ICPR,k+1,j  ,i) - dprl_cond(ICPR,k  ,j  ,i),
                                    dprl_cond(ICPR,k  ,j  ,i) - dprl_cond(ICPR,k-1,j  ,i),
                                    dprl_cond(ICPR,k+1,j-1,i) - dprl_cond(ICPR,k  ,j-1,i),
                                    dprl_cond(ICPR,k  ,j-1,i) - dprl_cond(ICPR,k-1,j-1,i));
          dplrho_dz /= dz;
          Real dpprho_dz = limiter4(dprl_cond(ICPP,k+1,j  ,i) - dprl_cond(ICPP,k  ,j  ,i),
                                    dprl_cond(ICPP,k  ,j  ,i) - dprl_cond(ICPP,k-1,j  ,i),
                                    dprl_cond(ICPP,k+1,j-1,i) - dprl_cond(ICPP,k  ,j-1,i),
                                    dprl_cond(ICPP,k  ,j-1,i) - dprl_cond(ICPP,k-1,j-1,i));
          dpprho_dz /= dz;
          Real dbmag_dz = limiter4(bmagcc_(k+1,j  ,i) - bmagcc_(k  ,j  ,i),
                                   bmagcc_(k  ,j  ,i) - bmagcc_(k-1,j  ,i),
                                   bmagcc_(k+1,j-1,i) - bmagcc_(k  ,j-1,i),
                                   bmagcc_(k  ,j-1,i) - bmagcc_(k-1,j-1,i));
          dbmag_dz /= dz;
          
          Real qprp = -2.*csprl2/(sqrt_twopi_*csprl*kl_lf + nu_c) *
                    (rho*(bhx*dpprho_dx + bhy*dpprho_dy + bhz*dpprho_dz) -
                        pprp*(1.-pprp/pprl)*
                            (bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz)/bmag);
          Real qprl = -8.*csprl2/(2.*sqrt_twopi_*csprl*kl_lf + threepi_m_eight_*nu_c) *
          rho * (bhx*dplrho_dx + bhy*dplrho_dy + bhz*dplrho_dz);
          
          x2flux(FLXEN,k,j,i) += bhy*(qprp + qprl/2.);
          x2flux(FLXMU,k,j,i) += bhy * qprp / bmag;
          
        }}}
  } else if (pmb_->block_size.nx2 > 1) {// 2D
    int k=0;
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        // Compute bhat at the j-1/2 interface
        Real bx = 0.5*(bcc(IB1,k,j-1,i) + bcc(IB1,k,j,i));
        Real by = b.x2f(k,j,i);
        Real bz = 0.5*(bcc(IB3,k,j-1,i) + bcc(IB3,k,j,i));
        Real bmag = std::sqrt(bx*bx + by*by + bz*bz);
        bmag = (bmag > bfloor) ?  bmag : bfloor;
        Real bhx = bx/bmag;
        Real bhy = by/bmag;
        Real bhz = bz/bmag;
        
        // rho, pprl and prp to multiply derivatives
        Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j-1,i));
        Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k,j-1,i));
        Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k,j-1,i));
        Real csprl2 = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k,j-1,i));
        Real csprl = std::sqrt(csprl2);
        
        Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
        Real dy = pco_->dx2v(j-1);
        
        // x gradients at i j-1/2 k
        Real dplrho_dx = limiter4(dprl_cond(ICPR,k,j  ,i+1) - dprl_cond(ICPR,k,j  ,i  ),
                                  dprl_cond(ICPR,k,j  ,i  ) - dprl_cond(ICPR,k,j  ,i-1),
                                  dprl_cond(ICPR,k,j-1,i+1) - dprl_cond(ICPR,k,j-1,i  ),
                                  dprl_cond(ICPR,k,j-1,i  ) - dprl_cond(ICPR,k,j-1,i-1));
        dplrho_dx /= dx;
        Real dpprho_dx = limiter4(dprl_cond(ICPP,k,j  ,i+1) - dprl_cond(ICPP,k,j  ,i),
                                  dprl_cond(ICPP,k,j  ,i  ) - dprl_cond(ICPP,k,j  ,i-1),
                                  dprl_cond(ICPP,k,j-1,i+1) - dprl_cond(ICPP,k,j-1,i),
                                  dprl_cond(ICPP,k,j-1,i  ) - dprl_cond(ICPP,k,j-1,i-1));
        dpprho_dx /= dx;
        Real dbmag_dx = limiter4( bmagcc_(k,j  ,i+1) - bmagcc_(k,j  ,i),
                                 bmagcc_(k,j  ,i  ) - bmagcc_(k,j  ,i-1),
                                 bmagcc_(k,j-1,i+1) - bmagcc_(k,j-1,i),
                                 bmagcc_(k,j-1,i  ) - bmagcc_(k,j-1,i-1));
        dbmag_dx /= dx;
        
        // x gradients at i-1/2 j k
        Real dplrho_dy = (dprl_cond(ICPR,k,j,i) - dprl_cond(ICPR,k,j-1,i) ) / dy;
        Real dpprho_dy = (dprl_cond(ICPP,k,j,i) - dprl_cond(ICPP,k,j-1,i) ) / dy;
        Real dbmag_dy = (bmagcc_(k,j,i) - bmagcc_(k,j-1,i) ) / dy;
        
        
        Real qprp = -2.*csprl2/(sqrt_twopi_*csprl*kl_lf + nu_c) *
                (rho*(bhx*dpprho_dx + bhy*dpprho_dy ) -
                       pprp*(1.-pprp/pprl)*
                            (bhx*dbmag_dx + bhy*dbmag_dy)/bmag);
        Real qprl = -8.*csprl2/(2.*sqrt_twopi_*csprl*kl_lf + threepi_m_eight_*nu_c) *
        rho * (bhx*dplrho_dx + bhy*dplrho_dy );
        
        x2flux(FLXEN,k,j,i) += bhy*(qprp + qprl/2.);
        x2flux(FLXMU,k,j,i) += bhy * qprp / bmag;
        
      }}
  } // Zero for 1D
  
  
  // ----------------------------------------------------------------- //
  // k-direction
  il=is-1, iu=ie+1, jl=js-1, ju=je+1;
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          // Compute bhat at the i-1/2 interface
          Real bx = 0.5*(bcc(IB1,k,j,i) + bcc(IB1,k-1,j,i));
          Real by = 0.5*(bcc(IB2,k,j,i) + bcc(IB2,k-1,j,i));
          Real bz = b.x3f(k,j,i);
          Real bmag = std::sqrt(bx*bx + by*by + bz*bz);
          bmag = (bmag > bfloor) ?  bmag : bfloor;
          Real bhx = bx/bmag;
          Real bhy = by/bmag;
          Real bhz = bz/bmag;
          
          // rho, pprl and prp to multiply derivatives
          Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k-1,j,i));
          Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k-1,j,i));
          Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k-1,j,i));
          Real csprl2 = 0.5*(dprl_cond(ICPR,k,j,i) + dprl_cond(ICPR,k-1,j,i));
          Real csprl = std::sqrt(csprl2);
          
          Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
          Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
          Real dz = pco_->dx3v(k-1);
          
          // x gradients at i j k-1/2
          Real dplrho_dx = limiter4(dprl_cond(ICPR,k  ,j,i+1) - dprl_cond(ICPR,k  ,j,i  ),
                                    dprl_cond(ICPR,k  ,j,i  ) - dprl_cond(ICPR,k  ,j,i-1),
                                    dprl_cond(ICPR,k-1,j,i+1) - dprl_cond(ICPR,k-1,j,i  ),
                                    dprl_cond(ICPR,k-1,j,i  ) - dprl_cond(ICPR,k-1,j,i-1));
          dplrho_dx /= dx;
          Real dpprho_dx = limiter4(dprl_cond(ICPP,k  ,j,i+1) - dprl_cond(ICPP,k  ,j,i  ),
                                    dprl_cond(ICPP,k  ,j,i  ) - dprl_cond(ICPP,k  ,j,i-1),
                                    dprl_cond(ICPP,k-1,j,i+1) - dprl_cond(ICPP,k-1,j,i  ),
                                    dprl_cond(ICPP,k-1,j,i  ) - dprl_cond(ICPP,k-1,j,i-1));
          dpprho_dx /= dx;
          Real dbmag_dx = limiter4(bmagcc_(k  ,j,i+1) - bmagcc_(k  ,j,i  ),
                                   bmagcc_(k  ,j,i  ) - bmagcc_(k  ,j,i-1),
                                   bmagcc_(k-1,j,i+1) - bmagcc_(k-1,j,i  ),
                                   bmagcc_(k-1,j,i  ) - bmagcc_(k-1,j,i-1));
          dbmag_dx /= dx;
          
          
          // y gradients at i j k-1/2
          Real dplrho_dy = limiter4(dprl_cond(ICPR,k  ,j+1,i) - dprl_cond(ICPR,k  ,j  ,i),
                                    dprl_cond(ICPR,k  ,j  ,i) - dprl_cond(ICPR,k  ,j-1,i),
                                    dprl_cond(ICPR,k-1,j+1,i) - dprl_cond(ICPR,k-1,j  ,i),
                                    dprl_cond(ICPR,k-1,j  ,i) - dprl_cond(ICPR,k-1,j-1,i));
          dplrho_dy /= dy;
          Real dpprho_dy =limiter4(dprl_cond(ICPP,k  ,j+1,i) - dprl_cond(ICPP,k  ,j  ,i),
                                   dprl_cond(ICPP,k  ,j  ,i) - dprl_cond(ICPP,k  ,j-1,i),
                                   dprl_cond(ICPP,k-1,j+1,i) - dprl_cond(ICPP,k-1,j  ,i),
                                   dprl_cond(ICPP,k-1,j  ,i) - dprl_cond(ICPP,k-1,j-1,i));
          dpprho_dy /= dy;
          Real dbmag_dy = limiter4(bmagcc_(k  ,j+1,i) - bmagcc_(k  ,j  ,i),
                                   bmagcc_(k  ,j  ,i) - bmagcc_(k  ,j-1,i),
                                   bmagcc_(k-1,j+1,i) - bmagcc_(k-1,j  ,i),
                                   bmagcc_(k-1,j  ,i) - bmagcc_(k-1,j-1,i));
          dbmag_dy /= dy;
          
          // x gradients at i-1/2 j k
          Real dplrho_dz = (dprl_cond(ICPR,k,j,i) - dprl_cond(ICPR,k-1,j,i) ) / dz;
          Real dpprho_dz = (dprl_cond(ICPP,k,j,i) - dprl_cond(ICPP,k-1,j,i) ) / dz;
          Real dbmag_dz = (bmagcc_(k,j,i) - bmagcc_(k-1,j,i) ) / dz;
          
          
          Real qprp = -2.*csprl2/(sqrt_twopi_*csprl*kl_lf + nu_c) *
                        (rho*(bhx*dpprho_dx + bhy*dpprho_dy + bhz*dpprho_dz) -
                            pprp*(1.-pprp/pprl)*
                                  (bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz)/bmag);
          Real qprl = -8.*csprl2/(2.*sqrt_twopi_*csprl*kl_lf + threepi_m_eight_*nu_c) *
          rho * (bhx*dplrho_dx + bhy*dplrho_dy + bhz*dplrho_dz);
          
          x3flux(FLXEN,k,j,i) += bhz*(qprp + qprl/2.);
          x3flux(FLXMU,k,j,i) += bhz * qprp / bmag;
          
        }}}
  } // Zero for 2D or 1D
  
  return;
}



//----------------------------------------------------------------------------------------
// constant conduction

void ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, int is, int ie, int js, int je, int ks, int ke) {
  if (phdif->kappa_iso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ISO,k,j,i) = phdif->kappa_iso;
      }
    }
  }
  if (phdif->kappa_aniso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ANI,k,j,i) = phdif->kappa_aniso;
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
// limiter2 and limiter4: call slope limiters to preserve monotonicity
#pragma omp declare simd simdlen(SIMD_WIDTH) notinbranch  // JS note: I think this declares that the function should be able to be called from vectorized loop, should test
inline static Real limiter2(const Real A, const Real B)
{
  // van Leer slope limiter
  return vanleer(A,B);
  
  // monotonized central (MC) limiter
  //  return minmod(2.0*minmod(A,B),0.5*(A+B));
}

#pragma omp declare simd simdlen(SIMD_WIDTH) notinbranch
inline static Real limiter4(const Real A, const Real B, const Real C, const Real D)
{
  return limiter2(limiter2(A,B),limiter2(C,D));
}

//----------------------------------------------------------------------------------------
// vanleer: van Leer slope limiter
#pragma omp declare simd simdlen(SIMD_WIDTH) notinbranch
inline static Real vanleer(const Real A, const Real B)
{
  if (A*B > 0) {
    return 2.0*A*B/(A+B);
  } else {
    return 0.0;
  }
//  if (A*B > 0) {
//    if (A > 0) {
//      return MIN(A,B);
//    } else {
//      return MAX(A,B);
//    }
//  } else {
//    return 0.0;
//  }
}
