//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turbulence.cpp
//  \brief implementation of functions in class Turbulence

// C/C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>
#include <algorithm>

// Athena++ headers
#include "athena_fft.hpp"
#include "turbulence.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../utils/utils.hpp"

//----------------------------------------------------------------------------------------
//! \fn TurbulenceDriver::TurbulenceDriver(Mesh *pm, ParameterInput *pin)
//  \brief TurbulenceDriver constructor

TurbulenceDriver::TurbulenceDriver(Mesh *pm, ParameterInput *pin)
 : FFTDriver(pm, pin) {

  rseed = pin->GetOrAddInteger("problem","rseed",-1); // seed for random number.

  nlow = pin->GetOrAddInteger("problem","nlow",0); // cut-off wavenumber
  // cut-off wavenumber, high:
  nhigh = pin->GetOrAddInteger("problem","nhigh",pm->mesh_size.nx1/2);
  expo = pin->GetOrAddReal("problem","expo",2); // power-law exponent
  dedt = pin->GetReal("problem","dedt"); // turbulence amplitude
  dtdrive = pin->GetOrAddReal("problem","dtdrive",0.0); // driving interval
  // Note that currently this does not allow time-correlated driving!
  tdrive = pm->time;
   
  // If set to i=1,2,3, no perturbations in dv_i direction, 0 is isotropic
  no_energy_in_i_direction = pin->GetOrAddInteger("problem","no_energy_in_i_direction",0);
   
  // Time-correlated forcing, has to create 3 extra arrays.
  tcorr = pin->GetOrAddReal("problem","tcorr",0.0); // Correlation time
  time_correlated = 0;
  if ((tcorr != 0.0) && (pm->turb_flag==3)) time_correlated=1;

  if (pm->turb_flag == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "Turbulence flag is set to zero! Shouldn't reach here!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  } else {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
#endif
  }

  int nx1=pm->pblock->block_size.nx1+2*NGHOST;
  int nx2=pm->pblock->block_size.nx2+2*NGHOST;
  int nx3=pm->pblock->block_size.nx3+2*NGHOST;

  vel = new AthenaArray<Real>[3];
  for (int nv=0; nv<3; nv++) vel[nv].NewAthenaArray(nmb,nx3,nx2,nx1);
   
  InitializeFFTBlock(true);
  QuickCreatePlan();
  dvol = pmy_fb->dx1*pmy_fb->dx2*pmy_fb->dx3;
   
  fv_curr_ = new AthenaFFTComplex*[3];
  for (int nv=0; nv<3; nv++) fv_curr_[nv] = new AthenaFFTComplex[pmy_fb->cnt_];
   
  if (time_correlated==1) {
    fv_old_ = new AthenaFFTComplex*[3];
    for (int nv=0; nv<3; nv++) fv_old_[nv] = new AthenaFFTComplex[pmy_fb->cnt_];
    InitializeOU();
    if (Globals::my_rank==0)
      std::cout << "Using time-correlated driving with tcorr = " << tcorr <<"\n";
   }

}

// destructor
TurbulenceDriver::~TurbulenceDriver() {
  for (int nv=0; nv<3; nv++) vel[nv].DeleteAthenaArray();
  delete [] vel;
  for (int nv=0; nv<3; nv++) delete[] fv_curr_[nv];
  delete [] fv_curr_;
  if (time_correlated==1) {
    for (int nv=0; nv<3; nv++) delete[] fv_old_[nv];
    delete [] fv_old_;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Driving(void)
//  \brief Generate and Perturb the velocity field

void TurbulenceDriver::Driving(void) {
  Mesh *pm=pmy_mesh_;
  bool new_perturb = false;

  // check driving time interval to generate new perturbation
  if (pm->time >= tdrive) {
    if (Globals::my_rank==0){
      std::cout << "generating turbulence at " << pm->time << std::endl;
    }
    Generate(pm->dt);
    new_perturb = true;
  }
  
  switch(pm->turb_flag) {
    case 1: // turb_flag == 1 : decaying turbulence
      Perturb(0);
      break;
    case 2: // turb_flag == 2 : impulsively driven turbulence
      if (new_perturb) {
        Perturb(dtdrive);
        tdrive = pm->time + dtdrive;
      }
      break;
    case 3: // turb_flag == 3 : continuously driven turbulence
      if (new_perturb){
        Perturb(pm->dt);
        tdrive = pm->time + pm->dt;
      }
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in TurbulenceDriver::Driving" << std::endl
      << "Turbulence flag " << pm->turb_flag << " is not supported!" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  
  return;

}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Generate()
//  \brief Generate velocity pertubation.

void TurbulenceDriver::Generate(Real dt) {
  Mesh *pm=pmy_mesh_;
  FFTBlock *pfb = pmy_fb;
  AthenaFFTPlan *plan = pfb->bplan_;
  AthenaFFTIndex *idx = pfb->b_in_;

  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
  int knx1=pfb->knx[0],knx2=pfb->knx[1],knx3=pfb->knx[2];
  
  // Form power spectrum for vx,vy,vz
  for (int nv = 0; nv<3; nv++){
    AthenaFFTComplex *fv = pfb->in_;
    PowerSpectrum(fv);
    for (int k=0; k<knx3; k++) {
      for (int j=0; j<knx2; j++) {
        for (int i=0; i<knx1; i++) {
          int64_t kidx=pfb->GetIndex(i,j,k,idx);
          if (no_energy_in_i_direction != nv+1) {
            // Copy to fv_curr_
            fv_curr_[nv][kidx][0] = fv[kidx][0];
            fv_curr_[nv][kidx][1] = fv[kidx][1];
          } else {
            fv_curr_[nv][kidx][0] = 0.;
            fv_curr_[nv][kidx][1] = 0.;
          }
    }}}
  }
  
  // Remove the divergence
  Project(fv_curr_);
  

  // If time correlated driving, add onto fv_old_ and copy to fv_curr_
  // Forcing still normalized in the same way
  if (time_correlated==1) OUEvolve(dt);
  
  
  // For each of fv_curr_, copy back into fv, take inverse fft, and copy into dv
  for (int nv = 0; nv<3; nv++){
    AthenaArray<Real> &dv = vel[nv], dv_mb;
    AthenaFFTComplex *fv = pfb->in_;
    
    // Copy into fv and take fft
    for (int k=0; k<knx3; k++) {
      for (int j=0; j<knx2; j++) {
        for (int i=0; i<knx1; i++) {
          int64_t kidx=pfb->GetIndex(i,j,k,idx);
          fv[kidx][0] = fv_curr_[nv][kidx][0];
          fv[kidx][1] = fv_curr_[nv][kidx][1];
    }}}
    pfb->Execute(plan);

    for (int igid=nbs, nb=0;igid<=nbe;igid++, nb++) {
      MeshBlock *pmb=pm->FindMeshBlock(igid);
      if (pmb != NULL) {
        dv_mb.InitWithShallowSlice(dv, 4, nb, 1); //(dv, 4, nb, 1)
        pfb->RetrieveResult(dv_mb,1,NGHOST,pmb->loc,pmb->block_size);
      }
//      std::cout << "nv=" << nv << ",nb="<<nb<<"\n";
      
    }
  }
  
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::PowerSpectrum(AthenaFFTComplex *amp)
//  \brief Generate Power spectrum in Fourier space with power-law

void TurbulenceDriver::PowerSpectrum(AthenaFFTComplex *amp) {
  Real pcoeff;
  FFTBlock *pfb = pmy_fb;
  AthenaFFTIndex *idx = pfb->b_in_;
  int knx1=pfb->knx[0],knx2=pfb->knx[1],knx3=pfb->knx[2];
// set random amplitudes with gaussian deviation
  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        Real q1=ran2(&rseed);
        Real q2=ran2(&rseed);
        Real q3=std::sqrt(-2.0*std::log(q1+1.e-20))*std::cos(2.0*PI*q2);
        q1=ran2(&rseed);
        int64_t kidx=pfb->GetIndex(i,j,k,idx);
        amp[kidx][0] = q3*std::cos(2.0*PI*q1);
        amp[kidx][1] = q3*std::sin(2.0*PI*q1);
      }
    }
  }

// set power spectrum: only power-law
  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        int64_t nx=GetKcomp(i,pfb->kdisp[0],pfb->kNx[0]);
        int64_t ny=GetKcomp(j,pfb->kdisp[1],pfb->kNx[1]);
        int64_t nz=GetKcomp(k,pfb->kdisp[2],pfb->kNx[2]);
        Real nmag = std::sqrt(nx*nx+ny*ny+nz*nz);
        Real kx=nx*pfb->dkx[0];
        Real ky=ny*pfb->dkx[1];
        Real kz=nz*pfb->dkx[2];
        Real kmag = std::sqrt(kx*kx+ky*ky+kz*kz);

        int64_t gidx = pfb->GetGlobalIndex(i,j,k);

        if (gidx == 0) {
          pcoeff = 0.0;
        } else {
          if ((abs(nx) > nlow) && (abs(nx) < nhigh) &&
              (abs(ny) > nlow) && (abs(ny) < nhigh) &&
              (abs(nz) > nlow) && (abs(nz) < nhigh)) {
            pcoeff = 1.0/std::pow(kmag,(expo+2.0)/2.0);
          } else {
            pcoeff = 0.0;
          }
        }
        int64_t kidx=pfb->GetIndex(i,j,k,idx);
        
        amp[kidx][0] *= pcoeff;
        amp[kidx][1] *= pcoeff;
//        if (nx==3 && ny==1 && nz==0){
//          amp[kidx][0] = 1.;
//          amp[kidx][1] = 0.;
//        } else {
//          amp[kidx][0] = 0.;
//          amp[kidx][1] = 0.;
//        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Project(AthenaFFTComplex **fv)
//  \brief Makes perturbation divergence free

void TurbulenceDriver::Project(AthenaFFTComplex **fv) {
  FFTBlock *pfb = pmy_fb;
  AthenaFFTIndex *idx = pfb->b_in_;
  int knx1=pfb->knx[0],knx2=pfb->knx[1],knx3=pfb->knx[2];
  Real divsum = 0.0;
  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        // Get khat
        int64_t nx=GetKcomp(i,pfb->kdisp[0],pfb->kNx[0]);
        int64_t ny=GetKcomp(j,pfb->kdisp[1],pfb->kNx[1]);
        int64_t nz=GetKcomp(k,pfb->kdisp[2],pfb->kNx[2]);
        Real kx=nx*pfb->dkx[0]; // dkx[i] gets weirdly permuted (based on meshblocks) when using MPI!!!
        Real ky=ny*pfb->dkx[1]; // JS: I *think* I fixed, Lx[i] wasn't in the functions to swap and permute axes
        Real kz=nz*pfb->dkx[2];
        Real kmag = std::sqrt(kx*kx + ky*ky + kz*kz);
        // khat stored in kx,ky,kz
        kx /= kmag;
        ky /= kmag;
        kz /= kmag;
        
        int64_t kidx=pfb->GetIndex(i,j,k,idx);
        if (kmag != 0.0) {
          // Form (khat.f)
          Real kdotf_re = kx*fv[0][kidx][0] + ky*fv[1][kidx][0] + kz*fv[2][kidx][0];
          Real kdotf_im = kx*fv[0][kidx][1] + ky*fv[1][kidx][1] + kz*fv[2][kidx][1];
          
          fv[0][kidx][0] -= kdotf_re * kx;
          fv[1][kidx][0] -= kdotf_re * ky;
          fv[2][kidx][0] -= kdotf_re * kz;
          
          fv[0][kidx][1] -= kdotf_im * kx;
          fv[1][kidx][1] -= kdotf_im * ky;
          fv[2][kidx][1] -= kdotf_im * kz;
          
//          if (fv[1][kidx][0] != 0.){
//            std::cout << "n=" << nx << "," << ny << "," << nz << "\n";
//            std::cout << "k=" << kx << "," << ky << "," << kz << "\n";
//            std::cout << "dk=" << pfb->dkx[0] << "," << pfb->dkx[1] << "," << pfb->dkx[2] << "\n";
//            std::cout << "Nx=" << pfb->b_in_->Nx[0] << "," << pfb->b_in_->Nx[1] << "," << pfb->b_in_->nx[2] << "\n";
//            std::cout << "Lx=" << pfb->b_in_->Lx[0] << "," << pfb->b_in_->Lx[1] << "," << pfb->b_in_->Lx[2] << "\n";
//            std::cout << "kdot = " << kdotf_re << "," << kdotf_im <<"\n";
//            std::cout << "proc = " << Globals::my_rank << "\n";
//          }
          
        }
  }}}
//  std::cout <<"\n";
//  std::cout << divsum <<"\n";
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Perturb(Real dt)
//  \brief Add velocity perturbation to the hydro variables

void TurbulenceDriver::Perturb(Real dt) {
  Mesh *pm = pmy_mesh_;
  std::stringstream msg;
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;

  int is=pm->pblock->is, ie=pm->pblock->ie;
  int js=pm->pblock->js, je=pm->pblock->je;
  int ks=pm->pblock->ks, ke=pm->pblock->ke;

  int mpierr;
  Real aa, b, c, s, de, v1, v2, v3, den, M1, M2, M3;
  Real m[4] = {0}, gm[4];
  AthenaArray<Real> &dv1 = vel[0], &dv2 = vel[1], &dv3 = vel[2];
  

  for (int igid=nbs, nb=0; igid<=nbe; igid++, nb++) {
    MeshBlock *pmb=pm->FindMeshBlock(igid);
    if (pmb != NULL) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
//            // Set perturbation in i direction to zero
//            dv1(nb,k,j,i) *= dv1_mult;
//            dv2(nb,k,j,i) *= dv2_mult;
//            dv3(nb,k,j,i) *= dv3_mult;
            // Sum momentum of perturbations
            den = pmb->phydro->u(IDN,k,j,i);
            m[0] += den;
            m[1] += den*dv1(nb,k,j,i);
            m[2] += den*dv2(nb,k,j,i);
            m[3] += den*dv3(nb,k,j,i);
          }
        }
      }
    }
  }

#ifdef MPI_PARALLEL
// Sum the perturbations over all processors
  mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    msg << "[normalize]: MPI_Allreduce error = "
        << mpierr << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  // Ask Changgoo about this
  for (int n=0; n<4; n++) m[n]=gm[n];
#endif // MPI_PARALLEL

  for (int nb=0; nb<nmb; nb++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          dv1(nb,k,j,i) -= m[1]/m[0];
          dv2(nb,k,j,i) -= m[2]/m[0];
          dv3(nb,k,j,i) -= m[3]/m[0];
        }
      }
    }
  }

  // Calculate unscaled energy of perturbations
  m[0] = 0.0;
  m[1] = 0.0;
  for (int igid=nbs, nb=0;igid<=nbe;igid++, nb++) {
    MeshBlock *pmb=pm->FindMeshBlock(igid);
    if (pmb != NULL) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            v1 = dv1(nb,k,j,i);
            v2 = dv2(nb,k,j,i);
            v3 = dv3(nb,k,j,i);
            den = pmb->phydro->u(IDN,k,j,i);
            M1 = pmb->phydro->u(IM1,k,j,i);
            M2 = pmb->phydro->u(IM2,k,j,i);
            M3 = pmb->phydro->u(IM3,k,j,i);
            m[0] += den*(SQR(v1) + SQR(v2) + SQR(v3));
            m[1] += M1*v1 + M2*v2 + M3*v3;
          }
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  // Sum the perturbations over all processors
  mpierr = MPI_Allreduce(m, gm, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    msg << "[normalize]: MPI_Allreduce error = "
        << mpierr << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  //  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  m[0] = gm[0];  m[1] = gm[1];
#endif // MPI_PARALLEL

  // Rescale to give the correct energy injection rate
  if (pm->turb_flag > 1) {
    // driven turbulence
    de = dedt*dt;
    if (Globals::my_rank==0)
      std::cout << "driven turbulence with de = " << de << std::endl;
  } else {
    // decaying turbulence (all in one shot)
    de = dedt;
    if (Globals::my_rank==0)
      std::cout << "decaying turbulence with de = " << de << std::endl;
  }
  aa = 0.5*m[0];
  aa = std::max(aa,static_cast<Real>(1.0e-20));
  b = m[1];
  c = -de/dvol;
  if (b >= 0.0)
    s = (-2.0*c)/(b + std::sqrt(b*b - 4.0*aa*c));
  else
    s = (-b + std::sqrt(b*b - 4.0*aa*c))/(2.0*aa);

  if (std::isnan(s)) std::cout << "[perturb]: s is NaN!" << std::endl;

  // Apply momentum pertubations
  for (int igid=nbs, nb=0; igid<=nbe; igid++, nb++) {
    MeshBlock *pmb=pm->FindMeshBlock(igid);
    if (pmb != NULL) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            v1 = dv1(nb,k,j,i);
            v2 = dv2(nb,k,j,i);
            v3 = dv3(nb,k,j,i);
            den = pmb->phydro->u(IDN,k,j,i);
            M1 = pmb->phydro->u(IM1,k,j,i);
            M2 = pmb->phydro->u(IM2,k,j,i);
            M3 = pmb->phydro->u(IM3,k,j,i);

            if (NON_BAROTROPIC_EOS) {
              pmb->phydro->u(IEN,k,j,i) += s*(M1*v1+M2*v2+M3*v3)
                                         + 0.5*s*s*den*(SQR(v1)+SQR(v2)+SQR(v3));
            }
            if (CGL_EOS) {
              pmb->phydro->u(IEN,k,j,i) += s*(M1*v1+M2*v2+M3*v3)
                                    + 0.5*s*s*den*(SQR(v1)+SQR(v2)+SQR(v3));
            }
            pmb->phydro->u(IM1,k,j,i) += s*den*v1;
            pmb->phydro->u(IM2,k,j,i) += s*den*v2;
            pmb->phydro->u(IM3,k,j,i) += s*den*v3;
          }
        }
      }
    }
  }
  return;

}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::GetKcomp(int idx, int disp, int Nx)
//  \brief Get k index, which runs from 0, 1, ... Nx/2-1, -Nx/2, -Nx/2+1, ..., -1.

int64_t TurbulenceDriver::GetKcomp(int idx, int disp, int Nx) {
  return ((idx+disp) - static_cast<int64_t>(2*(idx+disp)/Nx)*Nx);
}


//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::OUEvolve(AthenaFFTComplex **fv_curr_, AthenaFFTComplex **fv_old_, Real dt)
//  \brief evolves fv_old_ by dt with correlation time tcorr, and copies result in fv_curr_ to use in forcing.

void TurbulenceDriver::OUEvolve(Real dt){
  FFTBlock *pfb = pmy_fb;
  AthenaFFTIndex *idx = pfb->b_in_;
  int knx1=pfb->knx[0],knx2=pfb->knx[1],knx3=pfb->knx[2];
  Real dt_tcor  = dt/tcorr, sqrt_dt = std::sqrt(dt);
  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        int64_t kidx=pfb->GetIndex(i,j,k,idx);
        for (int reim=0; reim<2; reim ++){
          fv_old_[0][kidx][reim] -= dt_tcor * fv_old_[0][kidx][reim] + sqrt_dt * fv_curr_[0][kidx][reim];
          fv_old_[1][kidx][reim] -= dt_tcor * fv_old_[1][kidx][reim] + sqrt_dt * fv_curr_[1][kidx][reim];
          fv_old_[2][kidx][reim] -= dt_tcor * fv_old_[2][kidx][reim] + sqrt_dt * fv_curr_[2][kidx][reim];
          // Keep result in fv_old_ for use in next time step
          // Then also copy to fv_curr_
          fv_curr_[0][kidx][reim] = fv_old_[0][kidx][reim];
          fv_curr_[1][kidx][reim] = fv_old_[1][kidx][reim];
          fv_curr_[2][kidx][reim] = fv_old_[2][kidx][reim];

        }
  }}}
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::InitializeOU()
//  \brief Sets fv_old_ to it's steady state amplitude, so that you don't get a transient period of forcing at the start. The forcing is normalized after evolving the OU process, so its amplitude at steady state is just ~sqrt(tcorr/2)
// Also deals with restarts.
void TurbulenceDriver::InitializeOU(void){
  FFTBlock *pfb = pmy_fb;
  AthenaFFTIndex *idx = pfb->b_in_;
  int knx1=pfb->knx[0],knx2=pfb->knx[1],knx3=pfb->knx[2];
  Real sqrt_tc = std::sqrt(tcorr/2);

  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        int64_t kidx=pfb->GetIndex(i,j,k,idx);
        for (int reim=0; reim<2; reim ++){
          fv_old_[0][kidx][reim] = 0.0;
          fv_old_[1][kidx][reim] = 0.0;
          fv_old_[2][kidx][reim] = 0.0;
        }
  }}}
  // With fv_old_=0.0 this sets both fv_curr_ and fv_old_ to the basic generated noise
  Generate(1.0);

  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        int64_t kidx=pfb->GetIndex(i,j,k,idx);
        for (int reim=0; reim<2; reim ++){
          fv_old_[0][kidx][reim] = sqrt_tc * fv_curr_[0][kidx][reim];
          fv_old_[1][kidx][reim] = sqrt_tc * fv_curr_[1][kidx][reim];
          fv_old_[2][kidx][reim] = sqrt_tc * fv_curr_[2][kidx][reim];
        }
  }}}
}
