//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cgl.cpp
//  \brief implements functions in class EquationOfState for CGL model
//  Added J Squire 2018

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;
//  gamma_ = pin->GetReal("hydro", "gamma");
  density_floor_  = pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*(FLT_MIN)));
  magnetic_mag_floor_ = pin->GetOrAddReal("hydro", "bmagfloor", std::sqrt(1024*(FLT_MIN)));
  fh_hlld_floor_ = pin->GetOrAddReal("hydro", "fh_hlld_floor", 0.1);
  // Collisions and limiters
  collision_freq_ = pin->GetOrAddReal("problem", "nu_coll", 0.0);
  firehose_limiter_ = pin->GetOrAddInteger("problem", "firehose_limiter", 0);
  mirror_limiter_ = pin->GetOrAddInteger("problem", "mirror_limiter", 0);
  limiting_collision_freq_ = pin->GetOrAddReal("problem", "limiter_nu_coll", 1.e8);
}

// destructor

EquationOfState::~EquationOfState() {
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//    const AthenaArray<Real> &prim_old, const FaceField &b,
//    AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//    int il, int iu, int jl, int ju, int kl, int ku);
// \brief For the Hydro, converts conserved into primitive variables in CGL MHD.
//  For the Field, computes cell-centered from face-centered magnetic field.
//  Checks field is not too small

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
                                           const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
                                           AthenaArray<Real> &bcc, Coordinates *pco,
                                           int il, int iu, int jl, int ju, int kl, int ku) {
  
  pmy_block_->pfield->CalculateCellCenteredField(b,bcc,pco,il,iu,jl,ju,kl,ku);
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IVX,k,j,i);
        Real& u_m2 = cons(IVY,k,j,i);
        Real& u_m3 = cons(IVZ,k,j,i);
        Real& u_e  = cons(IEN,k,j,i);
        Real& u_mu = cons(IMU,k,j,i);
        
        Real& w_d  = prim(IDN,k,j,i);
        Real& w_vx = prim(IVX,k,j,i);
        Real& w_vy = prim(IVY,k,j,i);
        Real& w_vz = prim(IVZ,k,j,i);
        Real& w_p  = prim(IPR,k,j,i);
        Real& w_pp = prim(IPP,k,j,i);
        
        // apply density floor, without changing momentum or energy
        u_d = (u_d > density_floor_) ?  u_d : density_floor_;
        w_d = u_d;
        
        Real di = 1.0/u_d;
        w_vx = u_m1*di;
        w_vy = u_m2*di;
        w_vz = u_m3*di;
        
        const Real& bcc1 = bcc(IB1,k,j,i);
        const Real& bcc2 = bcc(IB2,k,j,i);
        const Real& bcc3 = bcc(IB3,k,j,i);
        
        Real pb = 0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
        Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
        Real bmag = std::sqrt(2.0*pb);
        if (bmag<magnetic_mag_floor_) {bmag = magnetic_mag_floor_; }; 
        // Might want to throw an error here or something...
        
        w_pp = bmag*u_mu;
        w_p = 2.0*(u_e - ke - pb - w_pp); //0.66666667*(u_e - ke - pb); //
        
        // apply pressure floor, correct total energy
        w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;
        w_pp = (w_pp > pressure_floor_) ?  w_pp : pressure_floor_;
        
        u_e = (w_p > pressure_floor_) ?  u_e : (1.5*pressure_floor_ + ke + pb);
        u_mu = (w_pp > pressure_floor_) ?  u_mu : (pressure_floor_/bmag);
      }
    }}
//  std::cout << std::endl;
  
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables
//        Note that this function assumes cell-centered fields are already calculated

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
                                           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
                                           int il, int iu, int jl, int ju, int kl, int ku) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);
        Real& u_e  = cons(IEN,k,j,i);
        Real& u_mu = cons(IMU,k,j,i);
        
        const Real& w_d  = prim(IDN,k,j,i);
        const Real& w_vx = prim(IVX,k,j,i);
        const Real& w_vy = prim(IVY,k,j,i);
        const Real& w_vz = prim(IVZ,k,j,i);
        const Real& w_p  = prim(IPR,k,j,i);
        const Real& w_pp = prim(IPP,k,j,i);
        
        const Real& bcc1 = bc(IB1,k,j,i);
        const Real& bcc2 = bc(IB2,k,j,i);
        const Real& bcc3 = bc(IB3,k,j,i);
        
        Real bsqr = SQR(bcc1) + SQR(bcc2) + SQR(bcc3);
        Real bmag = std::sqrt(bsqr);
        if (bmag<magnetic_mag_floor_) {bmag = magnetic_mag_floor_; };
        // Might want to throw an error here or something...
        
        u_d = w_d;
        u_m1 = w_vx*w_d;
        u_m2 = w_vy*w_d;
        u_m3 = w_vz*w_d;
//        u_e = w_p*1.5 + 0.5*(w_d*(SQR(w_vx) + SQR(w_vy) + SQR(w_vz))  + bsqr );
        u_e = 0.5*w_p + w_pp + 0.5*(w_d*(SQR(w_vx) + SQR(w_vy) + SQR(w_vz)) + bsqr);
        u_mu = w_pp/bmag;

      }
    }}
  
  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables
//  for CGL, this is taken as parallel sound speed p_parr/rho
Real EquationOfState::SoundSpeed(const Real prim[NHYDRO]) {
  std::cout << "Shouldn't be here!!" << std::endl;
  return std::sqrt(prim[IPR]/prim[IDN]);
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::FastMagnetosonicSpeed(const Real prim[], const Real bx)
// \brief returns fast magnetosonic speed given vector of primitive variables
// In CGL the fast and slow magnetosonic modes are different and complicated, although
// generally have a similar magnitude to stanard MHD (e.g., fast ~sqrt(v_A^2+cs^2). Expressions
//  can be found in e.g., Baranov 1970  10.1007/BF01080231. Note Meng et al. JCP 2012, equation
// (39) is incorrect, pprp^2(1-bx^2)/4 should be pprp^2(1-bx^2)*bx^2
Real EquationOfState::FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx) {
  Real bx2 = bx*bx;
  Real bprp2 = (prim[IBY]*prim[IBY] + prim[IBZ]*prim[IBZ]);
  Real bhatx2 = bx2/(bx2+bprp2);
  Real pprp = prim[IPP];
  Real pprl = prim[IPR];
  Real qsq = bx2 + bprp2 + 2*pprp + (2.*pprl - pprp)*bhatx2;
  return std::sqrt( 0.5*(qsq + std::sqrt(fabs(qsq*qsq + 4.*pprp*pprp*(1. - bhatx2)*bhatx2
                                  - 12.*pprl*pprp*bhatx2*(2.-bhatx2)
                                + 12.*pprl*pprl*bhatx2*bhatx2 - 12.*bx2*pprl)))/prim[IDN] );
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::Collisions(const Real prim[], const Real bx)
// \brief Operates on primitive variables to collisionally relax delta_p using implicit method.
//  Also enforces microinstability limiters. See Sharma et al. 2006 App. A3 for details.
void EquationOfState::Collisions(AthenaArray<Real> &prim, const AthenaArray<Real> &prim_old,
                                 const AthenaArray<Real> &bc,
                                 Real dt, int il, int iu, int jl, int ju, int kl, int ku) {
  // Compute basic collisions from solving dt(pprp,pprl) = 1/3*(-nu*Dp,2*nu*Dp) with ICs as current pprp,pprl
  if (collision_freq_ != 0.0) {
    Real expdtnu = std::exp(-collision_freq_*dt);
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        
        Real& w_pl = prim(IPR,k,j,i);
        Real& w_pp = prim(IPP,k,j,i);
        Real w_pp_tmp;
        
        w_pp_tmp = (ONE_3RD*expdtnu + TWO_3RD)*w_pp + (ONE_3RD - ONE_3RD*expdtnu)*w_pl;
        w_pl = (TWO_3RD - TWO_3RD*expdtnu)*w_pp + (TWO_3RD*expdtnu + ONE_3RD)*w_pl;
        w_pp = w_pp_tmp;
      }
    }}
  }
  
  // Mirror/firehose limiters.
  if (firehose_limiter_ || mirror_limiter_) {
    Real nudt = limiting_collision_freq_*dt;
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
  #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real& w_pl = prim(IPR,k,j,i);
          Real& w_pp = prim(IPP,k,j,i);
          Real w_pp_tmp;
          
          const Real& bcc1 = bc(IB1,k,j,i);
          const Real& bcc2 = bc(IB2,k,j,i);
          const Real& bcc3 = bc(IB3,k,j,i);
          Real bsqr = SQR(bcc1) + SQR(bcc2) + SQR(bcc3);
          
          if (firehose_limiter_ && (w_pp - w_pl < -bsqr)) {
            w_pp_tmp = (3.*w_pp + nudt*(2.*w_pp + w_pl - bsqr ))/(3.+3.*nudt);
            w_pl = (3.*w_pl + nudt*(2.*w_pp + w_pl + 2.*bsqr ))/(3.+3.*nudt);
            w_pp = w_pp_tmp;
          }
          
          if (mirror_limiter_ && (w_pp - w_pl > 0.5*bsqr)) {
            w_pp_tmp = (3.*w_pp + nudt*(2.*w_pp + w_pl + 0.5*bsqr ))/(3.+3.*nudt);
            w_pl = (3.*w_pl + nudt*(2.*w_pp + w_pl - bsqr ))/(3.+3.*nudt);
            w_pp = w_pp_tmp;
          }
        }
    }}
  }
  
  
  return;
  
}


//---------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim,
//           int k, int j, int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface states
void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i) {
  Real& w_d  = prim(IDN,k,j,i);
  Real& w_p  = prim(IPR,k,j,i);
  Real& w_pp  = prim(IPP,k,j,i);
  
  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floors
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;
  w_pp = (w_pp > pressure_floor_) ?  w_pp : pressure_floor_;
  
  return;
}
