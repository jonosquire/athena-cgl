//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C++ headers
#include <cmath>

// Athena++ headers
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

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
// heat fluxes in CGL, see Sharma et al. 2006
void HydroDiffusion::HeatFlux_LandauFluid(const AthenaArray<Real> &prim,const AthenaArray<Real> &cons,
                          AthenaArray<Real> *flx,
                          const FaceField &b, const AthenaArray<Real> &bcc){
  
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real bfloor = pmb_->peos->GetBFieldFloor();
  Real sqrt_twopi = std::sqrt(2.*PI);
  Real threepi_m_eight = 3.*PI - 8.;
  Real nu_c = pmb_->peos->GetCollisionFreq();
  
  // Store pprp/rho and pprl/rho to save lots of computation
  il=is-1, iu=ie+1, jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        pprl_rho_(k,j,i) = prim(IPR,k,j,i) / prim(IDN,k,j,i);
        pprp_rho_(k,j,i) = prim(IPP,k,j,i) / prim(IDN,k,j,i);
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
        
        // If you want to compute B derivatives along the i-1/2 line instead.
//        // B at j+1,j-1,k+1, and k-1 for taking d|B|/dy,z
//        by = 0.5*(bcc(IB2,k,j+1,i) + bcc(IB2,k,j+1,i-1));
//        bz = 0.5*(bcc(IB3,k,j+1,i) + bcc(IB3,k,j+1,i-1));
//        Real bmagjp = std::sqrt(b.x1f(k,j+1,i)*b.x1f(k,j+1,i) + by*by + bz*bz);
//        bmagjp = std::max(bmagjp, bfloor);
//
//        by = 0.5*(bcc(IB2,k,j-1,i) + bcc(IB2,k,j-1,i-1));
//        bz = 0.5*(bcc(IB3,k,j-1,i) + bcc(IB3,k,j-1,i-1));
//        Real bmagjm = std::sqrt(b.x1f(k,j-1,i)*b.x1f(k,j-1,i) + by*by + bz*bz);
//        bmagjm = std::max(bmagjm, bfloor);
//
//        by = 0.5*(bcc(IB2,k+1,j,i) + bcc(IB2,k+1,j,i-1));
//        bz = 0.5*(bcc(IB3,k+1,j,i) + bcc(IB3,k+1,j,i-1));
//        Real bmagkp = std::sqrt(b.x1f(k+1,j,i)*b.x1f(k+1,j,i) + by*by + bz*bz);
//        bmagkp = std::max(bmagkp, bfloor);
//
//        by = 0.5*(bcc(IB2,k-1,j,i) + bcc(IB2,k-1,j,i-1));
//        bz = 0.5*(bcc(IB3,k-1,j,i) + bcc(IB3,k-1,j,i-1));
//        Real bmagkm = std::sqrt(b.x1f(k-1,j,i)*b.x1f(k-1,j,i) + by*by + bz*bz);
//        bmagkm = std::max(bmagkm, bfloor);
//
//        // B at i-1,j,k and i,j,k to take x derivatives
//        bx = 0.5*(b.x1f(k,j,i) + b.x1f(k,j,i+1));
//        Real bmagi = std::sqrt(bx*bx + bcc(IB2,k,j,i)*bcc(IB2,k,j,i) +
//                               bcc(IB3,k,j,i)*bcc(IB3,k,j,i));
//        bmagi = std::max(bmagi, bfloor);
//        bx = 0.5*(b.x1f(k,j,i-1) + b.x1f(k,j,i));
//        Real bmagim = std::sqrt(bx*bx + bcc(IB2,k,j,i-1)*bcc(IB2,k,j,i-1) +
//                                 bcc(IB3,k,j,i-1)*bcc(IB3,k,j,i-1));
//        bmagim = std::max(bmagim, bfloor);
        
        // rho, pprl and prp to multiply derivatives
        Real rho = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
        Real pprl = 0.5*(prim(IPR,k,j,i) + prim(IPR,k,j,i-1));
        Real pprp = 0.5*(prim(IPP,k,j,i) + prim(IPP,k,j,i-1));
        Real csprl2 = 0.5*(pprl_rho_(k,j,i) + pprl_rho_(k,j,i-1));
        Real csprl = std::sqrt(csprl2);

        Real dx = pco_->dx1v(i-1);
        Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
        Real dz = 0.5*(pco_->dx3v(k-1) + pco_->dx3v(k));

        // x gradients at i-1/2 j k
        Real dplrho_dx = (pprl_rho_(k,j,i) - pprl_rho_(k,j,i-1) ) / dx;
        Real dpprho_dx = (pprp_rho_(k,j,i) - pprp_rho_(k,j,i-1) ) / dx;
        Real dbmag_dx = (bmagcc_(k,j,i) - bmagcc_(k,j,i-1) ) / dx;
        
        // y gradients at i-1/2 j k
        Real dplrho_dy = limiter4(pprl_rho_(k,j+1,i  ) - pprl_rho_(k,j  ,i  ),
                                  pprl_rho_(k,j  ,i  ) - pprl_rho_(k,j-1,i  ),
                                  pprl_rho_(k,j+1,i-1) - pprl_rho_(k,j  ,i-1),
                                  pprl_rho_(k,j  ,i-1) - pprl_rho_(k,j-1,i-1));
        dplrho_dy /= dy;
        Real dpprho_dy = limiter4(pprp_rho_(k,j+1,i  ) - pprp_rho_(k,j  ,i  ),
                                  pprp_rho_(k,j  ,i  ) - pprp_rho_(k,j-1,i  ),
                                  pprp_rho_(k,j+1,i-1) - pprp_rho_(k,j  ,i-1),
                                  pprp_rho_(k,j  ,i-1) - pprp_rho_(k,j-1,i-1));
        dpprho_dy /= dy;
        Real dbmag_dy = limiter4( bmagcc_(k,j+1,i  ) - bmagcc_(k,j  ,i  ),
                                  bmagcc_(k,j  ,i  ) - bmagcc_(k,j-1,i  ),
                                  bmagcc_(k,j+1,i-1) - bmagcc_(k,j  ,i-1),
                                  bmagcc_(k,j  ,i-1) - bmagcc_(k,j-1,i-1));
        dbmag_dy /= dy;
        
        // z gradients at i-1/2 j k
        Real dplrho_dz = limiter4(pprl_rho_(k+1,j,i  ) - pprl_rho_(k  ,j,i  ),
                                  pprl_rho_(k  ,j,i  ) - pprl_rho_(k-1,j,i  ),
                                  pprl_rho_(k+1,j,i-1) - pprl_rho_(k  ,j,i-1),
                                  pprl_rho_(k  ,j,i-1) - pprl_rho_(k-1,j,i-1));
        dplrho_dz /= dz;
        Real dpprho_dz = limiter4(pprp_rho_(k+1,j,i  ) - pprp_rho_(k  ,j,i  ),
                                  pprp_rho_(k  ,j,i  ) - pprp_rho_(k-1,j,i  ),
                                  pprp_rho_(k+1,j,i-1) - pprp_rho_(k  ,j,i-1),
                                  pprp_rho_(k  ,j,i-1) - pprp_rho_(k-1,j,i-1));
        dpprho_dz /= dz;
        Real dbmag_dz = limiter4( bmagcc_(k+1,j,i  ) - bmagcc_(k  ,j,i  ),
                                  bmagcc_(k  ,j,i  ) - bmagcc_(k-1,j,i  ),
                                  bmagcc_(k+1,j,i-1) - bmagcc_(k  ,j,i-1),
                                  bmagcc_(k  ,j,i-1) - bmagcc_(k-1,j,i-1));
        dbmag_dz /= dz;
        
        // Use computation in the i direction to compute kappa for time step later
        Real kappa_prp = 2.*csprl2/(sqrt_twopi*csprl*kl_lf + nu_c);
        Real kappa_prl = 8.*csprl2/(2.*sqrt_twopi*csprl*kl_lf + threepi_m_eight*nu_c);
        
        Real qprp = - kappa_prp * (rho*(bhx*dpprho_dx + bhy*dpprho_dy + bhz*dpprho_dz) -
                          - pprp*(1.-pprp/pprl)*
                              (bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz)/bmag);
        Real qprl = - kappa_prl * rho * (bhx*dplrho_dx + bhy*dplrho_dy + bhz*dplrho_dz);
        
        if (i != ie+1) kappa(k,j,i) = std::max(kappa_prp,kappa_prl);
        
        x1flux(FLXEN,k,j,i) += bhx*(qprp + qprl/2.);
        x1flux(FLXMU,k,j,i) += bhx * qprp / bmag;
        
  }}}

  // ----------------------------------------------------------------- //
  // j direction
  il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  for (int k=kl; k<=ku; ++k) {
    for (int j=js; j<=je+1; ++j) {
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
        Real csprl2 = 0.5*(pprl_rho_(k,j,i) + pprl_rho_(k,j-1,i));
        Real csprl = std::sqrt(csprl2);
        
        Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
        Real dy = pco_->dx2v(j-1);
        Real dz = 0.5*(pco_->dx3v(k-1) + pco_->dx3v(k));
        
        // x gradients at i j-1/2 k
        Real dplrho_dx = limiter4(pprl_rho_(k,j  ,i+1) - pprl_rho_(k,j  ,i  ),
                                  pprl_rho_(k,j  ,i  ) - pprl_rho_(k,j  ,i-1),
                                  pprl_rho_(k,j-1,i+1) - pprl_rho_(k,j-1,i  ),
                                  pprl_rho_(k,j-1,i  ) - pprl_rho_(k,j-1,i-1));
        dplrho_dx /= dx;
        Real dpprho_dx = limiter4(pprp_rho_(k,j  ,i+1) - pprp_rho_(k,j  ,i),
                                  pprp_rho_(k,j  ,i  ) - pprp_rho_(k,j  ,i-1),
                                  pprp_rho_(k,j-1,i+1) - pprp_rho_(k,j-1,i),
                                  pprp_rho_(k,j-1,i  ) - pprp_rho_(k,j-1,i-1));
        dpprho_dx /= dx;
        Real dbmag_dx = limiter4( bmagcc_(k,j  ,i+1) - bmagcc_(k,j  ,i),
                                  bmagcc_(k,j  ,i  ) - bmagcc_(k,j  ,i-1),
                                  bmagcc_(k,j-1,i+1) - bmagcc_(k,j-1,i),
                                  bmagcc_(k,j-1,i  ) - bmagcc_(k,j-1,i-1));
        dbmag_dx /= dx;
        
        // x gradients at i-1/2 j k
        Real dplrho_dy = (pprl_rho_(k,j,i) - pprl_rho_(k,j-1,i) ) / dy;
        Real dpprho_dy = (pprp_rho_(k,j,i) - pprp_rho_(k,j-1,i) ) / dy;
        Real dbmag_dy = (bmagcc_(k,j,i) - bmagcc_(k,j-1,i) ) / dy;
        
        // z gradients at i j-1/2 k
        Real dplrho_dz = limiter4(pprl_rho_(k+1,j  ,i) - pprl_rho_(k  ,j  ,i),
                                  pprl_rho_(k  ,j  ,i) - pprl_rho_(k-1,j  ,i),
                                  pprl_rho_(k+1,j-1,i) - pprl_rho_(k  ,j-1,i),
                                  pprl_rho_(k  ,j-1,i) - pprl_rho_(k-1,j-1,i));
        dplrho_dz /= dz;
        Real dpprho_dz = limiter4(pprp_rho_(k+1,j  ,i) - pprp_rho_(k  ,j  ,i),
                                  pprp_rho_(k  ,j  ,i) - pprp_rho_(k-1,j  ,i),
                                  pprp_rho_(k+1,j-1,i) - pprp_rho_(k  ,j-1,i),
                                  pprp_rho_(k  ,j-1,i) - pprp_rho_(k-1,j-1,i));
        dpprho_dz /= dz;
        Real dbmag_dz = limiter4(bmagcc_(k+1,j  ,i) - bmagcc_(k  ,j  ,i),
                                 bmagcc_(k  ,j  ,i) - bmagcc_(k-1,j  ,i),
                                 bmagcc_(k+1,j-1,i) - bmagcc_(k  ,j-1,i),
                                 bmagcc_(k  ,j-1,i) - bmagcc_(k-1,j-1,i));
        dbmag_dz /= dz;
        
        Real qprp = -2.*csprl2/(sqrt_twopi*csprl*kl_lf + nu_c) *
                      (rho*(bhx*dpprho_dx + bhy*dpprho_dy + bhz*dpprho_dz) -
                            - pprp*(1.-pprp/pprl)*
                        (bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz)/bmag);
        Real qprl = -8.*csprl2/(2.*sqrt_twopi*csprl*kl_lf + threepi_m_eight*nu_c) *
                      rho * (bhx*dplrho_dx + bhy*dplrho_dy + bhz*dplrho_dz);
        
        x2flux(FLXEN,k,j,i) += bhy*(qprp + qprl/2.);
        x2flux(FLXMU,k,j,i) += bhy * qprp / bmag;
        
  }}}
  
  
  // ----------------------------------------------------------------- //
  // k-direction
  il=is-1, iu=ie+1, jl=js-1, ju=je+1;
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
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
        Real csprl2 = 0.5*(pprl_rho_(k,j,i) + pprl_rho_(k-1,j,i));
        Real csprl = std::sqrt(csprl2);
        
        Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
        Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
        Real dz = pco_->dx3v(k-1);
        
        // x gradients at i j k-1/2
        Real dplrho_dx = limiter4(pprl_rho_(k  ,j,i+1) - pprl_rho_(k  ,j,i  ),
                                  pprl_rho_(k  ,j,i  ) - pprl_rho_(k  ,j,i-1),
                                  pprl_rho_(k-1,j,i+1) - pprl_rho_(k-1,j,i  ),
                                  pprl_rho_(k-1,j,i  ) - pprl_rho_(k-1,j,i-1));
        dplrho_dx /= dx;
        Real dpprho_dx = limiter4(pprp_rho_(k  ,j,i+1) - pprp_rho_(k  ,j,i  ),
                                  pprp_rho_(k  ,j,i  ) - pprp_rho_(k  ,j,i-1),
                                  pprp_rho_(k-1,j,i+1) - pprp_rho_(k-1,j,i  ),
                                  pprp_rho_(k-1,j,i  ) - pprp_rho_(k-1,j,i-1));
        dpprho_dx /= dx;
        Real dbmag_dx = limiter4(bmagcc_(k  ,j,i+1) - bmagcc_(k  ,j,i  ),
                                 bmagcc_(k  ,j,i  ) - bmagcc_(k  ,j,i-1),
                                 bmagcc_(k-1,j,i+1) - bmagcc_(k-1,j,i  ),
                                 bmagcc_(k-1,j,i  ) - bmagcc_(k-1,j,i-1));
        dbmag_dx /= dx;
        
        
        // y gradients at i j k-1/2
        Real dplrho_dy = limiter4(pprl_rho_(k  ,j+1,i) - pprl_rho_(k  ,j  ,i),
                                  pprl_rho_(k  ,j  ,i) - pprl_rho_(k  ,j-1,i),
                                  pprl_rho_(k-1,j+1,i) - pprl_rho_(k-1,j  ,i),
                                  pprl_rho_(k-1,j  ,i) - pprl_rho_(k-1,j-1,i));
        dplrho_dy /= dy;
        Real dpprho_dy =limiter4(pprp_rho_(k  ,j+1,i) - pprp_rho_(k  ,j  ,i),
                                 pprp_rho_(k  ,j  ,i) - pprp_rho_(k  ,j-1,i),
                                 pprp_rho_(k-1,j+1,i) - pprp_rho_(k-1,j  ,i),
                                 pprp_rho_(k-1,j  ,i) - pprp_rho_(k-1,j-1,i));
        dpprho_dy /= dy;
        Real dbmag_dy = limiter4(bmagcc_(k  ,j+1,i) - bmagcc_(k  ,j  ,i),
                                 bmagcc_(k  ,j  ,i) - bmagcc_(k  ,j-1,i),
                                 bmagcc_(k-1,j+1,i) - bmagcc_(k-1,j  ,i),
                                 bmagcc_(k-1,j  ,i) - bmagcc_(k-1,j-1,i));
        dbmag_dy /= dy;
        
        // x gradients at i-1/2 j k
        Real dplrho_dz = (pprl_rho_(k,j,i) - pprl_rho_(k-1,j,i) ) / dz;
        Real dpprho_dz = (pprp_rho_(k,j,i) - pprp_rho_(k-1,j,i) ) / dz;
        Real dbmag_dz = (bmagcc_(k,j,i) - bmagcc_(k-1,j,i) ) / dz;
        
        
        Real qprp = -2.*csprl2/(sqrt_twopi*csprl*kl_lf + nu_c) *
                      (rho*(bhx*dpprho_dx + bhy*dpprho_dy + bhz*dpprho_dz) -
                          - pprp*(1.-pprp/pprl)*
                            (bhx*dbmag_dx + bhy*dbmag_dy + bhz*dbmag_dz)/bmag);
        Real qprl = -8.*csprl2/(2.*sqrt_twopi*csprl*kl_lf + threepi_m_eight*nu_c) *
                        rho * (bhx*dplrho_dx + bhy*dplrho_dy + bhz*dplrho_dz);
        
        x3flux(FLXEN,k,j,i) += bhz*(qprp + qprl/2.);
        x3flux(FLXMU,k,j,i) += bhz * qprp / bmag;

  }}}
  
  
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
}
