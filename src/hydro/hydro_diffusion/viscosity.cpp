//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief functions to calculate viscous stresses

// Athena++ headers
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

// Declarations
static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);
static Real minmod  (const Real A, const Real B);


//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ViscousFlux_iso
//  \brief Calculate isotropic viscous stress as fluxes

void HydroDiffusion::ViscousFlux_iso(const AthenaArray<Real> &prim,
                                     const AthenaArray<Real> &cons,
                                     AthenaArray<Real> *visflx) {
  AthenaArray<Real> &x1flux=visflx[X1DIR];
  AthenaArray<Real> &x2flux=visflx[X2DIR];
  AthenaArray<Real> &x3flux=visflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real nu1, denf, flx1, flx2, flx3;
  Real nuiso2 = 0.0;//- TWO_3RD;

  Divv(prim, divv_);

  // Calculate the flux across each face.
  // i-direction
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb_->block_size.nx2 > 1) {
      if (pmb_->block_size.nx3 == 1) // 2D MHD limits
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D MHD limits
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      FaceXdx(k,j,is,ie+1,prim,fx_);
      FaceXdy(k,j,is,ie+1,prim,fy_);
      FaceXdz(k,j,is,ie+1,prim,fz_);
      for (int i=is; i<=ie+1; ++i) {
        nu1  = 0.5*(nu(ISO,k,j,i)   + nu(ISO,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
        flx1 = -denf*nu1*(fx_(i) + nuiso2*0.5*(divv_(k,j,i) + divv_(k,j,i-1)));
        flx2 = -denf*nu1*fy_(i);
        flx3 = -denf*nu1*fz_(i);
        x1flux(IM1,k,j,i) += flx1;
        x1flux(IM2,k,j,i) += flx2;
        x1flux(IM3,k,j,i) += flx3;
        if (NON_BAROTROPIC_EOS) // Could add CGL_EOS here.. Depends whether viscosity heats or not
          x1flux(IEN,k,j,i) += 0.5*((prim(IM1,k,j,i-1) + prim(IM1,k,j,i))*flx1 +
                                    (prim(IM2,k,j,i-1) + prim(IM2,k,j,i))*flx2 +
                                    (prim(IM3,k,j,i-1) + prim(IM3,k,j,i))*flx3);
      }
    }
  }

  // j-direction
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb_->block_size.nx3 == 1) // 2D MHD limits
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D MHD limits
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if (pmb_->block_size.nx2 > 1) { // modify x2flux for 2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
        // compute fluxes
        FaceYdx(k,j,is,ie,prim,fx_);
        FaceYdy(k,j,is,ie,prim,fy_);
        FaceYdz(k,j,is,ie,prim,fz_);
        // store fluxes
        for(int i=il; i<=iu; i++) {
          nu1  = 0.5*(nu(ISO,k,j,i)    + nu(ISO,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j-1,i)+ prim(IDN,k,j,i));
          flx1 = -denf*nu1*fx_(i);
          flx2 = -denf*nu1*(fy_(i) + nuiso2*0.5*(divv_(k,j-1,i) + divv_(k,j,i)));
          flx3 = -denf*nu1*fz_(i);
          x2flux(IM1,k,j,i) += flx1;
          x2flux(IM2,k,j,i) += flx2;
          x2flux(IM3,k,j,i) += flx3;
          if (NON_BAROTROPIC_EOS)// Could add CGL_EOS here.. Depends whether viscosity heats or not
            x2flux(IEN,k,j,i) += 0.5*((prim(IM1,k,j,i) + prim(IM1,k,j-1,i))*flx1 +
                                      (prim(IM2,k,j,i) + prim(IM2,k,j-1,i))*flx2 +
                                      (prim(IM3,k,j,i) + prim(IM3,k,j-1,i))*flx3);
        }
      }
    }
  } else { // modify x2flux for 1D
    // compute fluxes
    FaceYdx(ks,js,is,ie,prim,fx_);
    FaceYdy(ks,js,is,ie,prim,fy_);
    FaceYdz(ks,js,is,ie,prim,fz_);
    // store fluxes
    for(int i=il; i<=iu; i++) {
      nu1  = nu(ISO,ks,js,i);
      denf = prim(IDN,ks,js,i);
      flx1 = -denf*nu1*fx_(i);
      flx2 = -denf*nu1*(fy_(i) + nuiso2*divv_(ks,js,i));
      flx3 = -denf*nu1*fz_(i);
      x2flux(IM1,ks,js,i) += flx1;
      x2flux(IM2,ks,js,i) += flx2;
      x2flux(IM3,ks,js,i) += flx3;
      if (NON_BAROTROPIC_EOS)// Could add CGL_EOS here.. Depends whether viscosity heats or not
        x2flux(IEN,ks,js,i) += prim(IM1,ks,js,i)*flx1 +
                               prim(IM2,ks,js,i)*flx2 +
                               prim(IM3,ks,js,i)*flx3;
    }
    for(int i=il; i<=iu; i++) {
      x2flux(IM1,ks,je+1,i) = x2flux(IM1,ks,js,i);
      x2flux(IM2,ks,je+1,i) = x2flux(IM2,ks,js,i);
      x2flux(IM3,ks,je+1,i) = x2flux(IM3,ks,js,i);
      if (NON_BAROTROPIC_EOS)
        x2flux(IEN,ks,je+1,i) = x2flux(IEN,ks,js,i);
    }
  }
  // k-direction
  // set the loop limits
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb_->block_size.nx2 > 1) // 2D or 3D MHD limits
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D MHD limits
      il=is-1, iu=ie+1;
  }
  if (pmb_->block_size.nx3 > 1) { // modify x3flux for 3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // compute fluxes
        FaceZdx(k,j,is,ie,prim,fx_);
        FaceZdy(k,j,is,ie,prim,fy_);
        FaceZdz(k,j,is,ie,prim,fz_);
        // store fluxes
        for(int i=il; i<=iu; i++) {
          nu1  = 0.5*(nu(ISO,k,j,i)     + nu(ISO,k-1,j,i));
          denf = 0.5*(prim(IDN,k-1,j,i) + prim(IDN,k,j,i));
          flx1 = -denf*nu1*fx_(i);
          flx2 = -denf*nu1*fy_(i);
          flx3 = -denf*nu1*(fz_(i) + nuiso2*0.5*(divv_(k-1,j,i) + divv_(k,j,i)));
          x3flux(IM1,k,j,i) += flx1;
          x3flux(IM2,k,j,i) += flx2;
          x3flux(IM3,k,j,i) += flx3;
          if (NON_BAROTROPIC_EOS)
            x3flux(IEN,k,j,i) += 0.5*((prim(IM1,k,j,i) + prim(IM1,k-1,j,i))*flx1 +
                                      (prim(IM2,k,j,i) + prim(IM2,k-1,j,i))*flx2 +
                                      (prim(IM3,k,j,i) + prim(IM3,k-1,j,i))*flx3);
        }
      }
    }
  } else { // modify x2flux for 1D or 2D
    for (int j=jl; j<=ju; ++j) {
      // compute fluxes
      FaceZdx(ks,j,is,ie,prim,fx_);
      FaceZdy(ks,j,is,ie,prim,fy_);
      FaceZdz(ks,j,is,ie,prim,fz_);
      // store fluxes
      for(int i=il; i<=iu; i++) {
        nu1 = nu(ISO,ks,j,i);
        denf = prim(IDN,ks,j,i);
        flx1 = -denf*nu1*fx_(i);
        flx2 = -denf*nu1*fy_(i);
        flx3 = -denf*nu1*(fz_(i) + nuiso2*divv_(ks,j,i));
        x3flux(IM1,ks,j,i) += flx1;
        x3flux(IM2,ks,j,i) += flx2;
        x3flux(IM3,ks,j,i) += flx3;
        x3flux(IM1,ke+1,j,i) = x3flux(IM1,ks,j,i);
        x3flux(IM2,ke+1,j,i) = x3flux(IM2,ks,j,i);
        x3flux(IM3,ke+1,j,i) = x3flux(IM3,ks,j,i);
        if (NON_BAROTROPIC_EOS) {
          x3flux(IEN,ks,j,i) += prim(IM1,ks,j,i)*flx1 +
                                prim(IM2,ks,j,i)*flx2 +
                                prim(IM3,ks,j,i)*flx3;
          x3flux(IEN,ke+1,j,i) = x3flux(IEN,ks,j,i);
        }
      }
    }
  }
  return;
}


//-------------------------------------------------------------------------------------
// Calculate divergence of momenta

void HydroDiffusion::Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv) {
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  int il = is-1; int iu = ie+1;
  int jl, ju, kl, ku;
  Real area_p1, area;
  Real vel_p1, vel;

  if (pmb_->block_size.nx2 == 1) // 1D
    jl=js, ju=je, kl=ks, ku=ke;
  else if (pmb_->block_size.nx3 == 1) // 2D
    jl=js-1, ju=je+1, kl=ks, ku=ke;
  else // 3D
    jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // calculate x1-flux divergence
      pmb_->pcoord->Face1Area(k, j, il, iu+1, x1area_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(prim(IM1,k,j,i+1) + prim(IM1,k,j,i  ));
        vel     = 0.5*(prim(IM1,k,j,i  ) + prim(IM1,k,j,i-1));
        divv(k,j,i) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      if (pmb_->block_size.nx2 > 1) {
        pmb_->pcoord->Face2Area(k, j  , il, iu, x2area_);
        pmb_->pcoord->Face2Area(k, j+1, il, iu, x2area_p1_);
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          area_p1 = x2area_p1_(i);
          area    = x2area_(i);
          vel_p1  = 0.5*(prim(IM2,k,j+1,i) + prim(IM2,k,j  ,i));
          vel     = 0.5*(prim(IM2,k,j  ,i) + prim(IM2,k,j-1,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      if (pmb_->block_size.nx3 > 1) {
        pmb_->pcoord->Face3Area(k  , j, il, iu, x3area_);
        pmb_->pcoord->Face3Area(k+1, j, il, iu, x3area_p1_);
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          area_p1 = x3area_p1_(i);
          area    = x3area_(i);
          vel_p1  = 0.5*(prim(IM3,k+1,j,i) + prim(IM3, k  ,j,i));
          vel     = 0.5*(prim(IM3,k  ,j,i) + prim(IM3, k-1,j,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      pmb_->pcoord->CellVolume(k,j,il,iu,vol_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        divv(k,j,i) = divv(k,j,i)/vol_(i);
      }
    }
  }

  return;
}

// v_{x1;x1}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdx(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    len(i) = 2.0*(prim(IM1,k,j,i) - prim(IM1,k,j,i-1)) / pco_->dx1v(i-1);
  }
  return;
}

// v_{x2;x1}+v_{x1;x2}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdy(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx2 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h2f(i)
          * (prim(IM2,k,j,i)/pco_->h2v(i) - prim(IM2,k,j,i-1)/pco_->h2v(i-1))
          / pco_->dx1v(i-1)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          + 0.5*(   (prim(IM1,k,j+1,i) + prim(IM1,k,j+1,i-1))
                  - (prim(IM1,k,j-1,i) + prim(IM1,k,j-1,i-1)) )
          / pco_->h2f(i)
          / (pco_->dx2v(j-1) + pco_->dx2v(j));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h2f(i)
          * ( prim(IM2,k,j,i)/pco_->h2v(i) - prim(IM2,k,j,i-1)/pco_->h2v(i-1) )
          / pco_->dx1v(i-1);
    }
  }
  return;
}

// v_{x3;x1}+v_{x1;x3}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdz(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx3 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h31f(i)
          * (prim(IM3,k,j,i)/pco_->h31v(i) - prim(IM3,k,j,i-1)/pco_->h31v(i-1))
          / pco_->dx1v(i-1)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          + 0.5*(   (prim(IM1,k+1,j,i) + prim(IM1,k+1,j,i-1))
                  - (prim(IM1,k-1,j,i) + prim(IM1,k-1,j,i-1)) )
          / pco_->h31f(i)/pco_->h32v(j) // note, more terms than FaceXdy() line
          / (pco_->dx3v(k-1) + pco_->dx3v(k));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = pco_->h31f(i)
          * ( prim(IM3,k,j,i)/pco_->h31v(i) - prim(IM3,k,j,i-1)/pco_->h31v(i-1) )
          / pco_->dx1v(i-1);
  }
  return;
}

// v_{x1;x2}+v_{x2;x1}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdx(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx2 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = (prim(IM1,k,j,i) - prim(IM1,k,j-1,i)) / pco_->h2v(i) / pco_->dx2v(j-1)
          + pco_->h2v(i)*0.5*(  (prim(IM2,k,j,i+1) + prim(IM2,k,j-1,i+1)) /pco_->h2v(i+1)
                              - (prim(IM2,k,j,i-1) + prim(IM2,k,j-1,i-1)) /pco_->h2v(i-1)
                                ) / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h2v(i)
          * ( prim(IM2,k,j,i+1)/pco_->h2v(i+1) - prim(IM2,k,j,i-1)/pco_->h2v(i-1) )
          / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  }
  return;
}

// v_{x2;x2}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx2 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = 2.0*(prim(IM2,k,j,i) - prim(IM2,k,j-1,i)) / pco_->h2v(i) / pco_->dx2v(j-1)
          + (prim(IM1,k,j,i) + prim(IM1,k,j-1,i)) / pco_->h2v(i) * pco_->dh2vd1(i);
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = 2.0*prim(IM1,k,j,i) / pco_->h2v(i) * pco_->dh2vd1(i);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx3 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h32f(j)
          * ( prim(IM3,k,j,i)/pco_->h32v(j) - prim(IM3,k,j-1,i)/pco_->h32v(j-1) )
          / pco_->h2v(i) / pco_->dx2v(j-1)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          + 0.5*(    (prim(IM2,k+1,j,i) + prim(IM2,k+1,j-1,i))
                   - (prim(IM2,k-1,j,i) + prim(IM2,k-1,j-1,i)) )
          / pco_->h31v(i)
          / pco_->h32f(j) / (pco_->dx3v(k-1) + pco_->dx3v(k));
    }
  } else if (pmb_->block_size.nx2 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h32f(j)
          * ( prim(IM3,k,j,i)/pco_->h32v(j) - prim(IM3,k,j-1,i)/pco_->h32v(j-1) )
          / pco_->h2v(i) / pco_->dx2v(j-1);
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = 0.0;
  }
  return;
}

// v_{x1;x3}+v_{x3;x1}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdx(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx3 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = (prim(IM1,k,j,i) - prim(IM1,k-1,j,i))/pco_->dx3v(k-1)
          + 0.5*pco_->h31v(i)*( (prim(IM3,k,j,i+1) + prim(IM3,k-1,j,i+1))/pco_->h31v(i+1)
                               -(prim(IM3,k,j,i-1) + prim(IM3,k-1,j,i-1))/pco_->h31v(i-1)
                                ) / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h31v(i)
          * ( prim(IM3,k,j,i+1)/pco_->h31v(i+1) - prim(IM3,k,j,i-1)/pco_->h31v(i-1) )
          / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  }
  return;
}

// v_{x2;x3}+v_{x3;x2}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdy(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx3 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = (prim(IM2,k,j,i) - prim(IM2,k-1,j,i))
          / pco_->h31v(i) / pco_->h32v(j) / pco_->dx3v(k-1)
          + 0.5*pco_->h32v(j)
          * ( (prim(IM3,k,j+1,i) + prim(IM3,k-1,j+1,i))/pco_->h32v(j+1)
             -(prim(IM3,k,j-1,i) + prim(IM3,k-1,j-1,i))/pco_->h32v(j-1) )
          / pco_->h2v(i) / (pco_->dx2v(j-1) + pco_->dx2v(j));
    }
  } else if (pmb_->block_size.nx2 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h32v(j)
          * ( prim(IM3,k,j+1,i)/pco_->h32v(j+1) - prim(IM3,k,j-1,i)/pco_->h32v(j-1) )
          / pco_->h2v(i) / (pco_->dx2v(j-1) + pco_->dx2v(j));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = 0.0;
  }
  return;
}

// v_{x3;x3}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdz(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->block_size.nx3 > 1) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = 2.0*(prim(IM3,k,j,i) - prim(IM3,k-1,j,i))
          / pco_->dx3v(k-1) / pco_->h31v(i) / pco_->h32v(j)
          + ((prim(IM1,k,j,i) + prim(IM1,k-1,j,i))
             * pco_->dh31vd1(i)/pco_->h31v(i))
          + ((prim(IM2,k,j,i) + prim(IM2,k-1,j,i))
             * pco_->dh32vd2(j)/pco_->h32v(j)/pco_->h2v(i));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = 2.0*prim(IM1,k,j,i)*pco_->dh31vd1(i)/pco_->h31v(i)
          +    2.0*prim(IM2,k,j,i)*pco_->dh32vd2(j)/pco_->h32v(j)/pco_->h2v(i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ViscousFlux_aniso
//  \brief Calculate anisotropic viscous stress as fluxes
//
//   Note: difference with athena, delta_p is nu_aniso*BBdV, rather than 3*nu_aniso*BBdV
//   WARNING!!!! Does not include coordinate scale factors for now
//   WARNING!!!! 2D not implemented yet
void HydroDiffusion::ViscousFlux_aniso(const AthenaArray<Real> &prim,
                                       const AthenaArray<Real> &cons,
                                       AthenaArray<Real> *visflx,
                                       const FaceField &b,
                                       const AthenaArray<Real> &bcc) {
  
  AthenaArray<Real> &x1flux=visflx[X1DIR];
  AthenaArray<Real> &x2flux=visflx[X2DIR];
  AthenaArray<Real> &x3flux=visflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  // Add these below for vectorization
//  Real nu1, denf;
//  Real Bx, By, Bz, B02, BBdV, delta_p;
//  Real dx, dy, dz;
//  Real dVxdx, dVydx, dVzdx, dVxdy, dVydy, dVzdy, dVxdz, dVydz, dVzdz;
  
//  Divv(prim, divv_); // Should compare this to using calculated dVdx etc.
  
  // Calculate the flux across each face.
  // i-direction
//  if (pmb_->block_size.nx3 == 1) // 2D MHD limits
//    jl=js-1, ju=je+1, kl=ks, ku=ke;
//  else // 3D MHD limits
  
  jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // One possibility #pragma omp simd private(Bx, By, Bz, B02, BBdV, delta_p, dVxdx, dVydx, dVzdx, dVxdy, dVydy, dVzdy, dVxdz, dVydz, dVzdz,dx, dy, dz, nu1, denf)
      // Better is probably to initialize them in the loop, e.g., Real dVxdy
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        Real Bx = b.x1f(k,j,i);
        Real By = 0.5*(bcc(IB2,k,j,i) + bcc(IB2,k,j,i-1));
        Real Bz = 0.5*(bcc(IB3,k,j,i) + bcc(IB3,k,j,i-1));
        Real B02 = Bx*Bx + By*By + Bz*Bz;
        B02 = MAX(B02,TINY_NUMBER);  // In case B02=0
        
        Real dx = pco_->dx1v(i-1);
        Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
        Real dz = 0.5*(pco_->dx3v(k-1) + pco_->dx3v(k));
        
        // Gradient dVx/dx at i-1/2, j, k
        Real dVxdx = (prim(IM1,k,j,i) - prim(IM1,k,j,i-1)) / dx;
        // Gradient dVy/dx at i-1/2, j, k
        Real dVydx = (prim(IM2,k,j,i) - prim(IM2,k,j,i-1)) / dx;
        // Gradient dVz/dx at i-1/2, j, k
        Real dVzdx = (prim(IM3,k,j,i) - prim(IM3,k,j,i-1)) / dx;
        
        // Gradient dVx/dy at i-1/2, j, k
        Real dVxdy = limiter4(prim(IM1,k,j+1,i  ) - prim(IM1,k,j  ,i  ),
                         prim(IM1,k,j  ,i  ) - prim(IM1,k,j-1,i  ),
                         prim(IM1,k,j+1,i-1) - prim(IM1,k,j  ,i-1),
                         prim(IM1,k,j  ,i-1) - prim(IM1,k,j-1,i-1));
        dVxdy /= dy;
        // Gradient dVy/dy at i-1/2, j, k
        Real dVydy = limiter4(prim(IM2,k,j+1,i  ) - prim(IM2,k,j  ,i  ),
                         prim(IM2,k,j  ,i  ) - prim(IM2,k,j-1,i  ),
                         prim(IM2,k,j+1,i-1) - prim(IM2,k,j  ,i-1),
                         prim(IM2,k,j  ,i-1) - prim(IM2,k,j-1,i-1));
        dVydy /= dy;
        // Gradient dVz/dy at i-1/2, j, k
        Real dVzdy = limiter4(prim(IM3,k,j+1,i  ) - prim(IM3,k,j  ,i  ),
                         prim(IM3,k,j  ,i  ) - prim(IM3,k,j-1,i  ),
                         prim(IM3,k,j+1,i-1) - prim(IM3,k,j  ,i-1),
                         prim(IM3,k,j  ,i-1) - prim(IM3,k,j-1,i-1));
        dVzdy /= dy;
        
        // Gradient dVx/dz at i-1/2, j, k
        Real dVxdz = limiter4(prim(IM1,k+1,j,i  ) - prim(IM1,k  ,j,i  ),
                         prim(IM1,k  ,j,i  ) - prim(IM1,k-1,j,i  ),
                         prim(IM1,k+1,j,i-1) - prim(IM1,k  ,j,i-1),
                         prim(IM1,k  ,j,i-1) - prim(IM1,k-1,j,i-1));
        dVxdz /=dz;
        // Gradient dVy/dz at i-1/2, j, k
        Real dVydz = limiter4(prim(IM2,k+1,j,i  ) - prim(IM2,k  ,j,i  ),
                         prim(IM2,k  ,j,i  ) - prim(IM2,k-1,j,i  ),
                         prim(IM2,k+1,j,i-1) - prim(IM2,k  ,j,i-1),
                         prim(IM2,k  ,j,i-1) - prim(IM2,k-1,j,i-1));
        dVydz /=dz;
        // Gradient dVx/dz at i-1/2, j, k
        Real dVzdz = limiter4(prim(IM3,k+1,j,i  ) - prim(IM3,k  ,j,i  ),
                         prim(IM3,k  ,j,i  ) - prim(IM3,k-1,j,i  ),
                         prim(IM3,k+1,j,i-1) - prim(IM3,k  ,j,i-1),
                         prim(IM3,k  ,j,i-1) - prim(IM3,k-1,j,i-1));
        dVzdz /=dz;
        
        
        // Compute BB:GV
        Real BBdV = Bx*(Bx*dVxdx + By*dVxdy + Bz*dVxdz) +
              By*(Bx*dVydx + By*dVydy + Bz*dVydz) +
              Bz*(Bx*dVzdx + By*dVzdy + Bz*dVzdz);
        BBdV /= B02;
        
        // NB: nu_aniso is a kinematic viscosity, since we multiply by density
        Real nu1  = 0.5*(nu(ANI,k,j,i)   + nu(ANI,k,j,i-1));
        Real denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
        
        Real delta_p = nu1*denf*(BBdV - ONE_3RD*(dVxdx + dVydy + dVzdz));
        
        // Apply mirror and/or firehose limiters
        if (mirror_limit && delta_p>0.5*B02)
          delta_p = 0.5*B02;
        if (firehose_limit && delta_p<-1.0*B02)
          delta_p = -1.0*B02;
        
        // NB: Flux sign seems to be opposite to old athena
        x1flux(IM1,k,j,i) += -delta_p*(Bx*Bx/B02 - ONE_3RD);
        x1flux(IM2,k,j,i) += -delta_p*(Bx*By/B02);
        x1flux(IM3,k,j,i) += -delta_p*(Bx*Bz/B02);
        
        
        if (NON_BAROTROPIC_EOS)
          x1flux(IEN,k,j,i) +=
              0.5*((prim(IM1,k,j,i) + prim(IM1,k,j,i-1))*x1flux(IM1,k,j,i) +
                (prim(IM2,k,j,i) + prim(IM2,k,j,i-1))*x1flux(IM2,k,j,i) +
                (prim(IM3,k,j,i) + prim(IM3,k,j,i-1))*x1flux(IM3,k,j,i));
        
      }
    }
  }
  
//  // j-direction
//  if (pmb_->block_size.nx3 == 1) // 2D MHD limits
//    il=is-1, iu=ie+1, kl=ks, ku=ke;
//  else // 3D MHD limits
    il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real Bx = 0.5*(bcc(IB1,k,j,i) + bcc(IB1,k,j-1,i));
        Real By = b.x2f(k,j,i);
        Real Bz = 0.5*(bcc(IB3,k,j,i) + bcc(IB3,k,j-1,i));
        Real B02 = Bx*Bx + By*By + Bz*Bz;
        B02 = MAX(B02,TINY_NUMBER);  // In case B02=0
        
        Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
        Real dy = pco_->dx2v(j-1);
        Real dz = 0.5*(pco_->dx3v(k-1) + pco_->dx3v(k));
        
        // Gradient dVx/dx at i, j-1/2, k
        Real dVxdx = limiter4(prim(IM1,k,j  ,i+1) - prim(IM1,k,j  ,i  ),
                         prim(IM1,k,j  ,i  ) - prim(IM1,k,j  ,i-1),
                         prim(IM1,k,j-1,i+1) - prim(IM1,k,j-1,i  ),
                         prim(IM1,k,j-1,i  ) - prim(IM1,k,j-1,i-1));
        dVxdx /= dx;
        // Gradient dVx/dx at i, j-1/2, k
        Real dVydx = limiter4(prim(IM2,k,j  ,i+1) - prim(IM2,k,j  ,i  ),
                         prim(IM2,k,j  ,i  ) - prim(IM2,k,j  ,i-1),
                         prim(IM2,k,j-1,i+1) - prim(IM2,k,j-1,i  ),
                         prim(IM2,k,j-1,i  ) - prim(IM2,k,j-1,i-1));
        dVydx /= dx;
        // Gradient dVx/dx at i, j-1/2, k
        Real dVzdx = limiter4(prim(IM3,k,j  ,i+1) - prim(IM3,k,j  ,i  ),
                         prim(IM3,k,j  ,i  ) - prim(IM3,k,j  ,i-1),
                         prim(IM3,k,j-1,i+1) - prim(IM3,k,j-1,i  ),
                         prim(IM3,k,j-1,i  ) - prim(IM3,k,j-1,i-1));
        dVzdx /= dx;
        
        // Gradient dVx/dy at i, j-1/2, k
        Real dVxdy = (prim(IM1,k,j,i) - prim(IM1,k,j-1,i)) / dy;
        // Gradient dVx/dy at i, j-1/2, k
        Real dVydy = (prim(IM2,k,j,i) - prim(IM2,k,j-1,i)) / dy;
        // Gradient dVx/dy at i, j-1/2, k
        Real dVzdy = (prim(IM3,k,j,i) - prim(IM3,k,j-1,i)) / dy;
        
        // Gradient dVx/dz at i, j-1/2, k
        Real dVxdz = limiter4(prim(IM1,k+1,j  ,i) - prim(IM1,k  ,j  ,i),
                         prim(IM1,k  ,j  ,i) - prim(IM1,k-1,j  ,i),
                         prim(IM1,k+1,j-1,i) - prim(IM1,k  ,j-1,i),
                         prim(IM1,k  ,j-1,i) - prim(IM1,k-1,j-1,i));
        dVxdz /= dz;
        // Gradient dVx/dz at i, j-1/2, k
        Real dVydz = limiter4(prim(IM2,k+1,j  ,i) - prim(IM2,k  ,j  ,i),
                         prim(IM2,k  ,j  ,i) - prim(IM2,k-1,j  ,i),
                         prim(IM2,k+1,j-1,i) - prim(IM2,k  ,j-1,i),
                         prim(IM2,k  ,j-1,i) - prim(IM2,k-1,j-1,i));
        dVydz /= dz;
        // Gradient dVx/dz at i, j-1/2, k
        Real dVzdz = limiter4(prim(IM3,k+1,j  ,i) - prim(IM3,k  ,j  ,i),
                         prim(IM3,k  ,j  ,i) - prim(IM3,k-1,j  ,i),
                         prim(IM3,k+1,j-1,i) - prim(IM3,k  ,j-1,i),
                         prim(IM3,k  ,j-1,i) - prim(IM3,k-1,j-1,i));
        dVzdz /= dz;
        
        // Compute BB:GV
        Real BBdV = Bx*(Bx*dVxdx + By*dVxdy + Bz*dVxdz) +
              By*(Bx*dVydx + By*dVydy + Bz*dVydz) +
              Bz*(Bx*dVzdx + By*dVzdy + Bz*dVzdz);
        BBdV /= B02;
        
        Real nu1  = 0.5*(nu(ANI,k,j,i)   + nu(ANI,k,j-1,i));
        Real denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j-1,i));
        
        Real delta_p = nu1*denf*(BBdV - ONE_3RD*(dVxdx + dVydy + dVzdz));
        
        // Apply mirror and/or firehose limiters
        if (mirror_limit && delta_p>0.5*B02)
          delta_p = 0.5*B02;
        if (firehose_limit && delta_p<-1.0*B02)
          delta_p = -1.0*B02;
                
        x2flux(IM1,k,j,i) += -delta_p*(By*Bx/B02);
        x2flux(IM2,k,j,i) += -delta_p*(By*By/B02 - ONE_3RD);
        x2flux(IM3,k,j,i) += -delta_p*(By*Bz/B02);
        
        if (NON_BAROTROPIC_EOS)
          x2flux(IEN,k,j,i) +=
          0.5*((prim(IM1,k,j,i) + prim(IM1,k,j-1,i))*x2flux(IM1,k,j,i) +
               (prim(IM2,k,j,i) + prim(IM2,k,j-1,i))*x2flux(IM2,k,j,i) +
               (prim(IM3,k,j,i) + prim(IM3,k,j-1,i))*x2flux(IM3,k,j,i));
        
      }
    }
  }
  
  // k-direction
  il=is-1, iu=ie+1, jl=js-1, ju=je+1;
  
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real Bx = 0.5*(bcc(IB1,k,j,i) + bcc(IB1,k-1,j,i));
        Real By = 0.5*(bcc(IB2,k,j,i) + bcc(IB2,k-1,j,i));
        Real Bz = b.x3f(k,j,i);
        Real B02 = Bx*Bx + By*By + Bz*Bz;
        B02 = MAX(B02,TINY_NUMBER);  // In case B02=0
        
        Real dx = 0.5*(pco_->dx1v(i-1) + pco_->dx1v(i));
        Real dy = 0.5*(pco_->dx2v(j-1) + pco_->dx2v(j));
        Real dz = pco_->dx3v(k-1);
        
        // Gradient dVx/dx at i, j, k-1/2
        Real dVxdx = limiter4(prim(IM1,k  ,j,i+1) - prim(IM1,k  ,j,i  ),
                         prim(IM1,k  ,j,i  ) - prim(IM1,k  ,j,i-1),
                         prim(IM1,k-1,j,i+1) - prim(IM1,k-1,j,i  ),
                         prim(IM1,k-1,j,i  ) - prim(IM1,k-1,j,i-1));
        dVxdx /= dx;
        /// Gradient dVz/dx at i, j, k-1/2
        Real dVydx = limiter4(prim(IM2,k  ,j,i+1) - prim(IM2,k  ,j,i  ),
                         prim(IM2,k  ,j,i  ) - prim(IM2,k  ,j,i-1),
                         prim(IM2,k-1,j,i+1) - prim(IM2,k-1,j,i  ),
                         prim(IM2,k-1,j,i  ) - prim(IM2,k-1,j,i-1));
        dVydx /= dx;
        // Gradient dVz/dx at i, j, k-1/2
        Real dVzdx = limiter4(prim(IM3,k  ,j,i+1) - prim(IM3,k  ,j,i  ),
                         prim(IM3,k  ,j,i  ) - prim(IM3,k  ,j,i-1),
                         prim(IM3,k-1,j,i+1) - prim(IM3,k-1,j,i  ),
                         prim(IM3,k-1,j,i  ) - prim(IM3,k-1,j,i-1));
        dVzdx /= dx;
        
        // Gradient dVx/dy at i, j, k-1/2
        Real dVxdy = limiter4(prim(IM1,k  ,j+1,i) - prim(IM1,k  ,j  ,i),
                         prim(IM1,k  ,j  ,i) - prim(IM1,k  ,j-1,i),
                         prim(IM1,k-1,j+1,i) - prim(IM1,k-1,j  ,i),
                         prim(IM1,k-1,j  ,i) - prim(IM1,k-1,j-1,i));
        dVxdy /= dy;
        // Gradient dVx/dy at i, j, k-1/2
        Real dVydy = limiter4(prim(IM2,k  ,j+1,i) - prim(IM2,k  ,j  ,i),
                         prim(IM2,k  ,j  ,i) - prim(IM2,k  ,j-1,i),
                         prim(IM2,k-1,j+1,i) - prim(IM2,k-1,j  ,i),
                         prim(IM2,k-1,j  ,i) - prim(IM2,k-1,j-1,i));
        dVydy /= dy;
        // Gradient dVx/dy at i, j, k-1/2
        Real dVzdy = limiter4(prim(IM3,k  ,j+1,i) - prim(IM3,k  ,j  ,i),
                         prim(IM3,k  ,j  ,i) - prim(IM3,k  ,j-1,i),
                         prim(IM3,k-1,j+1,i) - prim(IM3,k-1,j  ,i),
                         prim(IM3,k-1,j  ,i) - prim(IM3,k-1,j-1,i));
        dVzdy /= dy;
        
        // Gradient dVx/dz at i, j, k-1/2
        Real dVxdz = (prim(IM1,k,j,i) - prim(IM1,k-1,j,i)) / dz;
        // Gradient dVy/dz at i, j, k-1/2
        Real dVydz = (prim(IM2,k,j,i) - prim(IM2,k-1,j,i)) / dz;
        // Gradient dVz/dz at i, j, k-1/2
        Real dVzdz = (prim(IM3,k,j,i) - prim(IM3,k-1,j,i)) / dz;
        
        // Compute BB:GV
        Real BBdV = Bx*(Bx*dVxdx + By*dVxdy + Bz*dVxdz) +
              By*(Bx*dVydx + By*dVydy + Bz*dVydz) +
              Bz*(Bx*dVzdx + By*dVzdy + Bz*dVzdz);
        BBdV /= B02;
        
        Real nu1  = 0.5*(nu(ANI,k,j,i)   + nu(ANI,k-1,j,i));
        Real denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k-1,j,i));
        
        Real delta_p = nu1*denf*(BBdV - ONE_3RD*(dVxdx + dVydy + dVzdz));
        
        // Apply mirror and/or firehose limiters
        if (mirror_limit && delta_p>0.5*B02)
          delta_p = 0.5*B02;
        if (firehose_limit && delta_p<-1.0*B02)
          delta_p = -1.0*B02;
        
        x3flux(IM1,k,j,i) += -delta_p*(Bz*Bx/B02);
        x3flux(IM2,k,j,i) += -delta_p*(Bz*By/B02);
        x3flux(IM3,k,j,i) += -delta_p*(Bz*Bz/B02 - ONE_3RD);
        
        if (NON_BAROTROPIC_EOS)
          x3flux(IEN,k,j,i) +=
          0.5*((prim(IM1,k,j,i) + prim(IM1,k-1,j,i))*x3flux(IM1,k,j,i) +
               (prim(IM2,k,j,i) + prim(IM2,k-1,j,i))*x3flux(IM2,k,j,i) +
               (prim(IM3,k,j,i) + prim(IM3,k-1,j,i))*x3flux(IM3,k,j,i));
        
      }
    }
  }
  
  
  return;
}

//----------------------------------------------------------------------------------------
// constant viscosity

void ConstViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &bcc, int is, int ie, int js, int je,
                    int ks, int ke) {
  if (phdif->nu_iso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->nu(ISO,k,j,i) = phdif->nu_iso;
      }
    }
  }
  if (phdif->nu_aniso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->nu(ANI,k,j,i) = phdif->nu_aniso;
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

static Real minmod(const Real A, const Real B)
{
  if (A*B > 0) {
    if (A > 0) {
      return MIN(A,B);
    } else {
      return MAX(A,B);
    }
  } else {
    return 0.0;
  }
}




