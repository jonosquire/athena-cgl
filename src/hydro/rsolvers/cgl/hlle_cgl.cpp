//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlle_cgl.cpp
//  \brief HLLE Riemann solver for CGL equation of state (with MHD).  See the hydro version for details.
//  Added J Squire 2018

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn

void Hydro::RiemannSolver(const int kl, const int ku, const int jl, const int ju,
                          const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NWAVE)],wri[(NWAVE)],wroe[(NWAVE)],fl[(NWAVE)],fr[(NWAVE)],flxi[(NWAVE)];
  Real gm1 = pmy_block->peos->GetGamma() - 1.0;
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd private(wli,wri,wroe,fl,fr,flxi)
      for (int i=il; i<=iu; ++i) {
        
        //--- Step 1.  Load L/R states into local variables
        
        wli[IDN]=wl(IDN,k,j,i);
        wli[IVX]=wl(ivx,k,j,i);
        wli[IVY]=wl(ivy,k,j,i);
        wli[IVZ]=wl(ivz,k,j,i);
        wli[IPR]=wl(IPR,k,j,i);
        wli[IPP]=wl(IPP,k,j,i);
        wli[IBY]=wl(IBY,k,j,i);
        wli[IBZ]=wl(IBZ,k,j,i);
        
        wri[IDN]=wr(IDN,k,j,i);
        wri[IVX]=wr(ivx,k,j,i);
        wri[IVY]=wr(ivy,k,j,i);
        wri[IVZ]=wr(ivz,k,j,i);
        wri[IPR]=wr(IPR,k,j,i);
        wri[IPP]=wr(IPP,k,j,i);
        wri[IBY]=wr(IBY,k,j,i);
        wri[IBZ]=wr(IBZ,k,j,i);
        
        Real bxi = bx(k,j,i);
        
        //--- For now, we will just use the max/min to estimate wave speeds 10.48 of Toro
        
        Real pbl = 0.5*(bxi*bxi + SQR(wli[IBY]) + SQR(wli[IBZ]));
        Real pbr = 0.5*(bxi*bxi + SQR(wri[IBY]) + SQR(wri[IBZ]));
        
        Real dpobsql = 1. + (wli[IPP] - wli[IPR])/(2.*pbl);
        Real dpobsqr = 1. + (wri[IPP] - wri[IPR])/(2.*pbr);
        
        Real el = wli[IPR]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX])+SQR(wli[IVY])+SQR(wli[IVZ])) + pbl;
        Real er = wri[IPR]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX])+SQR(wri[IVY])+SQR(wri[IVZ])) + pbr;
        
        Real mul = wli[IPP]/std::sqrt(2.*pbl);
        Real mur = wri[IPP]/std::sqrt(2.*pbr);
        
        //--- Step 3.  Compute fast magnetosonic speed in L,R, and Roe-averaged states
        
        Real cl = pmy_block->peos->FastMagnetosonicSpeed(wli,bxi);
        Real cr = pmy_block->peos->FastMagnetosonicSpeed(wri,bxi);
        
        //--- Step 4.  Compute the max/min wave speeds based on L/R and Roe-averaged values
        
        Real al = std::min((wri[IVX] - cr),(wli[IVX] - cl));
        Real ar = std::max((wri[IVX] + cr),(wli[IVX] + cl));
        
        Real bp = ar > 0.0 ? ar : 0.0;
        Real bm = al < 0.0 ? al : 0.0;
        
        //--- Step 5.  Compute L/R fluxes along the lines bm/bp: F_L - (S_L)U_L; F_R - (S_R)U_R
        
        Real vxl = wli[IVX] - bm;
        Real vxr = wri[IVX] - bp;
        
        fl[IDN] = wli[IDN]*vxl;
        fr[IDN] = wri[IDN]*vxr;
        
        fl[IVX] = wli[IDN]*wli[IVX]*vxl + pbl + wli[IPR] - SQR(bxi)*dpobsql;
        fr[IVX] = wri[IDN]*wri[IVX]*vxr + pbr + wri[IPR] - SQR(bxi)*dpobsqr;
        
        fl[IVY] = wli[IDN]*wli[IVY]*vxl - bxi*wli[IBY]*dpobsql;
        fr[IVY] = wri[IDN]*wri[IVY]*vxr - bxi*wri[IBY]*dpobsqr;
        
        fl[IVZ] = wli[IDN]*wli[IVZ]*vxl - bxi*wli[IBZ]*dpobsql;
        fr[IVZ] = wri[IDN]*wri[IVZ]*vxr - bxi*wri[IBZ]*dpobsqr;
        
        fl[IEN] = el*vxl + wli[IVX]*(wli[IPR] + pbl);
        fr[IEN] = er*vxr + wri[IVX]*(wri[IPR] + pbr);
        fl[IEN] -= bxi*(bxi*wli[IVX] + wli[IBY]*wli[IVY] + wli[IBZ]*wli[IVZ])*dpobsql;
        fr[IEN] -= bxi*(bxi*wri[IVX] + wri[IBY]*wri[IVY] + wri[IBZ]*wri[IVZ])*dpobsqr;
        
        fl[IMU] = mul*vxl;
        fr[IMU] = mur*vxr;
        
        fl[IBY] = wli[IBY]*vxl - bxi*wli[IVY];
        fr[IBY] = wri[IBY]*vxr - bxi*wri[IVY];
        
        fl[IBZ] = wli[IBZ]*vxl - bxi*wli[IVZ];
        fr[IBZ] = wri[IBZ]*vxr - bxi*wri[IVZ];
        
        //--- Step 6.  Compute the HLLE flux at interface.
        
        Real tmp=0.0;
        if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);
        
        flxi[IDN] = 0.5*(fl[IDN]+fr[IDN]) + (fl[IDN]-fr[IDN])*tmp;
        flxi[IVX] = 0.5*(fl[IVX]+fr[IVX]) + (fl[IVX]-fr[IVX])*tmp;
        flxi[IVY] = 0.5*(fl[IVY]+fr[IVY]) + (fl[IVY]-fr[IVY])*tmp;
        flxi[IVZ] = 0.5*(fl[IVZ]+fr[IVZ]) + (fl[IVZ]-fr[IVZ])*tmp;
        flxi[IEN] = 0.5*(fl[IEN]+fr[IEN]) + (fl[IEN]-fr[IEN])*tmp;
        flxi[IMU] = 0.5*(fl[IMU]+fr[IMU]) + (fl[IMU]-fr[IMU])*tmp;
        flxi[IBY] = 0.5*(fl[IBY]+fr[IBY]) + (fl[IBY]-fr[IBY])*tmp;
        flxi[IBZ] = 0.5*(fl[IBZ]+fr[IBZ]) + (fl[IBZ]-fr[IBZ])*tmp;
        
        flx(IDN,k,j,i) = flxi[IDN];
        flx(ivx,k,j,i) = flxi[IVX];
        flx(ivy,k,j,i) = flxi[IVY];
        flx(ivz,k,j,i) = flxi[IVZ];
        flx(IEN,k,j,i) = flxi[IEN];
        flx(IMU,k,j,i) = flxi[IMU];
        ey(k,j,i) = -flxi[IBY];
        ez(k,j,i) =  flxi[IBZ];
      }
    }}
  return;
}
