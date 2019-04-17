//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlld_cgl.cpp
//  \brief HLLD/E Riemann solver for CGL.
//
// REFERENCES:
// To come... similar to
// - A. Mignone, "A simple and accurate Riemann solver for isothermal MHD", JPC, 225,
//   1427 (2007)

// C/C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"

// container to store (density, momentum, tranverse magnetic field)
// minimizes changes required to adopt hlld_iso solver
typedef struct Cons1D {
  Real d,mx,my,mz,e,mu,by,bz;
} Cons1D;

#define SMALL_NUMBER 1.0e-8

//----------------------------------------------------------------------------------------
//! \fn

void Hydro::RiemannSolver(const int kl, const int ku, const int jl, const int ju,
                          const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real flxi[(NWAVE)];             // temporary variable to store flux
  Real wli[(NWAVE)],wri[(NWAVE)]; // L/R states, primitive variables (input)
  Real spd[5];                    // signal speeds, left to right
  
  Real dfloor = pmy_block->peos->GetDensityFloor();
  Real pfloor = pmy_block->peos->GetPressureFloor();
  Real bmag_floor = pmy_block->peos->GetBFieldFloor();
  Real bsq_floor = SQR(bmag_floor);
  // Value of 1+∆p/B^2 below which we interpolate to HLLE solver
  Real fh_floor = pmy_block->peos->GetFHFloor();
  //
  Real anti_diff_alpha = pmy_block->peos->GetAntiDiffAlpha();
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd private(flxi,wli,wri,spd)
      for (int i=il; i<=iu; ++i) {
        Cons1D ul,ur;                   // L/R states, conserved variables (computed)
        Cons1D hll;            // hll variables
        Cons1D  ulst, urst;           // L and R starred variables
        Cons1D fl,fr;                   // Fluxes for left & right states
        Cons1D fhll;                   // HLL flux
        
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
        
        // Compute L/R states for conserved variables, and other useful auxilliaries
        Real bxsq = bxi*bxi;
        Real bsql = bxsq + SQR(wli[IBY]) + SQR(wli[IBZ]);
        Real bsqr = bxsq + SQR(wri[IBY]) + SQR(wri[IBZ]);
        bsql = bsql > bsq_floor ? bsql : bsq_floor;
        bsqr = bsqr > bsq_floor ? bsqr : bsq_floor;
        Real bmagl = std::sqrt(bsql);
        Real bmagr = std::sqrt(bsqr);
        
        Real kel = 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
        Real ker = 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
        
        Real fhl = 1. + (wli[IPP] - wli[IPR])/bsql;
        Real fhr = 1. + (wri[IPP] - wri[IPR])/bsqr;
        
        ul.d  = wli[IDN];
        ul.mx = wli[IVX]*ul.d;
        ul.my = wli[IVY]*ul.d;
        ul.mz = wli[IVZ]*ul.d;
        ul.e  = wli[IPP] + 0.5*wli[IPR] + 0.5*bsql + kel;
        ul.mu = wli[IPP]/bmagl;
        ul.by = wli[IBY];
        ul.bz = wli[IBZ];
        
        ur.d  = wri[IDN];
        ur.mx = wri[IVX]*ur.d;
        ur.my = wri[IVY]*ur.d;
        ur.mz = wri[IVZ]*ur.d;
        ur.e  = wri[IPP] + 0.5*wri[IPR] + 0.5*bsqr + ker;
        ur.mu = wri[IPP]/bmagr;
        ur.by = wri[IBY];
        ur.bz = wri[IBZ];
        
        //--- Step 2.  Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)
        
        Real cfl = pmy_block->peos->FastMagnetosonicSpeed(wli,bxi);
        Real cfr = pmy_block->peos->FastMagnetosonicSpeed(wri,bxi);
        
        spd[0] = std::min( wli[IVX]-cfl, wri[IVX]-cfr );
        spd[4] = std::max( wli[IVX]+cfl, wri[IVX]+cfr );
        
        //--- Step 3.  Compute L/R fluxes
        
        Real bxdpl = bxi*fhl;
        Real bxdpr = bxi*fhr;
        
        fl.d  = ul.mx;
        fl.mx = ul.mx*wli[IVX] + wli[IPP] + 0.5*bsql - bxdpl*bxi;
        fl.my = ul.my*wli[IVX] - ul.by*bxdpl;
        fl.mz = ul.mz*wli[IVX] - ul.bz*bxdpl;
        fl.e  = (ul.e + wli[IPP] + 0.5*bsql)*wli[IVX] -
                  (bxi*wli[IVX] + ul.by*wli[IVY] + ul.bz*wli[IVZ])*bxdpl;
        fl.mu = ul.mu*wli[IVX];
        fl.by = ul.by*wli[IVX] - bxi*wli[IVY];
        fl.bz = ul.bz*wli[IVX] - bxi*wli[IVZ];
        
        fr.d  = ur.mx;
        fr.mx = ur.mx*wri[IVX] + wri[IPP] + 0.5*bsqr - bxdpr*bxi;
        fr.my = ur.my*wri[IVX] - ur.by*bxdpr;
        fr.mz = ur.mz*wri[IVX] - ur.bz*bxdpr;
        fr.e  = (ur.e + wri[IPP] + 0.5*bsqr)*wri[IVX] -
                  (bxi*wri[IVX] + ur.by*wri[IVY] + ur.bz*wri[IVZ])*bxdpr;
        fr.mu = ur.mu*wri[IVX];
        fr.by = ur.by*wri[IVX] - bxi*wri[IVY];
        fr.bz = ur.bz*wri[IVX] - bxi*wri[IVZ];
        
        //--- Step 4.  Compute hll fluxes for all variables,
        
        // inverse of difference between right and left signal speeds
        Real idspd = 1.0/(spd[4]-spd[0]);
        
        // rho component of U^{hll} from Mignone eqn. (15); uses F_L and F_R from eqn. (6)
        Real dhll = (spd[4]*ur.d - spd[0]*ul.d - fr.d + fl.d)*idspd;
        dhll = std::max(dhll, dfloor);
        Real sqrtdhll = std::sqrt(dhll);
        
        // rho, mx, e, and mu components of F^{hll} like from Mignone eqn. (17)
        fhll.d  = (spd[4]*fl.d  - spd[0]*fr.d  + spd[4]*spd[0]*(ur.d -ul.d ))*idspd;
        fhll.mx = (spd[4]*fl.mx - spd[0]*fr.mx + spd[4]*spd[0]*(ur.mx-ul.mx))*idspd;
        fhll.e  = (spd[4]*fl.e  - spd[0]*fr.e  + spd[4]*spd[0]*(ur.e -ul.e ))*idspd;
        fhll.mu = (spd[4]*fl.mu - spd[0]*fr.mu + spd[4]*spd[0]*(ur.mu-ul.mu))*idspd;
        fhll.my = (spd[4]*fl.my - spd[0]*fr.my + spd[4]*spd[0]*(ur.my-ul.my))*idspd;
        fhll.mz = (spd[4]*fl.mz - spd[0]*fr.mz + spd[4]*spd[0]*(ur.mz-ul.mz))*idspd;
        fhll.by = (spd[4]*fl.by - spd[0]*fr.by + spd[4]*spd[0]*(ur.by-ul.by))*idspd;
        fhll.bz = (spd[4]*fl.bz - spd[0]*fr.bz + spd[4]*spd[0]*(ur.bz-ul.bz))*idspd;
        
        // ustar from paragraph between eqns. (23) and (24)
        Real ustar = fhll.d/dhll;
        
        // other components of U^{hll} like Mignone eqn. (15)
        Real mxhll = (spd[4]*ur.mx - spd[0]*ul.mx - fr.mx + fl.mx)*idspd;
        Real ehll  = (spd[4]*ur.e  - spd[0]*ul.e  - fr.e  + fl.e )*idspd;
        Real muhll = (spd[4]*ur.mu - spd[0]*ul.mu - fr.mu + fl.mu)*idspd;
        
        // Compute fhstar = 1 + Dpstar/B^2_star from hll averages
        Real myhll = (spd[4]*ur.my - spd[0]*ul.my - fr.my + fl.my)*idspd;
        Real mzhll = (spd[4]*ur.mz - spd[0]*ul.mz - fr.mz + fl.mz)*idspd;
        Real byhll = (spd[4]*ur.by - spd[0]*ul.by - fr.by + fl.by)*idspd;
        Real bzhll = (spd[4]*ur.bz - spd[0]*ul.bz - fr.bz + fl.bz)*idspd;
        Real bsqhll = bxsq + byhll*byhll + bzhll*bzhll;
        Real dphll = 3.*muhll*std::sqrt(bsqhll) - 2.0*ehll + bsqhll +
                        (fhll.d*ustar + myhll*myhll/dhll + mzhll*mzhll/dhll);
        Real fhstar = 1. + dphll / bsqhll; // Initial guess
//          // another option is to compute from mx flux (effectively like isothermal version)
//          Real fhstar = (-fhll.mx + fhll.d*ustar + pphll + 0.5*bsqhll)/bxsq;

        //--- Step 4.9  Compute S*_L and S*_R, Alfven discontinuity speeds
        Real sqrt_fhstar = fhstar > 0.0 ? std::sqrt(fhstar) : 0.0;
        spd[1] = ustar - fabs(bxi)*sqrt_fhstar/sqrtdhll;
        spd[3] = ustar + fabs(bxi)*sqrt_fhstar/sqrtdhll;
        
        //--- Step 5. Compute left and right star states for my,mz,by,bz
        Real csmean = 0.5*(cfl+cfr);
        Real tmp = (spd[0]-spd[1])*(spd[0]-spd[3]);
        if (fabs(spd[0]-spd[1]) < (SMALL_NUMBER)*csmean) {
          // degenerate case described below eqn. (39)
          ulst.my = ul.my;
          ulst.mz = ul.mz;
          ulst.by = ul.by;
          ulst.bz = ul.bz;
        } else {
          Real mfact = bxi*(ustar*fhl - wli[IVX]*fhstar + spd[0]*(fhstar-fhl))/tmp;
          Real bfact = (ul.d*SQR(spd[0]-wli[IVX]) - bxsq*fhl)/(dhll*tmp);
          
          ulst.my = dhll*wli[IVY] - ul.by*mfact; // eqn. (30) of Mignone
          ulst.mz = dhll*wli[IVZ] - ul.bz*mfact; // eqn. (31) of Mignone
          ulst.by = ul.by*bfact; // eqn. (32) of Mignone
          ulst.bz = ul.bz*bfact; // eqn. (33) of Mignone
        }
        
        tmp = (spd[4]-spd[1])*(spd[4]-spd[3]);
        if (fabs(spd[4]-spd[3]) < (SMALL_NUMBER)*csmean) {
          // degenerate case described below eqn. (39)
          urst.my = ur.my;
          urst.mz = ur.mz;
          urst.by = ur.by;
          urst.bz = ur.bz;
        } else {
          Real mfact = bxi*(ustar*fhr - wri[IVX]*fhstar + spd[4]*(fhstar-fhr))/tmp;
          Real bfact = (ur.d*SQR(spd[4]-wri[IVX]) - bxsq*fhr)/(dhll*tmp);
          
          urst.my = dhll*wri[IVY] - ur.by*mfact; // eqn. (30) of Mignone
          urst.mz = dhll*wri[IVZ] - ur.bz*mfact; // eqn. (31) of Mignone
          urst.by = ur.by*bfact; // eqn. (32) of Mignone
          urst.bz = ur.bz*bfact; // eqn. (33) of Mignone
        }
        
        //--- Step 6. Compute compute fluxes in L,R star regions.
        // We do this even if we want the central flux, because FsL and FsR used to
        // compute Fsc
        
        //--- Step 6.1 Compute "Anti-Diffusive Control" (Simon & Mandell 2018).
        // Our solver is more HLLE-like if ∆p/B2 becomes oscillatory, or ∆p/B^2~-1.
        Real bsqltmp, bsqrtmp, ddpm1, ddpp1;
        if (ivx==IVX) {
          // Compute nearby variation of Dp/B^2
          bsqltmp = SQR(bx(k,j,i-1)) + SQR(wl(IBY,k,j,i-1)) + SQR(wl(IBZ,k,j,i-1));
          bsqrtmp = SQR(bx(k,j,i-1)) + SQR(wr(IBY,k,j,i-1)) + SQR(wr(IBZ,k,j,i-1));
          ddpm1 = (wl(IPP,k,j,i-1) - wl(IPR,k,j,i-1))/bsqltmp -
                  (wr(IPP,k,j,i-1) - wr(IPR,k,j,i-1))/bsqrtmp;
          
          bsqltmp = SQR(bx(k,j,i+1)) + SQR(wl(IBY,k,j,i+1)) + SQR(wl(IBZ,k,j,i+1));
          bsqrtmp = SQR(bx(k,j,i+1)) + SQR(wr(IBY,k,j,i+1)) + SQR(wr(IBZ,k,j,i+1));
          ddpp1 = (wl(IPP,k,j,i+1) - wl(IPR,k,j,i+1))/bsqltmp -
                  (wr(IPP,k,j,i+1) - wr(IPR,k,j,i+1))/bsqrtmp;
        } else if (ivx==IVY) {
          // Compute nearby variation of Dp/B^2
          bsqltmp = SQR(bx(k,j-1,i)) + SQR(wl(IBY,k,j-1,i)) + SQR(wl(IBZ,k,j-1,i));
          bsqrtmp = SQR(bx(k,j-1,i)) + SQR(wr(IBY,k,j-1,i)) + SQR(wr(IBZ,k,j-1,i));
          ddpm1 = (wl(IPP,k,j-1,i) - wl(IPR,k,j-1,i))/bsqltmp -
                  (wr(IPP,k,j-1,i) - wr(IPR,k,j-1,i))/bsqrtmp;
          
          bsqltmp = SQR(bx(k,j+1,i)) + SQR(wl(IBY,k,j+1,i)) + SQR(wl(IBZ,k,j+1,i));
          bsqrtmp = SQR(bx(k,j+1,i)) + SQR(wr(IBY,k,j+1,i)) + SQR(wr(IBZ,k,j+1,i));
          ddpp1 = (wl(IPP,k,j+1,i) - wl(IPR,k,j+1,i))/bsqltmp -
                  (wr(IPP,k,j+1,i) - wr(IPR,k,j+1,i))/bsqrtmp;
        } else if (ivx==IVZ) {
          // Compute nearby variation of Dp/B^2
          bsqltmp = SQR(bx(k-1,j,i)) + SQR(wl(IBY,k-1,j,i)) + SQR(wl(IBZ,k-1,j,i));
          bsqrtmp = SQR(bx(k-1,j,i)) + SQR(wr(IBY,k-1,j,i)) + SQR(wr(IBZ,k-1,j,i));
          ddpm1 = (wl(IPP,k-1,j,i) - wl(IPR,k-1,j,i))/bsqltmp -
                  (wr(IPP,k-1,j,i) - wr(IPR,k-1,j,i))/bsqrtmp;
          
          bsqltmp = SQR(bx(k+1,j,i)) + SQR(wl(IBY,k+1,j,i)) + SQR(wl(IBZ,k+1,j,i));
          bsqrtmp = SQR(bx(k+1,j,i)) + SQR(wr(IBY,k+1,j,i)) + SQR(wr(IBZ,k+1,j,i));
          ddpp1 = (wl(IPP,k+1,j,i) - wl(IPR,k+1,j,i))/bsqltmp -
                  (wr(IPP,k+1,j,i) - wr(IPR,k+1,j,i))/bsqrtmp;
        }
        // Switch to HLLE if fh becomes oscillatory
        Real ddp = std::max(fabs(fhl-fhr), fabs(ddpm1));
        ddp = std::max(ddp, fabs(ddpp1));
        // Switch to HLLE for fh close to zero
        Real smallfh_limit = (fhstar < fh_floor) ? 0. : fhstar;
        // If omega~0, it's HLLE, if omega~1 it's HLLD
        Real omega = std::exp(-anti_diff_alpha*ddp) * smallfh_limit;
        
        //--- Step 6.2 Compute star region fluxes of my,mz,by,bz
        // Left-hand star region fluxes
        Real flstmy = fhll.my + omega*spd[0]*(ulst.my - myhll);
        Real flstmz = fhll.mz + omega*spd[0]*(ulst.mz - mzhll);
        Real flstby = fhll.by + omega*spd[0]*(ulst.by - byhll);
        Real flstbz = fhll.bz + omega*spd[0]*(ulst.bz - bzhll);
        // Right-hand star region fluxes
        Real frstmy = fhll.my + omega*spd[4]*(urst.my - myhll);
        Real frstmz = fhll.mz + omega*spd[4]*(urst.mz - mzhll);
        Real frstby = fhll.by + omega*spd[4]*(urst.by - byhll);
        Real frstbz = fhll.bz + omega*spd[4]*(urst.bz - bzhll);
        
        //--- Step 7.  Compute flux
          
        if (spd[0] >= 0.0) {
          // return Fl if flow is supersonic, eqn. (38a) of Mignone
          flxi[IDN] = fl.d;
          flxi[IVX] = fl.mx;
          flxi[IVY] = fl.my;
          flxi[IVZ] = fl.mz;
          flxi[IPR] = fl.e;
          flxi[IPP] = fl.mu;
          flxi[IBY] = fl.by;
          flxi[IBZ] = fl.bz;
        } else if (spd[4] <= 0.0) {
          // return Fr if flow is supersonic, eqn. (38e) of Mignone
          flxi[IDN] = fr.d;
          flxi[IVX] = fr.mx;
          flxi[IVY] = fr.my;
          flxi[IVZ] = fr.mz;
          flxi[IPR] = fr.e;
          flxi[IPP] = fr.mu;
          flxi[IBY] = fr.by;
          flxi[IBZ] = fr.bz;
        } else if (spd[1] >= 0.0) {
          // return (Fl+Sl*(Ulst-Ul))=Fhll+SL(ULst-Uhll), eqn. (38b) of Mignone
          flxi[IDN] = fhll.d;
          flxi[IVX] = fhll.mx;
          flxi[IVY] = flstmy;
          flxi[IVZ] = flstmz;
          flxi[IPR] = fhll.e ;
          flxi[IPP] = fhll.mu;
          flxi[IBY] = flstby;
          flxi[IBZ] = flstbz;
        } else if (spd[3] <= 0.0) {
          // return (Fr+Sr*(Urst-Ur))=Fhll+SR*(URst-Uhll), eqn. (38d) of Mignone
          flxi[IDN] = fhll.d;
          flxi[IVX] = fhll.mx;
          flxi[IVY] = frstmy;
          flxi[IVZ] = frstmz;
          flxi[IPR] = fhll.e ;
          flxi[IPP] = fhll.mu;
          flxi[IBY] = frstby;
          flxi[IBZ] = frstbz;
        } else {
          Real isrsl = 1.0/(spd[3]-spd[1]);
          Real srtsl = spd[3]*spd[1];
          // return Fcst, but compute directly from consistency, eqn. (16) of Mignone
          flxi[IDN] = fhll.d;
          flxi[IVX] = fhll.mx;
          flxi[IVY] = (spd[3]*flstmy - spd[1]*frstmy + srtsl*(urst.my-ulst.my) )*isrsl;
          flxi[IVZ] = (spd[3]*flstmz - spd[1]*frstmz + srtsl*(urst.mz-ulst.mz) )*isrsl;
          flxi[IPR] = fhll.e;
          flxi[IPP] = fhll.mu;
          flxi[IBY] = (spd[3]*flstby - spd[1]*frstby + srtsl*(urst.by-ulst.by) )*isrsl;
          flxi[IBZ] = (spd[3]*flstbz - spd[1]*frstbz + srtsl*(urst.bz-ulst.bz) )*isrsl;
        }

        flx(IDN,k,j,i) = flxi[IDN];
        flx(ivx,k,j,i) = flxi[IVX];
        flx(ivy,k,j,i) = flxi[IVY];
        flx(ivz,k,j,i) = flxi[IVZ];
        flx(IPR,k,j,i) = flxi[IPR];
        flx(IPP,k,j,i) = flxi[IPP];
        ey(k,j,i) = -flxi[IBY];
        ez(k,j,i) =  flxi[IBZ];

      }
    }}
//  std::cout <<"\n\n";
  
  return;
}
