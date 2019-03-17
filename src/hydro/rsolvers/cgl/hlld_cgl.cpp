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
// minimizes changes required to adopt athena4.2 version of this solver
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
  
  // Number of Newton-Raphson iterations to solve for (pp-pl)/B^2
  const int num_newton_iter = 4;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd private(flxi,wli,wri,spd)
      for (int i=il; i<=iu; ++i) {
        Cons1D ul,ur;                   // L/R states, conserved variables (computed)
        Cons1D ulst,urst,ucst;          // Conserved variable for all states
        Cons1D fl,fr;                   // Fluxes for left & right states
        
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
        bmagl = bmagl > bmag_floor ? bmagl : bmag_floor;
        bmagr = bmagr > bmag_floor ? bmagr : bmag_floor;
        
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
        
        //--- Step 4.  Compute hll averages
        
        // inverse of difference between right and left signal speeds
        Real idspd = 1.0/(spd[4]-spd[0]);
        
        // rho component of U^{hll} from Mignone eqn. (15); uses F_L and F_R from eqn. (6)
        Real dhll = (spd[4]*ur.d - spd[0]*ul.d - fr.d + fl.d)*idspd;
        dhll = std::max(dhll, dfloor);
        Real sqrtdhll = std::sqrt(dhll);
        
        // rho, mx, e, and mu components of F^{hll} like from Mignone eqn. (17)
        Real fdhll  = (spd[4]*fl.d  - spd[0]*fr.d  + spd[4]*spd[0]*(ur.d -ul.d ))*idspd;
        Real fmxhll = (spd[4]*fl.mx - spd[0]*fr.mx + spd[4]*spd[0]*(ur.mx-ul.mx))*idspd;
        Real fehll  = (spd[4]*fl.e  - spd[0]*fr.e  + spd[4]*spd[0]*(ur.e -ul.e ))*idspd;
        Real fmuhll = (spd[4]*fl.mu - spd[0]*fr.mu + spd[4]*spd[0]*(ur.mu-ul.mu))*idspd;
        
        // ustar from paragraph between eqns. (23) and (24)
        Real ustar = fdhll/dhll;
        
        // other components of U^{hll} like Mignone eqn. (15)
        Real mxhll = (spd[4]*ur.mx - spd[0]*ul.mx - fr.mx + fl.mx)*idspd;
        Real ehll  = (spd[4]*ur.e  - spd[0]*ul.e  - fr.e  + fl.e )*idspd;
        Real muhll = (spd[4]*ur.mu - spd[0]*ul.mu - fr.mu + fl.mu)*idspd;
        
        //--- Step 4.5  Compute dpobsq = 1 + Dp/B^2 in the star region. This is assumed
        // to be constant across the central region.
        // Compute by solving nonlinear equation from fmxhll, using B^2=sqrt(B*L^2*B*R^2)
        // from Mignone (32)-(33)
        
        // Only run nonlinear solve if left and right states very different
        Real fhstar = 0.5*(fhr+fhl); // Initial guess (add handling of small dpobsq)
        if ( fabs(fhr - fhl) > 1e10){
          // Start by computing unchanging auxilliary variables for the iteration
          Real slul = spd[0] - wli[IVX], srur = spd[4] - wri[IVX];
          Real slus = spd[0] - ustar, srus = spd[4] - ustar;
          Real al = dhll*SQR(slus) / bxsq;
          Real ar = dhll*SQR(srus) / bxsq;
          Real a_full = std::sqrt((bsql - bxsq)*(bsqr - bxsq))/( SQR(bxsq) ) *
                  fabs( SQR(slul)*ul.d - bxsq*fhl) * (SQR(srur)*ur.d - bxsq*fhr );
          Real frhou_fmx_bx = fdhll*ustar - fmxhll + 0.5*bxsq;
          // Newton Raphson iteration for 1 + Dp/B^2
//          std::cout << "NR init L=" << fhl << ": NR init R=" << fhr <<"\n";
          for (int n = 0; n < num_newton_iter; ++n) {
            Real aldpardp = fabs( (al - fhstar) * (ar - fhstar) );
            Real aral_aldpardp = 0.5 * (ar + al - 2*fhstar) / SQR(aldpardp);
            Real bmag_guess = std::sqrt( bxsq + a_full/aldpardp );
            // Function to minimize, divided by its derivative
            Real f0 = (frhou_fmx_bx + muhll*bmag_guess - bxsq*fhstar + 0.5*a_full/aldpardp);
            Real delta_dp =  -f0 /
                  (-bxsq + 0.5*a_full*bxsq*aral_aldpardp + a_full*muhll*aral_aldpardp/bmag_guess);
            fhstar += delta_dp;
//            std::cout << n << " f0="<< f0 <<", dp=" << delta_dp << ": ";
          }
//          std::cout << "\nNR final " << fhstar << "\n";
        }
        //--- Step 4.9  Compute S*_L and S*_R, Alfven discontinuity speeds
        Real sqrt_dpstar = std::sqrt(fhstar);
        spd[1] = ustar - fabs(bxi)*sqrt_dpstar/sqrtdhll;
        spd[3] = ustar + fabs(bxi)*sqrt_dpstar/sqrtdhll;
        
        //--- Step 5. Compute intermediate states
        
        // Ul* - eqn. (20) of Mignone
        ulst.d  = dhll;
        ulst.mx = mxhll; // eqn. (24) of Mignone
        ulst.e  = ehll;
        ulst.mu = muhll;
        
        Real tmp = (spd[0]-spd[1])*(spd[0]-spd[3]);
        if (fabs(spd[0]-spd[1]) < (SMALL_NUMBER)*0.0) { //FIX THIS
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
        
        // Ur* - eqn. (20) of Mignone */
        urst.d  = dhll;
        urst.mx = mxhll; // eqn. (24) of Mignone
        urst.e  = ehll;
        urst.mu = muhll;
        
        tmp = (spd[4]-spd[1])*(spd[4]-spd[3]);
        if (fabs(spd[4]-spd[3]) < (SMALL_NUMBER)*0.0) { //FIX THIS
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
        
        // Uc*
        Real x = sqrt_dpstar*sqrtdhll*(bxi > 0.0 ? 1.0 : -1.0); // from below eqn. (37) of Mignone
        ucst.d  = dhll;  // eqn. (20) of Mignone
        ucst.mx = mxhll; // eqn. (24) of Mignone
        ucst.e  = ehll;
        ucst.mu = muhll;
        ucst.my = 0.5*(ulst.my + urst.my + (urst.by-ulst.by)*x); // eqn. (34) of Mignone
        ucst.mz = 0.5*(ulst.mz + urst.mz + (urst.bz-ulst.bz)*x); // eqn. (35) of Mignone
        ucst.by = 0.5*(ulst.by + urst.by + (urst.my-ulst.my)/x); // eqn. (36) of Mignone
        ucst.bz = 0.5*(ulst.bz + urst.bz + (urst.mz-ulst.mz)/x); // eqn. (37) of Mignone
        
//        // FOR TESTING
//        Real fmyhll = (spd[4]*fl.my - spd[0]*fr.my + spd[4]*spd[0]*(ur.my-ul.my))*idspd;
//        Real fmzhll = (spd[4]*fl.mz - spd[0]*fr.mz + spd[4]*spd[0]*(ur.mz-ul.mz))*idspd;
//        Real fbyhll = (spd[4]*fl.by - spd[0]*fr.by + spd[4]*spd[0]*(ur.by-ul.by))*idspd;
//        Real fbzhll = (spd[4]*fl.bz - spd[0]*fr.bz + spd[4]*spd[0]*(ur.bz-ul.bz))*idspd;
        
        //--- Step 6.  Compute flux
        
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
          // return (Fl+Sl*(Ulst-Ul)), eqn. (38b) of Mignone
          flxi[IDN] = fl.d  + spd[0]*(ulst.d  - ul.d);
          flxi[IVX] = fl.mx + spd[0]*(ulst.mx - ul.mx);
          flxi[IVY] = fl.my + spd[0]*(ulst.my - ul.my);
          flxi[IVZ] = fl.mz + spd[0]*(ulst.mz - ul.mz);
          flxi[IPR] = fl.e  + spd[0]*(ulst.e  - ul.e);
          flxi[IPP] = fr.mu + spd[0]*(ulst.mu - ul.mu);
          flxi[IBY] = fl.by + spd[0]*(ulst.by - ul.by);
          flxi[IBZ] = fl.bz + spd[0]*(ulst.bz - ul.bz);
        } else if (spd[3] <= 0.0) {
          // return (Fr+Sr*(Urst-Ur)), eqn. (38d) of Mignone
          flxi[IDN] = fr.d  + spd[4]*(urst.d  - ur.d);
          flxi[IVX] = fr.mx + spd[4]*(urst.mx - ur.mx);
          flxi[IVY] = fr.my + spd[4]*(urst.my - ur.my);
          flxi[IVZ] = fr.mz + spd[4]*(urst.mz - ur.mz);
          flxi[IPR] = fr.e  + spd[0]*(urst.e  - ur.e);
          flxi[IPP] = fr.mu + spd[0]*(urst.mu - ur.mu);
          flxi[IBY] = fr.by + spd[4]*(urst.by - ur.by);
          flxi[IBZ] = fr.bz + spd[4]*(urst.bz - ur.bz);
        } else {
          // return Fcst, eqn. (38c) of Mignone, using eqn. (24)
          flxi[IDN] = dhll*ustar;
          flxi[IVX] = fmxhll;
          flxi[IVY] = ucst.my*ustar - bxi*ucst.by*fhstar;
          flxi[IVZ] = ucst.mz*ustar - bxi*ucst.bz*fhstar;
          flxi[IPR] = fehll;
          flxi[IPP] = fmuhll;
          flxi[IBY] = ucst.by*ustar - bxi*ucst.my/ucst.d;
          flxi[IBZ] = ucst.bz*ustar - bxi*ucst.mz/ucst.d;
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
  
  return;
}
