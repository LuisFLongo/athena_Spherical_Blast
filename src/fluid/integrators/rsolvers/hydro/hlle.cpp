//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 

//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../../../athena.hpp"         // enums, macros, Real
#include "../../../../athena_arrays.hpp"  // AthenaArray
#include "../../../fluid.hpp"             // Fluid
#include "../../../eos/eos.hpp"           // GetGamma

//======================================================================================
//! \file hlle.cpp
//  \brief HLLE Riemann solver for hydrodynamics
//
//  Computes 1D fluxes using the Harten-Lax-van Leer (HLL) Riemann solver.  This flux is
//  very diffusive, especially for contacts, and so it is not recommended for use in
//  applications.  However, as shown by Einfeldt et al.(1991), it is positively
//  conservative (cannot return negative densities or pressure), so it is a useful
//  option when other approximate solvers fail and/or when extra dissipation is needed.
//
// REFERENCES:
// - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//   Springer-Verlag, Berlin, (1999) chpt. 10.
// - Einfeldt et al., "On Godunov-type methods near low densities", JCP, 92, 273 (1991)
// - A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and Godunov-type
//   schemes for hyperbolic conservation laws", SIAM Review 25, 35-61 (1983).
//======================================================================================

void FluidIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[NFLUID],wri[NFLUID],wroe[NFLUID],fl[NFLUID],fr[NFLUID],flxi[NFLUID];
  Real gamma_m1 = pmy_fluid->pf_eos->GetGamma() - 1.0;

#pragma simd
  for (int i=il; i<=iu; ++i){

//--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    if (NON_BAROTROPIC_EOS) wli[IEN]=wl(IEN,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    if (NON_BAROTROPIC_EOS) wri[IEN]=wr(IEN,i);

//--- Step2.  Compute Roe-averaged state

    Real sqrtdl = sqrt(wli[IDN]);
    Real sqrtdr = sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN] = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;

// Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
// rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl

    Real el,er,hroe;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IEN]/gamma_m1 + 0.5*wli[IDN]*
        (wli[IVX]*wli[IVX] + wli[IVY]*wli[IVY] + wli[IVZ]*wli[IVZ]);
      er = wri[IEN]/gamma_m1 + 0.5*wri[IDN]*
        (wri[IVX]*wri[IVX] + wri[IVY]*wri[IVY] + wri[IVZ]*wri[IVZ]);
      hroe = ((el + wli[IEN])/sqrtdl + (er + wri[IEN])/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute sound speed in L,R, and Roe-averaged states

    Real cl = pmy_fluid->pf_eos->SoundSpeed(wli);
    Real cr = pmy_fluid->pf_eos->SoundSpeed(wri);
    Real a;
    if (NON_BAROTROPIC_EOS) {
      Real q = hroe - 0.5*(wroe[IVX]*wroe[IVX]+wroe[IVY]*wroe[IVY]+wroe[IVZ]*wroe[IVZ]);
      if (q < 0.0) q=0.0;
      a = sqrt(gamma_m1*q);
    } else {
      a = pmy_fluid->pf_eos->SoundSpeed(wroe);
    }

//--- Step 4. Compute the max/min wave speeds based on L/R and Roe-averaged values

    Real al = std::min((wroe[IVX] - a),(wli[IVX] - cl));
    Real ar = std::max((wroe[IVX] + a),(wri[IVX] + cr));

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

//-- Step 5. Compute L/R fluxes along the lines bm/bp: F_L - (S_L)U_L; F_R - (S_R)U_R

    fl[IDN] = wli[IDN]*wli[IVX] - bm*wli[IDN];
    fr[IDN] = wri[IDN]*wri[IVX] - bp*wri[IDN];

    fl[IVX] = wli[IDN]*wli[IVX]*(wli[IVX] - bm);
    fr[IVX] = wri[IDN]*wri[IVX]*(wri[IVX] - bp);

    fl[IVY] = wli[IDN]*wli[IVY]*(wli[IVX] - bm);
    fr[IVY] = wri[IDN]*wri[IVY]*(wri[IVX] - bp);

    fl[IVZ] = wli[IDN]*wli[IVZ]*(wli[IVX] - bm);
    fr[IVZ] = wri[IDN]*wri[IVZ]*(wri[IVX] - bp);

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IEN];
      fr[IVX] += wri[IEN];
      fl[IEN] = el*(wli[IVX] - bm) + wli[IEN]*wli[IVX];
      fr[IEN] = er*(wri[IVX] - bp) + wri[IEN]*wri[IVX];
    } else {
      Real iso_cs = pmy_fluid->pf_eos->SoundSpeed(wli);
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

//--- Step 6. Compute the HLLE flux at interface.

    Real tmp;
    if (bp == bm) {
      tmp = 0.0;
    } else {
      tmp = 0.5*(bp + bm)/(bp - bm);
    }

    flxi[IDN] = 0.5*(fl[IDN]+fr[IDN]) + (fl[IDN]-fr[IDN])*tmp;
    flxi[IVX] = 0.5*(fl[IVX]+fr[IVX]) + (fl[IVX]-fr[IVX])*tmp;
    flxi[IVY] = 0.5*(fl[IVY]+fr[IVY]) + (fl[IVY]-fr[IVY])*tmp;
    flxi[IVZ] = 0.5*(fl[IVZ]+fr[IVZ]) + (fl[IVZ]-fr[IVZ])*tmp;
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN]+fr[IEN]) + (fl[IEN]-fr[IEN])*tmp;

    flx(IDN,i) = flxi[IDN];
    flx(ivx,i) = flxi[IVX];
    flx(ivy,i) = flxi[IVY];
    flx(ivz,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,i) = flxi[IEN];
  }

  return;
}