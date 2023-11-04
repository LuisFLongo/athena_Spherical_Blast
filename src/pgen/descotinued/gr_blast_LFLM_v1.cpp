//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_blast.cpp
//! \brief Problem generator for GRMHD spherical blast wave in flat spacetime.


// Comment : This version seats on top of gr_blast.cpp as out of the box
// 	     the generalization here is to create the blast wave in the BoyerLindquist coordinates
// 	     and including a central compact object with mass M and spin a.
// 	     Further versions will include the time dependent boundary condition in order to include ejecta evolution.

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <cstring>    // strcmp()

// Athena++ headers
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"          // ParameterInput

// Configuration checking
//#if not GENERAL_RELATIVITY
//#error "This problem generator must be used with general relativity"
//#endif

// Declarations
namespace {
//void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
//                             Real *px, Real *py, Real *pz);
void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                  Real *ptheta, Real *pphi);


void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3);
Real DistanceBetweenPoints(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3);
Real m, a;



void InnerBoundary(MeshBlock *pmb, Coordinates *pcoord,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku,
                  int ngh);

} // namespace


void Mesh::InitUserMeshData(ParameterInput *pin) {


  EnrollUserBoundaryFunction(BoundaryFace::inner_x1,InnerBoundary);


}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }


  // Get mass and spin of black hole
  m = pcoord->GetMass();
  a = pcoord->GetSpin();

  // Get ratio of specific heats
  // Real gamma_adi = peos->GetGamma();
  // Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read problem parameters
  Real num_x = pin->GetReal("problem", "num_x");
  Real num_y = pin->GetReal("problem", "num_y");
  Real x_spacing = pin->GetReal("problem", "x_spacing");
  Real y_spacing = pin->GetReal("problem", "y_spacing");
  Real radius = pin->GetReal("problem", "radius");
  Real rho_inner = pin->GetReal("problem", "rho_inner");
  Real pgas_inner = pin->GetReal("problem", "pgas_inner");
  Real rho_outer = pin->GetReal("problem", "rho_outer");
  Real pgas_outer = pin->GetReal("problem", "pgas_outer");
  Real bx = 0.0, by = 0.0, bz = 0.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    bx = pin->GetReal("problem", "bx");
    by = pin->GetReal("problem", "by");
    bz = pin->GetReal("problem", "bz");
  }

  // Prepare auxiliary arrays
  //AthenaArray<Real> b, g, gi;
  //b.NewAthenaArray(3, ncells3, ncells2, ncells1);
  AthenaArray<Real> g,gi ;
  g.NewAthenaArray(NMETRIC,iu+1);// ncells1);
  gi.NewAthenaArray(NMETRIC,iu+1);// ncells1);

  // Initialize hydro variables
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i=il; i<=iu; ++i) {
        Real r(0.0), theta(0.0), phi(0.0);


        GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j), pcoord->x3v(k), &r,
                                     &theta, &phi);
        Real rho, pgas, ut, ur;


        // Set pressure and density
        
        if (r < radius) {
          rho = rho_inner;
          pgas = pgas_inner;
        } else {
          rho = rho_outer;
          pgas = pgas_outer;
        }

        Real u0(0.0), u1(0.0), u2(0.0), u3(0.0);
        TransformVector(ut, ur, 0.0, 0.0, r, theta, phi, &u0, &u1, &u2, &u3);
        Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
        Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
        Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IM1,k,j,i) = phydro->w1(IM1,k,j,i) = uu1;
        phydro->w(IM2,k,j,i) = phydro->w1(IM2,k,j,i) = uu2;
        phydro->w(IM3,k,j,i) = phydro->w1(IM3,k,j,i) = uu3;

        // Get Minkowski coordinates of point
      }
    }
  }


  // Initialize magnetic field
  //if (MAGNETIC_FIELDS_ENABLED) {
    // Find normalization
    //Real r, theta, phi;
    //GetBoyerLindquistCoordinates(pcoord->x1f(is), pcoord->x2v((jl+ju)/2),
    //                             pcoord->x3v((kl+ku)/2), &r, &theta, &phi);
    //Real rho, pgas, ut, ur;
    //CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
    //Real bbr = 1.0/SQR(r);
    //Real bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
    //Real br = (bbr + bt * ur) / ut;
    //Real bsq = -(1.0-2.0*m/r) * SQR(bt) + 1.0/(1.0-2.0*m/r) * SQR(br);
    //Real bsq_over_rho_actual = bsq/rho;
    //Real normalization = std::sqrt(bsq_over_rho/bsq_over_rho_actual);

    // Set face-centered field
    //for (int k=kl; k<=ku+1; ++k) {
    //  for (int j=jl; j<=ju+1; ++j) {
    //    for (int i=il; i<=iu+1; ++i) {
          // Set B^1
    //      if (j != ju+1 && k != ku+1) {
    //        GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3v(k),
    //                                     &r, &theta, &phi);
    //        CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
    //        bbr = normalization/SQR(r);
    //        bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
    //        br = (bbr + bt * ur) / ut;
    //        Real u0, u1, u2, u3;
    //        TransformVector(ut, ur, 0.0, 0.0, r, theta, phi, &u0, &u1, &u2, &u3);
    //        Real b0, b1, b2, b3;
    //        TransformVector(bt, br, 0.0, 0.0, r, theta, phi, &b0, &b1, &b2, &b3);
    //        pfield->b.x1f(k,j,i) = b1 * u0 - b0 * u1;
    //      }

          // Set B^2
    //      if (i != iu+1 && k != ku+1) {
    //        GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3v(k),
    //                                     &r, &theta, &phi);
    //        CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
    //        bbr = normalization/SQR(r);
    //        bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
    //        br = (bbr + bt * ur) / ut;
    //        Real u0, u1, u2, u3;
    //        TransformVector(ut, ur, 0.0, 0.0, r, theta, phi, &u0, &u1, &u2, &u3);
    //        Real b0, b1, b2, b3;
    //        TransformVector(bt, br, 0.0, 0.0, r, theta, phi, &b0, &b1, &b2, &b3);
    //        pfield->b.x2f(k,j,i) = b2 * u0 - b0 * u2;
    //      }

          // Set B^3
    //      if (i != iu+1 && j != ju+1) {
    //        GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j), pcoord->x3f(k),
    //                                     &r, &theta, &phi);
    //        CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
    //        bbr = normalization/SQR(r);
    //        bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
    //        br = (bbr + bt * ur) / ut;
    //        Real u0, u1, u2, u3;
    //        TransformVector(ut, ur, 0.0, 0.0, r, theta, phi, &u0, &u1, &u2, &u3);
    //        Real b0, b1, b2, b3;
    //        TransformVector(bt, br, 0.0, 0.0, r, theta, phi, &b0, &b1, &b2, &b3);
    //        pfield->b.x3f(k,j,i) = b3 * u0 - b0 * u3;
    //      }
    //    }
    //  }
    //}

    // Calculate cell-centered magnetic field
    //pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, il, iu, jl, ju, kl,
    //                                   ku);
 // }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Delete auxiliary array
  // Initialize magnetic field

  return;
         
           
}


void Mesh::InitUserMeshData(ParameterInput *pin) {                                                                                                          

     EnrollUserBoundaryFunction(BoundaryFace::inner_x1,InnerBoundary);
     return; 
}     

namespace {
//----------------------------------------------------------------------------------------
// Function for returning corresponding Minkowski coordinates of point
// Inputs:
//   x0,x1,x2,x3: global coordinates to be converted
// Outputs:
//   pt,px,py,pz: variables pointed to set to Minkowski coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

//void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
//                             Real *px, Real *py, Real *pz) {
//  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
//    *pt = x0;
//    *px = x1;
//    *py = x2;
//    *pz = x3;
//  }
//  return;
//}
//

void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                  Real *ptheta, Real *pphi) {
  if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0 ||
      std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    *pr = x1;
    *ptheta = x2;
    *pphi = x3;
  }
  return;
}



//----------------------------------------------------------------------------------------
// Function for transforming 4-vector from Minkowski to desired coordinates
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   x,y,z: Minkowski coordinates of point
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in desired coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

void TransformVector(Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, Real r,
                     Real theta, Real phi,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
    *pa0 = a0_bl;
    *pa1 = a1_bl;
    *pa2 = a2_bl;
    *pa3 = a3_bl;
  } else if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    Real delta = SQR(r) - 2.0*m*r + SQR(a);
    *pa0 = a0_bl + 2.0*m*r/delta * a1_bl;
    *pa1 = a1_bl;
    *pa2 = a2_bl;
    *pa3 = a3_bl + a/delta * a1_bl;
  }
  return;
}


//----------------------------------------------------------------------------------------
// Function for returning spatial separation between points at same time
// Inputs:
//   x1,x2,x3: spatial coordinates of one point
//   y1,y2,y3: spatial coordinates of other point
// Outputs:
//   returned value: spatial separation between x and y
// Notes:
//   distance function is Euclidean in Minkowski coordinates

Real DistanceBetweenPoints(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3) {
  Real distance = 0.0;
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    distance = std::sqrt(SQR(x1-y1) + SQR(x2-y2) + SQR(x3-y3));
  }
  return distance;
}

void InnerBoundary(MeshBlock *pmb, Coordinates *pcoord,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku,
                  int ngh) {


        printf("Entering BC condition \n");

        for (int k=kl; k<=ku; ++k) {
                for (int j=jl; j<=ju; ++j){
                        for (int i=1; i<=ngh; ++i){

                                Real rval   = pcoord->x1v(il-i) ;
                                Real thval  = pcoord->x2v(j) ;
                                Real phival = pcoord->x3v(k) ;

                                prim(IDN, k , j, il-i) = 10000 ;
                                prim(IPR, k , j, il-i) = 25000 ;
                                prim(IVX, k , j, il-i) = 0.3 ;
                                prim(IVY, k , j, il-i) = 0.0 ;
                                prim(IVZ, k , j, il-i) = 0.0 ;

				cons(IM1, k , j, il-i) = 0.0 ;
				cons(IM2, k , j, il-i) = 0.0 ;
				
                        }
                }
        }
//	printf("TESTING BC 1 \n");
        pmb->peos->PrimitiveToConserved(prim, pmb->pfield->bcc , cons, pcoord, il, iu, jl, ju, kl, ku);
        return;
}


} // namespace
