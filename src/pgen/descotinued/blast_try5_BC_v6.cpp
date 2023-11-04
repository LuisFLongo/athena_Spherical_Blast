
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//  \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//         cylindrical, and spherical coordinates.  Contains post-processing code
//         to check whether blast is spherical for regression tests
//
// REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// Comment : this version includes some corrections when moving from the cartesian coordinates 
// 	     to spherical coordinates in the boundary condition, specially in the velocity


// C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../parameter_input.hpp"
#include "../utils/interp_table3D.hpp"

int RefinementCondition(MeshBlock *pmb);

// LFLM inclusion starts

Real threshold, method, rbc, velfloor, BCThetamin, BCPhimin, BCTmin, BCThetamax, BCPhimax, BCTmax, BCMtheta, BCMphi,BCMt , BCtin , BCMghts ,BCDtheta, BCDphi, BCDt, Gammaval, Kval;

int BCNtheta, BCNphi,BCNt, BCNghts;

InterpTable3D *table;

void InnerBoundary(MeshBlock *pmb, Coordinates *pcoord,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku,
                  int ngh);

void OutflowReader(const char *filename, InterpTable3D *table);

Real DenProfile(const Real th, const Real ph, const Real tm);
Real PresProfile(const Real th, const Real ph, const Real tm);
Real VRProfile(const Real th, const Real ph, const Real tm);
Real VTProfile(const Real th, const Real ph, const Real tm);
Real VPProfile(const Real th, const Real ph, const Real tm);

// LFLM inclusion ends

void Mesh::InitUserMeshData(ParameterInput *pin) {
//LFLM inclusion start
  velfloor   = pin->GetReal("problem","vfloor");

  BCThetamin = pin->GetReal("problem","BCx2min");
  BCPhimin   = pin->GetReal("problem","BCx3min");
  BCTmin     = pin->GetReal("problem","BCtmin") ;

  BCThetamax = pin->GetReal("problem","BCx2max");
  BCPhimax   = pin->GetReal("problem","BCx3max");
  BCTmax     = pin->GetReal("problem","BCtmax") ;

  BCMtheta   = pin->GetReal("problem","BCNx2");
  BCMphi     = pin->GetReal("problem","BCNx3");
  BCMt       = pin->GetReal("problem","BCNt") ;

  BCMt       = pin->GetReal("problem","BCNt") ;
  BCtin      = pin->GetReal("problem","BCtin");

  BCMghts    = pin->GetReal("problem", "BCNghts");
  BCNghts  = BCMghts ;

  //printf("BCMtheta = %6.4lf\n", BCMtheta);
  //printf("BCMphi   = %6.4lf\n", BCMphi);
  //printf("BCMt     = %6.4lf\n", BCMt);

  BCDtheta   = (BCThetamax-BCThetamin)/(BCMtheta-1);
  BCDphi     = (BCPhimax-BCPhimin)/(BCMphi-1)      ; 
  BCDt       = (BCTmax-BCTmin)/(BCMt-1)            ; 

  Gammaval = pin->GetReal("hydro","gamma");
  Kval     = pin->GetReal("hydro","kappa");

  //printf("Test1 \n");
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1,InnerBoundary);
  //printf("Test2 \n");
  std::string fn;
  fn = pin->GetString("problem","outflow_table");
  //outflowtable = new LagrangeInterpND<2*NGHOST-1, 3>(origin, delta, size, coord); // is that right?
  table = new InterpTable3D();
  OutflowReader(fn.c_str(), table);
  if (adaptive) {
    //printf("Test3 \n");
    //EnrollUserBoundaryFunction(BoundaryFace::inner_x1,InnerBoundary);
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
    method = pin->GetReal("problem","mtd");
    rbc  = pin->GetReal("mesh","x1min");
   
    if ( method == 1 ) {
	if ( threshold <=1.0 ) {
		printf("Recommended threshold value is above 1.0 \n");
		printf("Your choice was = %6.4lf \n", threshold);
		printf("We automatically corrected it to 1.75 \n");
		Real threshold = 1.75;
	}
    }
    else if( method == 2 ) {
        if ( threshold >=1.0 ) {
                printf("Recommended threshold value is below 1.0 \n");
                printf("Your choice was = %6.4lf \n", threshold);
                printf("We automatically corrected it to 0.2 \n");
                Real threshold = 0.2;
        }
    }
    else{
	printf("mtd variable only takes the values of 1 or 2 \n");
	printf("We automatically corrected to mtd=1 and thr = 1.75 \n");
	Real method = 1 ;
	Real threshold = 1.75;
    }
   // EnrollUserBoundaryFunction(BoundaryFace::inner_x1,InnerBoundary);   
   // printf("Test4 \n");
  }
  
//LFLM inclusion ends
  return;
}


// LFLM inclusion start
//
int sign( double x ) {

        if ( x > 0.0 ) return 1;
        if ( x < 0.0 ) return -1;
        return 0;
}

//LFLM inclusion end
//

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //printf("Entering MeshBlock::ProblemGenerator \n");

  Real rout = pin->GetReal("problem", "radius");
  Real rin  = rout - pin->GetOrAddReal("problem", "ramp", 0.0);
  Real pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  Real da   = pin->GetOrAddReal("problem", "damb", 1.0);
  Real prat = pin->GetReal("problem", "prat");
  Real drat = pin->GetOrAddReal("problem", "drat", 1.0);
  Real b0, angle;
  //printf("Testing MeshBlock::ProblemGenerator 1\n");
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI/180.0)*pin->GetReal("problem", "angle");
  }
  //printf("Testing MeshBlock::ProblemGenerator 1.0\n");
  Real gamma = peos->GetGamma();
  //printf("Testing MeshBlock::ProblemGenerator 1.0.0\n");
  Real gm1 = gamma - 1.0;
  //printf("Testing MeshBlock::ProblemGenerator 1.1\n");
  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }
  //printf("Testing MeshBlock::ProblemGenerator 1.2\n");

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++){
        Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }

        Real den = da;
        if (rad < rout) {
          if (rad < rin) {
            den = drat*da;
          } else {   // add smooth ramp in density
            Real f = (rad-rin) / (rout-rin);
            Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
            den = std::exp(log_den);
          }
        }

        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          Real pres = pa;
          if (rad < rout) {
            if (rad < rin) {
              pres = prat*pa;
            } else {  // add smooth ramp in pressure
              Real f = (rad-rin) / (rout-rin);
              Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
              pres = std::exp(log_pres);
            }
          }
          phydro->u(IEN,k,j,i) = pres/gm1;
          if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
            phydro->u(IEN,k,j,i) += den;
        }
      }
    }
  }
  //printf("Testing MeshBlock::ProblemGenerator 2\n");
  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real phi = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) =
                b0 * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0 * std::abs(std::sin(theta))
                                   * (std::cos(angle) * std::cos(phi)
                                      + std::sin(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            pfield->b.x2f(k,j,i) = b0 * std::sin(angle);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real phi = pcoord->x2v(j);
            pfield->b.x2f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x2f(k,j,i) = b0 * std::cos(theta)
                                   * (std::cos(angle) * std::cos(phi)
                                      + std::sin(angle) * std::sin(phi));
            if (std::sin(theta) < 0.0)
              pfield->b.x2f(k,j,i) *= -1.0;
          }
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0
              || std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            pfield->b.x3f(k,j,i) = 0.0;
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real phi = pcoord->x3v(k);
            pfield->b.x3f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }
      }
    }
  }
  printf("Exiting MeshBlock::ProblemGenerator \n");
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  printf("Entering Mesh::UserWorkAfterLoop \n");

  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  // analysis - check shape of the spherical blast wave
  //int is = pblock->is, ie = pblock->ie;
  //int js = pblock->js, je = pblock->je;
  //int ks = pblock->ks, ke = pblock->ke;
  MeshBlock *pmb = my_blocks(0);
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  AthenaArray<Real> pr;
  //pr.InitWithShallowSlice(pblock->phydro->w, 4, IPR, 1);
  pr.InitWithShallowSlice(pmb->phydro->w, 4, IPR, 1);


  // get coordinate location of the center, convert to Cartesian
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // find indices of the center
  int ic, jc, kc;
  for (ic=is; ic<=ie; ic++)
    //if (pblock->pcoord->x1f(ic) > x1_0) break;
    if (pmb->pcoord->x1f(ic) > x1_0) break;
  ic--;
  //for (jc=pblock->js; jc<=pblock->je; jc++)
    //if (pblock->pcoord->x2f(jc) > x2_0) break;
  for (jc=pmb->js; jc<=pmb->je;jc++)
    if (pmb->pcoord->x2f(jc) > x2_0) break;
  jc--;
  //for (kc=pblock->ks; kc<=pblock->ke; kc++)
    //if (pblock->pcoord->x3f(kc) > x3_0) break;
  for (kc=pmb->ks; kc<pmb->ke; kc++)
    if (pmb->pcoord->x3f(kc) > x3_0) break;
  kc--;

  // search pressure maximum in each direction
  Real rmax = 0.0, rmin = 100.0, rave = 0.0;
  int nr = 0;
  for (int o=0; o<=6; o++) {
    int ios = 0, jos = 0, kos = 0;
    if (o == 1) ios=-10;
    else if (o == 2) ios =  10;
    else if (o == 3) jos = -10;
    else if (o == 4) jos =  10;
    else if (o == 5) kos = -10;
    else if (o == 6) kos =  10;
    for (int d=0; d<6; d++) {
      Real pmax = 0.0;
      int imax(0), jmax(0), kmax(0);
      if (d == 0) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i>=is; i--) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 1) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i<=ie; i++) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 2) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j>=js; j--) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 3) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j<=je; j++) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 4) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k>=ks; k--) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      } else { // if (d == 5) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k<=ke; k++) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      }

      Real xm, ym, zm;
      //Real x1m = pblock->pcoord->x1v(imax);
      //Real x2m = pblock->pcoord->x2v(jmax);
      //Real x3m = pblock->pcoord->x3v(kmax);

      Real x1m = pmb->pcoord->x1v(imax);
      Real x2m = pmb->pcoord->x2v(jmax);
      Real x3m = pmb->pcoord->x3v(kmax);

      if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
        xm = x1m*std::cos(x2m);
        ym = x1m*std::sin(x2m);
        zm = x3m;
      } else {  // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        xm = x1m*std::sin(x2m)*std::cos(x3m);
        ym = x1m*std::sin(x2m)*std::sin(x3m);
        zm = x1m*std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm-x0)+SQR(ym-y0)+SQR(zm-z0));
      if (rad > rmax) rmax = rad;
      if (rad < rmin) rmin = rad;
      rave += rad;
      nr++;
    }
  }
  rave /= static_cast<Real>(nr);

  // use physical grid spacing at center of blast
  Real dr_max;
  //Real  x1c = pblock->pcoord->x1v(ic);
  //Real dx1c = pblock->pcoord->dx1f(ic);
  //Real  x2c = pblock->pcoord->x2v(jc);
  //Real dx2c = pblock->pcoord->dx2f(jc);
  //Real dx3c = pblock->pcoord->dx3f(kc);

  Real  x1c = pmb->pcoord->x1v(ic);
  Real dx1c = pmb->pcoord->dx1f(ic);
  Real  x2c = pmb->pcoord->x2v(jc);
  Real dx2c = pmb->pcoord->dx2f(jc);
  Real dx3c = pmb->pcoord->dx3f(kc);

  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), dx3c);
  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), x1c*std::sin(x2c)*dx3c);
  }
  Real deform=(rmax-rmin)/dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
    }
    std::fprintf(pfile,"# Offset blast wave test in %s coordinates:\n",COORDINATE_SYSTEM);
    std::fprintf(pfile,"# Rmax       Rmin       Rave        Deformation\n");
    std::fprintf(pfile,"%e  %e  %e  %e \n",rmax,rmin,rave,deform);
    std::fclose(pfile);
  }
  //printf("Exiting Mesh::UserWorkAfterLoop \n");
  return;
}



// LFLM implementation of the reader for the outflow data start
// this should read the file
// import the data to an 3D AthenaArray, as in (var, theta, phi, t)
//

void OutflowReader(const char *filename, InterpTable3D *table) {

//first set some constants to translate the physical unities to code unities
//Some conversion factors to go between Lorene data and Athena data. Shamelessly stolen from
//the Einstein Toolkit's Mag_NS.cc.
//
	//printf("Test Reading Outflow 0 \n");
	printf(" Entering Outflow Reader\n");	

	Real const c_light = 299792458.0; // Speed of light [m/s]
	Real const mu0 = 4.0 * M_PI * 1.0e-7; // Vacuum permeability [N/A^2]
	Real const eps0 = 1.0 / (mu0 * pow(c_light, 2));
        Real const G_grav = 6.67428e-11; // Gravitational constant [m^3/kg/s^2]
        Real const M_sun = 1.98892e+30; // Solar mass [kg]
        Real const athenaM = M_sun;
        Real const athenaL = athenaM * G_grav / (c_light * c_light);
        Real const athenaT = athenaL / c_light;
        Real const athenaB = 1.0 / athenaL / sqrt(eps0 * G_grav / (c_light * c_light));
        Real const coord_unit = athenaL/1.0e3; // Convert to km for Lorene.
        Real const rho_unit = athenaM/(athenaL*athenaL*athenaL); // kg/m^3.
        Real const ener_unit = 1.0; // c^2
        Real const vel_unit = athenaL / athenaT / c_light; // c
        Real const B_unit = athenaB / 1.0e+9; // 10^9 T
        Real const mtodimensionless = c_light*c_light/(G_grav*M_sun);
        Real const kgbym3todimensionledd =pow(G_grav,3)*pow(M_sun,2)/pow(c_light,6);

	//printf("Test Reading outflow 0\n");
	
	std::ifstream fin;
	fin.open(filename);
	const int MAX_CHARS_PER_LINE  = 512;
	const int MAX_TOKENS_PER_LINE = 16;
	const char* const DELIMITER = "     ";
	int nline = 0;

        //printf("Test Reading outflow 0.0\n");

	
	if (!fin.good()){
		printf("Outflow file not found \n");
		return;
	}

        //printf("Test Reading outflow 0.0.0\n");

        //printf("BCMtheta = %6.4lf \n", BCMtheta);
        //printf("BCMphi   = %6.4lf \n", BCMphi)  ;
        //printf("BCMt     = %6.4lf \n", BCMt)    ;

        int BCNtheta = BCMtheta;
	int BCNphi   = BCMphi  ;
        int BCNt     = BCMt    ;
	
	//printf("BCNtheta = %d \n", BCNtheta);
	//printf("BCNphi   = %d \n", BCNphi)  ;
	//printf("BCNt     = %d \n", BCNt)    ;
	
	//printf("Type of BCNtheta is: \n");
	//std::cout << typeid(BCNtheta).name() << '\n';
        //printf("Type of BCNphi is: \n");
        //std::cout << typeid(BCNphi).name() << '\n';
	//printf("Type of BCNt is: \n");
	//std::cout << typeid(BCNt).name() << '\n';


	//table->SetSize(13,BCNtheta, BCNphi, BCNt);
	table->SetSize(5,BCNtheta+2*BCNghts,BCNphi+2*BCNghts,BCNt);
	
        //printf("Test Reading outflow 0.1\n");
	
	table -> SetY3lim(BCThetamin - BCNghts*BCDtheta ,BCThetamax + BCNghts*BCDtheta);
        table -> SetY2lim(BCPhimin   - BCNghts*BCDphi   ,BCPhimax   + BCNghts*BCDphi  );
        table -> SetY1lim(0                             ,BCTmax     - BCTmin          );
	//table -> data;	

	printf("thetamin = %6.40lf \n",BCThetamin - BCNghts*BCDtheta );
	printf("thetamax = %6.40lf \n",BCThetamax + BCNghts*BCDtheta );
	
	printf("phimin = %6.40lf \n", BCPhimin    - BCNghts*BCDphi);
        printf("phimax = %6.40lf \n", BCPhimax    + BCNghts*BCDphi);


	float Token[MAX_TOKENS_PER_LINE] ;
	
	//printf("Test Reading Outflow 1\n");
	while (!fin.eof()){
        	//printf("entering read loop\n");
        	// read an entire line into memory
        	char buf[MAX_CHARS_PER_LINE];
        	fin.getline(buf, MAX_CHARS_PER_LINE);
        	// parse the line into blank-delimited tokens
        	int n = 0; // a for-loop index
        	// array to store memory addresses of the tokens in buf
        	const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
        	// parse the line
        	token[0] = strtok(buf, DELIMITER); // first token
        	if (token[0]){ // zero if line is blank
        		for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
                        	//printf("entering token loop\n");
                                token[n] = strtok(0, DELIMITER); // subsequent tokens
                                //printf("going to transform the %d-th token\n",n);
                                Token[n] = std::stod(token[n]);
                                //Token[n] = (float)(token[n]);
                                //printf("Token= %6.40lf\n",Token[n]);
                                //if (!token[n]) break; // no more tokens
                        }
                }

		//printf("Test Reading Outflow 2 \n");

        	float BCtime      = Token[1] - BCTmin; //01-02 // time is already in code units it seems
						         // but I am shifting the time in order to 
						         // start the simulation of the ejecta from t=0
        	float BCx         = Token[2]; //02-03
        	float BCy         = Token[3]; //03-04
        	float BCz         = Token[4]; //04-05
        	float BCFluxd     = Token[5]; //05-06
        	float BCWlorentz  = Token[6]; //06-07
        	float BCeninf     = Token[7]; //07-08
        	float BCsufelem   = Token[8]; //08-09
       		float BCalp       = Token[9]; //09-10
        	float BCrho       = Token[10]; //10-11
        	float BCvx        = Token[11];// + velfloor; //11-12
        	float BCvy        = Token[12];// + velfloor; //12-13
        	float BCvz        = Token[13];// + velfloor; //13-14
        	float BCye        = Token[14]; //14-15
        	float BCentropy   = Token[15]; //15-16
	        float BCtemp      = Token[16]; //16-17
		float BCpress     = Token[17]; // pressure interpolated
        	
		//printf("BCx = %6.4lf \n", BCx);
		//printf("BCy = %6.4lf \n", BCy);
		//printf("BCz = %6.4lf \n", BCz);
		//printf("Test Reading Outflow 3 \n");

        	float BCradius = std::sqrt( pow(BCx ,2) + std::pow(BCy , 2) + std::pow(BCz , 2))                      ;
        	float BCtheta  = std::acos( BCz / std::sqrt(std::pow(BCx ,2) + std::pow(BCy , 2) + std::pow(BCz , 2)));
	        float BCphi    = sign(BCy) * std::acos(BCx /std::sqrt( std::pow(BCx ,2) + std::pow(BCy , 2))) + M_PI  ;

		//printf("In BC the coordinates are\n");
		//printf("    r = %6.4lf \n", BCradius);
		//printf("theta = %6.4lf \n", BCtheta);
		//printf("  phi = %6.4lf \n", BCphi);

	
                //float BCvp =-BCvx*std::sin(BCphi) + BCvy*std::cos(BCphi);
                //float BCvt = BCvx*std::cos(BCtheta)*std::cos(BCphi) + BCvy*std::cos(BCtheta)*std::sin(BCphi) - BCvz*std::sin(BCtheta);
                //float BCvr = BCvx*std::sin(BCtheta)*std::cos(BCphi) + BCvy*std::sin(BCtheta)*std::sin(BCphi) + BCvz*std::cos(BCtheta);
		
		float BCvr = ( BCx * BCvx + BCy * BCvy  + BCz * BCvz  ) / std::sqrt( BCx*BCx + BCy*BCy + BCz*BCz ) ;
		float BCvt = ( BCx * BCz * BCvx + BCy * BCz * BCvy - ( BCx*BCx + BCy*BCy ) * BCvz) / ( std::sqrt( BCx*BCx +BCy*BCy ) * ( BCx*BCx + BCy*BCy + BCz*BCz ) ) ;
		float BCvp = (- BCy * BCvx + BCx * BCvy )/(BCx*BCx + BCy*BCy);	       

		//printf("BCradius=%6.16lf \n", BCradius);
                //printf("BCtheta =%6.16lf \n", BCtheta );
                //printf("BCphi   =%6.16lf \n", BCphi   );

	        float BCmtheta = (BCtheta - BCThetamin)/BCDtheta;
		float BCmphi   = (BCphi   - BCPhimin)  /BCDphi  ;
		float BCmt     = (BCtime  - BCTmin)    /BCDt    ;

	//	printf("BCmtheta = %6.4lf\n", BCmtheta);
	//	printf("BCmphi   = %6.4lf\n", BCmphi);
        //      printf("BCmt     = %6.4lf\n", BCmt);

		//int BCntheta = std::floor(BCmtheta);// + 1;
		//int BCnphi   = std::floor(BCmphi)  ;// + 1;
		//int BCnt     = std::floor(BCmt)    ;// + 1;

		int BCntheta = BCmtheta +BCNghts;
		int BCnphi   = BCmphi   +BCNghts;
		int BCnt     = BCmt             ;

		if (BCntheta >BCNtheta+2*BCNghts or BCntheta <0) {
			printf(" Theta index is outside of range  \n");
			break;
			return;
		}

                if (BCnphi >BCNphi+2*BCNghts or BCnphi <0) {
                        printf(" Phi index is outside of range  \n");
                        break;
			return;
		}
		if (BCnt > BCNt or BCnt <0) {
			printf("Time index is outside of range \n");
			break;
			return;
		}

		//printf("*********************************\n");

                //printf("BCtheta = %6.10lf \n",BCtheta);
                //printf("  BCphi = %6.10lf \n",BCphi);
                //printf("    BCt = %6.10lf \n",BCtime);

		//printf("BCntheta = %d \n",BCntheta);
		//printf("  BCnphi = %d \n",BCnphi);
		//printf("    BCnt = %d \n",BCnt);

                if (BCtheta < 0) {
                        printf("wrong choice of parametrization for the interpolation in theta \n");
                        break;
                }

                if (BCphi < 0) {
                        printf("wrong choice of parametrization for the interpolation in phi \n");
                        break;
                }
                //printf("*********************************\n");
                //printf("BCntheta = %d \n",BCntheta);
                //printf("BCnphi = %d \n",BCnphi);
                //printf("BCnt = %d \n",BCnt);


	        //printf("Type of BCntheta is: \n");
       		//std::cout << typeid(BCntheta).name() << '\n';
	        //printf("Type of BCnphi is: \n");
	        //std::cout << typeid(BCnphi).name() << '\n';
	        //printf("Type of BCnt is: \n");
	        //std::cout << typeid(BCnt).name() << '\n';

		
		table->data(0 , BCntheta, BCnphi, BCnt ) = BCrho; //density 
		//printf("density was set \n"); 
		//printf("rho = %6.40lf \n", BCrho);
		//table->data(1 , BCntheta, BCnphi, BCnt ) = Kval * std::pow(BCrho,Gammaval); //pressure - here I am assuming polytropic EOS
                table->data(1,BCntheta, BCnphi, BCnt ) = BCpress; // pressure - here I am using pressure as read from the table with pressure been previous interpolated by the eos
		//printf("pressure was set \n");
                //printf("press = %6.40lf \n", Kval * std::pow(BCrho,Gammaval));
		table->data(2 , BCntheta, BCnphi, BCnt ) = BCvr;      // v_r
                //printf("vr was set \n");
                //printf("vr = %6.40lf \n", BCvr);
                table->data(3 , BCntheta, BCnphi, BCnt ) = BCvt;      // v_theta
                //printf("vt was set \n");
                //printf("vtheta = %6.40lf \n", BCvt);
		table->data(4 , BCntheta, BCnphi, BCnt ) = BCvp;      // v_phi
                //printf("vp was set \n");
		//printf("vphi = %6.40lf \n", BCvp);

                //table->data(5 , BCntheta, BCnphi, BCnt ) = BCtemp;    //temperature
                //printf("tmp was set \n");
		//table->data(6 , BCntheta, BCnphi, BCnt ) = BCentropy; // entropy 
                //printf("ent was set \n");                
		//table->data(7 , BCntheta, BCnphi, BCnt ) = BCye ;     // Ye
                //printf("ye was set \n");                
		//table->data(8 , BCntheta, BCnphi, BCnt ) = BCWlorentz;//W_lorentz
                //printf("w was set \n");		
		//table->data(9 , BCntheta, BCnphi, BCnt ) = BCeninf;   //eninf
                //printf("eninf was set \n");
		//table->data(10, BCntheta, BCnphi, BCnt ) = BCsufelem; //surface_element
                //printf("surf_elem was set \n");	        
		//table->data(11, BCntheta, BCnphi, BCnt ) = BCFluxd;   //fluxdens
                //printf("fluxdens was set \n");
		//table->data(12, BCntheta, BCnphi, BCnt ) = BCalp;     //alp
                //printf("alp was set \n");
              //  printf(" BCrho =%6.16lf \n",BCrho);
              //  printf(" BCpress = %6.16lf \n",std::pow(BCrho,Gammaval));
              //  printf(" BCvr = %6.16lf \n",BCvr);

               // printf("Done reading line %d \n", nline);
		//printf("Test Reading Outflow 4 \n");
		nline = nline+1;

	}
	//table -> SetX3lim(BCThetamin ,BCThetamax ) ;
	//table -> SetX2lim(BCPhimin   ,BCPhimax   ) ;
	//table -> SetX1lim(0     ,BCTmax - BCTmin ) ;

	//table -> data;

	printf("Filling ghost zones \n ");
	for (int kk = 0 ; kk <= BCNt ; kk++) {
		for (int ii =0 ; ii<= BCNghts-1 ; ii++){
			for (int jj =0; jj <= BCNghts-1 ; jj++){

				//printf("======================= \n");
				//printf("kk = %d \n",kk);
				//printf("jj = %d \n",jj);
				//printf("ii = %d \n",ii);				
				//printf("======================= \n");

				table->data(0 , ii, jj, kk ) = table->data(0 , BCNtheta+BCNghts-ii, BCNphi+BCNghts-jj, kk ); //density 
                		table->data(1 , ii, jj, kk ) = table->data(1 , BCNtheta+BCNghts-ii, BCNphi+BCNghts-jj, kk ); //pressure
				table->data(2 , ii, jj, kk ) = table->data(2 , BCNtheta+BCNghts-ii, BCNphi+BCNghts-jj, kk ); //vr
				table->data(3 , ii, jj, kk ) = table->data(3 , BCNtheta+BCNghts-ii, BCNphi+BCNghts-jj, kk ); //vt
				table->data(4 , ii, jj, kk ) = table->data(4 , BCNtheta+BCNghts-ii, BCNphi+BCNghts-jj, kk ); //vp
				//printf("AA \n");

                        	table->data(0 , ii, BCNphi+2*BCNghts-jj-1, kk ) = table->data(0 , BCNtheta+BCNghts-ii, BCNghts-jj, kk ); //density 
                        	table->data(1 , ii, BCNphi+2*BCNghts-jj-1, kk ) = table->data(1 , BCNtheta+BCNghts-ii, BCNghts+jj, kk ); //pressure
                        	table->data(2 , ii, BCNphi+2*BCNghts-jj-1, kk ) = table->data(2 , BCNtheta+BCNghts-ii, BCNghts+jj, kk ); //vr
                        	table->data(3 , ii, BCNphi+2*BCNghts-jj-1, kk ) = table->data(3 , BCNtheta+BCNghts-ii, BCNghts+jj, kk ); //vt
                        	table->data(4 , ii, BCNphi+2*BCNghts-jj-1, kk ) = table->data(4 , BCNtheta+BCNghts-ii, BCNghts+jj, kk ); //vp
				//printf("BB \n");

                        	table->data(0 , BCNtheta+2*BCNghts-ii-1, jj, kk ) = table->data(0 , BCNghts+ii, BCNphi+BCNghts-jj, kk ); //density 
                        	table->data(1 , BCNtheta+2*BCNghts-ii-1, jj, kk ) = table->data(1 , BCNghts+ii, BCNphi+BCNghts-jj, kk ); //pressure
                        	table->data(2 , BCNtheta+2*BCNghts-ii-1, jj, kk ) = table->data(2 , BCNghts+ii, BCNphi+BCNghts-jj, kk ); //vr
                        	table->data(3 , BCNtheta+2*BCNghts-ii-1, jj, kk ) = table->data(3 , BCNghts+ii, BCNphi+BCNghts-jj, kk ); //vt
                        	table->data(4 , BCNtheta+2*BCNghts-ii-1, jj, kk ) = table->data(4 , BCNghts+ii, BCNphi+BCNghts-jj, kk ); //vp
				//printf("CC \n");

                        	table->data(0 , BCNtheta+2*BCNghts-ii-1, BCNphi+2*BCNghts-jj-1, kk ) = table->data(0 , BCNghts+ii, BCNghts+jj, kk ); //density 
                        	table->data(1 , BCNtheta+2*BCNghts-ii-1, BCNphi+2*BCNghts-jj-1, kk ) = table->data(1 , BCNghts+ii, BCNghts+jj, kk ); //pressure
                        	table->data(2 , BCNtheta+2*BCNghts-ii-1, BCNphi+2*BCNghts-jj-1, kk ) = table->data(2 , BCNghts+ii, BCNghts+jj, kk ); //vr
                        	table->data(3 , BCNtheta+2*BCNghts-ii-1, BCNphi+2*BCNghts-jj-1, kk ) = table->data(3 , BCNghts+ii, BCNghts+jj, kk ); //vt
                        	table->data(4 , BCNtheta+2*BCNghts-ii-1, BCNphi+2*BCNghts-jj-1, kk ) = table->data(4 , BCNghts+ii, BCNghts+jj, kk ); //vp
				//printf("DD \n");
			}
		}
	}
	printf("Done filling ghost zones \n ");
	printf("Done with reading the data \n ");
	return ;
}


//LFLM implementation of the reader for the outflow data end

// LFLM implementation of the interpolation starts
//

Real DenProfile(const Real th, const Real ph, const Real tm) {
  Real den;
  den = table->interpolate(0, th, ph, tm);
  return den ;
}

Real PressProfile(const Real th, const Real ph, const Real tm) {
  Real press;
  press = table->interpolate(1, th, ph, tm);
  return press ;
}

Real VRProfile(const Real th, const Real ph, const Real tm) {
  Real velR;
  velR = table->interpolate(2, th, ph, tm);
  return  velR ;
}

Real VTProfile(const Real th, const Real ph, const Real tm) {
  Real velT;
  velT = table->interpolate(3, th, ph, tm);
  return velT ;
}


Real VPProfile(const Real th, const Real ph, const Real tm) {
  Real velP;
  velP = table->interpolate(4, th, ph, tm);
  return velP ;
}




// LFLM inplementation of the intepolation ends
//


// LFLM implementation of the inner boundary condition start

void InnerBoundary(MeshBlock *pmb, Coordinates *pcoord,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku,
                  int ngh) {
        //AthenaArray<Real> &w = pmb->phydro->w;
        //AthenaArray<Real> prim = pmb->phydro->w; 
        printf("Entering BC condition \n");
	//for (int n=0; n<=3; ++n) {
	//	printf("Mesh block =%d \n",n);
		     	
	float BCThetamaxval = BCThetamax + BCNghts*BCDtheta;
        float BCThetaminval = BCThetamin - BCNghts*BCDtheta;
        float BCPhimaxval   = BCPhimax   + BCNghts*BCDphi;
        float BCPhiminval   = BCPhimin   - BCNghts*BCDphi;
        float DeltaPhi      = BCNghts*BCDphi;
        float DeltaTheta    = BCNghts*BCDtheta;
        float BCTmaxval     = BCTmax + BCtin;
        
	//printf("BCThetaminval = %6.40lf \n", BCThetaminval );
	//printf("BCThetamaxval = %6.40lf \n", BCThetamaxval );

	//printf("BCPhiminval = %6.40lf\n", BCPhiminval );
	//printf("BCPhimaxval = %6.40lf\n", BCPhimaxval );


        //printf("Deltaphi   = %6.40lf \n", DeltaPhi);
        //printf("Deltatheta = %6.40lf \n", DeltaTheta);

        //printf("BCNghts = %6.40lf \n", BCNghts);

        //printf("BCDphi   = %6.40lf \n", BCDphi);
        //printf("BCDtheta = %6.40lf \n", BCDtheta);


	for (int k=kl; k<=ku; ++k) {
		//printf("Inside the BC loop \n");
               	for (int j=jl; j<=ju; ++j){
                       	for (int i=1; i<=ngh; ++i){

		                Real rval   = pcoord->x1v(il-i) ;
                               	Real thval  = pcoord->x2v(j) ;
                               	Real phival = pcoord->x3v(k) ;
				
				float timeshift     = time   + BCtin;

				
                                //printf("     r = %6.4lf \n", rval);
                                //printf(" theta = %6.4lf \n", thval);
                                //printf("phival = %6.4lf \n", phival);
                                //printf("  time = %6.4lf \n", timeshift);

				if (thval > BCThetamaxval or thval < BCThetaminval){
					printf(" Theta outside of interpolation range  \n");
				        printf("theta = %6.40lf \n", thval);
					printf("thetamin = %6.40lf \n", BCThetaminval );
        				printf("thetamax = %6.40lf \n", BCThetamaxval );
					break;
				}
                                if (phival > BCPhimaxval or phival < BCPhiminval){
                                        printf(" Phi outside of interpolation range  \n");
					printf("phi = %6.40lf \n", phival);
        				printf("phimin  = %6.40lf \n", BCPhiminval);
        				printf("phimax  = %6.40lf \n", BCPhimaxval);
                                        break;
                                }
                                if (timeshift > BCTmaxval or timeshift < 0 ){
                                        printf(" time outside of interpolation range  \n");
					printf("timeshift = %6.40lf \n", timeshift);
                                        printf("tmin = %6.40lf \n", 0);
        				printf("tmax = %6.40lf \n", BCTmaxval);

					break;
                                }
                                //printf("==================================================================== \n");
				if (std::isfinite(DenProfile(thval, phival, timeshift)) == 1 ) {
					prim(IDN, k ,j ,il-i) = DenProfile(thval, phival, timeshift);
					//printf("           density = %6.40lf \n", prim(IDN, k, j, il-i));
					//printf("from interpolation = %6.40lf \n", DenProfile(thval, phival, time));
					//printf("==================================================================== \n");
				}
				else {
					printf("Error in setting density! \n");
					break;
				}
                                if (std::isfinite(PressProfile(thval, phival, timeshift)) == 1 ) {

	                                prim(IPR, k, j, il-i) = PressProfile(thval, phival, timeshift);
					//printf("          pressure = %6.40lf \n", prim(IPR, k, j, il-i));
                	                //printf("from interpolation = %6.40lf \n", PressProfile(thval, phival, time));
                        	        //printf("==================================================================== \n");
                               	}
				else {
                                        printf("Error in setting pressure! \n");
                                        break;
				}

                                if (std::isfinite(VRProfile(thval, phival, time + BCtin)) == 1 ) {


					prim(IVX, k, j, il-i) = VRProfile(thval, phival, timeshift);
					//printf("                vr = %6.40lf \n", prim(IVX, k, j, il-i));
                                	//printf("from interpolation = %6.40lf \n", VRProfile(thval, phival, time));
                                	//printf("==================================================================== \n");
                               		if (std::abs(VRProfile(thval, phival, timeshift)) >= 1.0 ) {
						printf("!!!!!!!!!!!!!!!!!!  Radial  Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                		printf("     r = %6.4lf \n", rval);
                                                printf(" theta = %6.4lf \n", thval);
                                                printf("phival = %6.4lf \n", phival);
                                                printf("  time = %6.4lf \n", timeshift);
                                		printf("    Vr = %6.40lf \n", VRProfile(thval, phival, timeshift));			
						printf(" !!!!!!!!!!!!!!!!!!  Radial  Velocity is too big !!!!!!!!!!!!!!!!! \n");
						break;
					}
					//else {
					//	printf("All good \n");
					//}
					
				}
				else {
                                        printf("Error in setting radial velocity! \n");
                                        break;
				}


				if (std::isfinite(VTProfile(thval, phival, timeshift)) == 1 ) {

                                	prim(IVY, k, j, il-i) = VTProfile(thval, phival, timeshift);
					//printf("            vtheta = %6.40lf \n", prim(IVY, k, j, il-i));
                                	//printf("from interpolation = %6.40lf \n", VTProfile(thval, phival, time));
                                	//printf("==================================================================== \n");
                                        if (std::abs(VRProfile(thval, phival, timeshift)) >= 1.0 ) {

                                                printf(" !!!!!!!!!!!!!!!!!! Theta Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                                printf("     r = %6.4lf \n", rval);
                                                printf(" theta = %6.4lf \n", thval);
                                                printf("phival = %6.4lf \n", phival);
                                                printf("  time = %6.4lf \n", timeshift);
						printf("    Vt = %6.40lf \n", VTProfile(thval, phival, timeshift));
                                                printf(" !!!!!!!!!!!!!!!!!! Theta Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                break;
                                        }
                                        //else {
                                        //        printf("All good \n");
                                        //}
                               	

				}
                                else {
                                        printf("Error in setting theta velocity! \n");
                                        break;
                                }

                                
                                if (std::isfinite(VPProfile(thval, phival, timeshift)) == 1 ) {
					prim(IVZ, k , j, il-i) = VPProfile(thval, phival, timeshift);
					//printf("              vphi = %6.40lf \n", prim(IVZ, k, j, il-i));
                                        //printf("from interpolation = %6.40lf \n", VPProfile(thval, phival, time));
                                	//printf("==================================================================== \n");
                                        if (std::abs(VRProfile(thval, phival, timeshift)) >= 1.0 ) {

                                                printf(" !!!!!!!!!!!!!!!!!! Phi Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                                printf("     r = %6.4lf \n", rval);
                                                printf(" theta = %6.4lf \n", thval);
                                                printf("phival = %6.4lf \n", phival);
                                                printf("    Vp = %6.4lf \n", timeshift);
						printf("from interpolation = %6.40lf \n", VPProfile(thval, phival, timeshift)); 
                                                printf(" !!!!!!!!!!!!!!!!!! Phi Velocity is too big !!!!!!!!!!!!!!!!! \n");


                                                break;
                                        }
                                       // else {
                                       //         printf("All good \n");
                                       // }


				}
				else {
                                        printf("Error in setting phi velocity! \n");
                                        break;
                                }

                                //printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
                                //printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
                                //printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
                                //printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");


				//printf("Testing bc 2 \n");
				//printf("####################################################################### \n");
                       	}
               	}
        }
	//}
	//printf("Exiting BC condition \n");
        return;
}



// LFLM implementation of the inner boundary condition end


// LFLM condition for mesh refinement
// STARTS
//

int RefinementCondition(MeshBlock *pmb) {
	//printf("Entering RefinementCondition \n");
	if ( method == 1 ) {

	      AthenaArray<Real> &w = pmb->phydro->w;
	      int test   = -1;
	      int temp   = -1;
	      float rhosum = 0.0;
	      float psum   = 0.0;
	      float Npoints= 0.0;
//	      float threshold = 1.75;
	      for (int k=pmb->ks; k<=pmb->ke; k++) {
        	        for (int j=pmb->js; j<=pmb->je; j++) {
                	       for (int i=pmb->is; i<=pmb->ie; i++) {
                        	      rhosum = rhosum + w(IDN,k,j,i);
                             	      psum   = psum   + w(IPR,k,j,i);
                             	      Npoints= Npoints + 1.0;

                       		 }
                	}
        	}

	      float rhoavg = rhosum/Npoints;
	      float pavg   = psum/Npoints;

     	     // printf("We are inside the refinement condition \n");
      	      for (int k=pmb->ks; k<=pmb->ke; k++) {
              		for (int j=pmb->js; j<=pmb->je; j++) {
                      		for (int i=pmb->is; i<=pmb->ie; i++) {
					                                                         
                                	int drx1adv = sign( w(IDN,k,j,i+1) - w(IDN,k,j,i) );
                                	int drx2adv = sign( w(IDN,k,j+1,i) - w(IDN,k,j,i) );
                                	int drx3adv = sign( w(IDN,k+1,j,i) - w(IDN,k,j,i) );

                                	int dpx1adv = sign( w(IPR,k,j,i+1) - w(IPR,k,j,i) );
                                	int dpx2adv = sign( w(IPR,k,j+1,i) - w(IPR,k,j,i) );
                                	int dpx3adv = sign( w(IPR,k+1,j,i) - w(IPR,k,j,i) );


                                	int drx1ret = sign( w(IDN,k,j,i) - w(IDN,k,j,i-1) );
                                	int drx2ret = sign( w(IDN,k,j,i) - w(IDN,k,j-1,i) );
                                	int drx3ret = sign( w(IDN,k,j,i) - w(IDN,k-1,j,i) );

                                	int dpx1ret = sign( w(IPR,k,j,i) - w(IPR,k,j,i-1) );
                                	int dpx2ret = sign( w(IPR,k,j,i) - w(IPR,k,j-1,i) );
                                	int dpx3ret = sign( w(IPR,k,j,i) - w(IPR,k-1,j,i) );
					Real rval = pmb-> pcoord -> x1v(i) ;
					//if (rval <= 1.2 * rbc) {
					//	temp = -1;
					//}
                              		//else 
                              		if ( w(IDN,k,j,i) >= threshold * rhoavg || w(IPR,k,j,i) >= threshold * pavg ) {
                                        	 if (drx1ret != drx1adv || drx2ret != drx2adv || drx3ret != drx3adv || dpx1ret != dpx1adv || dpx2ret != dpx2adv || dpx3ret != dpx3adv) {
                                                	temp = 1;
                                             	}
                                       		 else {
                                               		temp = 0;
                                             	 }
                               		 }
                              		else {
                                        	 temp =-1 ;
                                  	 }
                              		test = std::max(test,temp);
                    		 }
              		}
      		}
     		return test;
      }
      if ( method ==2 ) {
      		AthenaArray<Real> &w = pmb->phydro->w;
        	Real maxeps=0.0;
        	for(int k=pmb->ks; k<=pmb->ke; k++) {
                	for(int j=pmb->js; j<=pmb->je; j++) {
                        	for(int i=pmb->is; i<=pmb->ie; i++) {
                                	 Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                                        	    +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i))
                                            	    +std::abs(w(IDN,k+1,j,i)-2.0*w(IDN,k,j,i)+w(IDN,k-1,j,i)))/w(IDN,k,j,i);
                                 	 Real epsp= (std::abs(w(IPR,k,j,i+1)-2.0*w(IPR,k,j,i)+w(IPR,k,j,i-1))
                                         	    +std::abs(w(IPR,k,j+1,i)-2.0*w(IPR,k,j,i)+w(IPR,k,j-1,i))
                                            	    +std::abs(w(IPR,k+1,j,i)-2.0*w(IPR,k,j,i)+w(IPR,k-1,j,i)))/w(IPR,k,j,i);
                                 	 Real eps = std::max(epsr, epsp);
                                 	 maxeps = std::max(maxeps, eps);
                        	}
                	}
         	}
         	if(maxeps > threshold) return 1;
         	if(maxeps < threshold/2) return -1;
         	return 0;

	}
	//printf("Exiting RefinementCondition \n");

}


// LFLM condition for mesh refinement
// END



// ORIGINAL REFINEMENT CONDITION STARTS

// refinement condition: check the maximum pressure gradient
//int RefinementCondition(MeshBlock *pmb) {
//  AthenaArray<Real> &w = pmb->phydro->w;
//  Real maxeps = 0.0;
//  if (pmb->pmy_mesh->f3) {
//    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
//      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
//        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
//          Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
//                               +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
//                               +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
//          maxeps = std::max(maxeps, eps);
//        }
//      }
//    }
//  } else if (pmb->pmy_mesh->f2) {
//    int k = pmb->ks;
//    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
//      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
//        Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
//                             + SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
//        maxeps = std::max(maxeps, eps);
//      }
//    }
//  } else {
//    return 0;
//  }

//  if (maxeps > threshold) return 1;
//  if (maxeps < 0.25*threshold) return -1;
//  return 0;
//}

// ORIGINAL REFINEMENT CONDITION ENDS
