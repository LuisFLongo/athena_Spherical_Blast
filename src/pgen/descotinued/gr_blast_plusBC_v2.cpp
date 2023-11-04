//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_blast.cpp
//! \brief Problem generator for GRMHD spherical blast wave in flat spacetime.


// COMMENT : this problem sits on top of gr_blast_plusBC.cpp
// 	     the only modification is the user defined BC 
// 	     while the coordinate system is maintained 
// 	     i will be reading GW170817 
//           The BC and data reading are copied from
//           blast_try5_BC_v9.cpp
//

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <cstring>    // strcmp()
#include <iostream>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>

// Athena++ headers
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../globals.hpp"
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../parameter_input.hpp"          // ParameterInput
#include "../utils/interp_table3D.hpp"


// Configuration checking
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif

// Declarations
namespace {
void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
                             Real *px, Real *py, Real *pz);
void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3);
Real DistanceBetweenPoints(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3);


Real threshold, method, da, pa;
float BCThetamin, BCPhimin, BCTmin, BCThetamax, BCPhimax, BCTmax, BCMtheta, BCMphi,BCMt , BCtin , BCDtheta, BCDphi, BCDt; //, BCMghts
Real Gammaval, Kval;
float Trhomax, Tpressmax, Tvrmax, Tvtmax, Tvpmax;
float Trhomin, Tpressmin, Tvrmin, Tvtmin, Tvpmin;
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
Real IntEnProfile(const Real th, const Real ph, const Real tm);
Real WProfile(const Real th, const Real ph, const Real tm);


int sign( double x ) {

        if ( x > 0.0 ) return 1;
        if ( x < 0.0 ) return -1;
        return 0;
}

} // namespace

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters

void Mesh::InitUserMeshData(ParameterInput *pin) {


  Real pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  Real da   = pin->GetOrAddReal("problem", "damb", 1.0);
  Gammaval  = pin->GetOrAddReal("hydro", "gamma", 1.0);

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

  BCDtheta   = (BCThetamax-BCThetamin)/(BCMtheta-1);
  BCDphi     = (BCPhimax-BCPhimin)/(BCMphi-1)      ;
  BCDt       = (BCTmax-BCTmin)/(BCMt-1)            ;

  float Nthetadown = ( BCThetamin -  (-2*M_PI/2) )/BCDtheta;
  float Nthetaup   = (  4*M_PI/2  - BCThetamax )/BCDtheta;

  float Nphidown   = ( BCPhimin - (-2*M_PI/2) )/BCDphi;
  float Nphiup     = ( 6*M_PI/2 -  BCPhimax )/BCDphi;
  BCNghts = std::ceil(std::max(std::max(Nthetadown, Nthetaup),std::max(Nphidown, Nphiup)));

 printf("Type of BCNghts is: \n");
  std::cout << typeid(BCNghts).name() << '\n';

  printf("Number of ghost set to = %d \n", BCNghts);
  float A1 = BCThetamax + BCNghts*BCDtheta;
  float A2 = BCThetamin - BCNghts*BCDtheta;
  printf("New theta max = %6.40lf \n", A1 );
  printf("New theta min = %6.40lf \n", A2 );

  float A3 = BCPhimax + BCNghts*BCDphi;
  float A4 = BCPhimin - BCNghts*BCDphi;
  printf("New phi max = %6.40lf \n", A3 );
  printf("New phi min = %6.40lf \n", A4 );

  Trhomin = 9.99908793 * std::pow(10,-15);
  Trhomax = 2.04474417 * std::pow(10,-10);

  Tpressmin = 8.254368620859124 * std::pow(10,-20);
  Tpressmax = 1.952452818229532 * std::pow(10,-14);

  Tvrmin = -0.0037080970266283662;
  Tvrmax =  0.5145687580274914;

  Tvtmin = -0.00025542644782544616;
  Tvtmax =  0.00025542644782544616;

  Tvpmin = -0.004348400171912313;
  Tvpmax =  0.002945209113444522;

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1,InnerBoundary);

  std::string fn;
  fn = pin->GetString("problem","outflow_table");

  table = new InterpTable3D();
  OutflowReader(fn.c_str(), table);


  return;
}



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
  AthenaArray<Real> b, g, gi;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);
  g.NewAthenaArray(NMETRIC, ncells1);
  gi.NewAthenaArray(NMETRIC, ncells1);

  // Initialize hydro variables
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i=il; i<=iu; ++i) {
        // Calculate distance to nearest blast center
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real min_separation = DistanceBetweenPoints(x1, x2, x3, 0.0, 0.0, 0.0);
        for (int x_index = -num_x; x_index <= num_x; ++x_index) {
          Real center_x = x_index * x_spacing;
          for (int y_index = -num_y; y_index <= num_y; ++y_index) {
            if (x_index == 0 && y_index == 0) {
              continue;
            }
            Real center_y = y_index * y_spacing;
            Real separation = DistanceBetweenPoints(x1, x2, x3, center_x, center_y, 0.0);
            min_separation = std::min(min_separation, separation);
          }
        }

        // Set pressure and density
        Real rho, pgas;
        if (min_separation < radius) {
          rho = rho_inner;
          pgas = pgas_inner;
        } else {
          rho = rho_outer;
          pgas = pgas_outer;
        }

        // Get Minkowski coordinates of point
        Real t, x, y, z;
        GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);

        // Set velocity
        Real ut = 1.0;
        Real ux = 0.0;
        Real uy = 0.0;
        Real uz = 0.0;
        Real u0, u1, u2, u3;
        TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = u1 - gi(I01,i)/gi(I00,i) * u0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = u2 - gi(I02,i)/gi(I00,i) * u0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = u3 - gi(I03,i)/gi(I00,i) * u0;

        // Calculate cell-centered magnetic fields given Minkowski values
        Real bcont = 0.0;
        Real bconx = bx;
        Real bcony = by;
        Real bconz = bz;
        Real bcon0, bcon1, bcon2, bcon3;
        TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                        &bcon3);
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Delete auxiliary array

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
        for (int i=il; i<=iu+1; ++i) {
          Real ut = 1.0;
          Real ux = 0.0;
          Real uy = 0.0;
          Real uz = 0.0;
          Real bcont = 0.0;
          Real bconx = bx;
          Real bcony = by;
          Real bconz = bz;
          Real u0, u1, u2, u3;
          Real bcon0, bcon1, bcon2, bcon3;
          if (j != ju+1 && k != ku+1) {
            Real x1 = pcoord->x1f(i);
            Real x2 = pcoord->x2v(j);
            Real x3 = pcoord->x3v(k);
            Real t, x, y, z;
            GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
            TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
            TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                            &bcon3);
            pfield->b.x1f(k,j,i) = bcon1 * u0 - bcon0 * u1;
          }
          if (i != iu+1 && k != ku+1) {
            Real x1 = pcoord->x1v(i);
            Real x2 = pcoord->x2f(j);
            Real x3 = pcoord->x3v(k);
            Real t, x, y, z;
            GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
            TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
            TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                            &bcon3);
            pfield->b.x2f(k,j,i) = bcon2 * u0 - bcon0 * u2;
          }
          if (i != iu+1 && j != ju+1) {
            Real x1 = pcoord->x1v(i);
            Real x2 = pcoord->x2v(j);
            Real x3 = pcoord->x3f(k);
            Real t, x, y, z;
            GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
            TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
            TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                            &bcon3);
            pfield->b.x3f(k,j,i) = bcon3 * u0 - bcon0 * u3;
          }
        }
      }
    }
  }
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

void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
                             Real *px, Real *py, Real *pz) {
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    *pt = x0;
    *px = x1;
    *py = x2;
    *pz = x3;
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

void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    *pa0 = at;
    *pa1 = ax;
    *pa2 = ay;
    *pa3 = az;
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

//***********************************************************

 void OutflowReader(const char *filename, InterpTable3D *table) {

        std::ifstream fin;
        fin.open(filename);
        const int MAX_CHARS_PER_LINE  = 512;
        const int MAX_TOKENS_PER_LINE = 19;
        const char* const DELIMITER = "     ";
        int nline = 0;

        if (!fin.good()){
                printf("Outflow file not found \n");
                return;
        }

        int BCNtheta = BCMtheta;
        int BCNphi   = BCMphi  ;
        int BCNt     = BCMt    ;

        table->SetSize(8,BCNtheta+2*BCNghts,BCNphi+2*BCNghts,BCNt);

        table -> SetY3lim(BCThetamin - BCNghts*BCDtheta ,BCThetamax + BCNghts*BCDtheta);
        table -> SetY2lim(BCPhimin   - BCNghts*BCDphi   ,BCPhimax   + BCNghts*BCDphi  );
        table -> SetY1lim(0                             ,BCTmax     - BCTmin          );

        float Token[MAX_TOKENS_PER_LINE] ;

        while (!fin.eof()){
                char buf[MAX_CHARS_PER_LINE];
                fin.getline(buf, MAX_CHARS_PER_LINE);
                int n = 0; // a for-loop index
                const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
                token[0] = strtok(buf, DELIMITER); // first token
                if (token[0]){ // zero if line is blank
                        for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
                                token[n] = strtok(0, DELIMITER); // subsequent tokens
                                Token[n] = std::stod(token[n]);
                        }
                }

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
                float BCinter     = Token[18];//internal energy interpolated

                float BCradius  = std::sqrt( pow(BCx ,2) + std::pow(BCy , 2) + std::pow(BCz , 2))                      ;
                float BCtheta1  = std::acos( BCz / std::sqrt(std::pow(BCx ,2) + std::pow(BCy , 2) + std::pow(BCz , 2)));
                float BCphi1    = sign(BCy) * std::acos(BCx /std::sqrt( std::pow(BCx ,2) + std::pow(BCy , 2))) + M_PI  ;

                float BCtheta  = std::max(std::min(BCtheta1,BCThetamax),BCThetamin);
                float BCphi    = std::max(std::min(BCphi1,BCPhimax),BCPhimin);

                float BCvr = ( BCx * BCvx + BCy * BCvy  + BCz * BCvz  ) / std::sqrt( BCx*BCx + BCy*BCy + BCz*BCz ) ;
                float BCvt = ( BCx * BCz * BCvx + BCy * BCz * BCvy - ( BCx*BCx + BCy*BCy ) * BCvz) / ( std::sqrt( BCx*BCx +BCy*BCy ) * ( BCx*BCx + BCy*BCy + BCz*BCz ) ) ;
                float BCvp = (- BCy * BCvx + BCx * BCvy )/(BCx*BCx + BCy*BCy);
                                                                                            
                int BCntheta = round( (BCtheta - BCThetamin) /BCDtheta + BCNghts);
                int BCnphi   = round( (BCphi   - BCPhimin  ) /BCDphi   + BCNghts);
                int BCnt     = round( (BCtime  - BCTmin    ) /BCDt    );                    
                if (BCntheta > BCNtheta+BCNghts+1 or BCntheta < BCNghts) {
                        printf(" Theta index is outside of range  \n");
                }

                if (BCnphi > BCNphi+BCNghts+1 or BCnphi < BCNghts) {
                        printf(" Phi index is outside of range  \n");
                }
                if (BCnt > BCNt+1 or BCnt < 0) {
                        printf("Time index is outside of range \n");
                }
                if (BCtheta < 0) {
                        printf("wrong choice of parametrization for the interpolation in theta \n");
                }
                if (BCphi < 0) {
                        printf("wrong choice of parametrization for the interpolation in phi \n");
                }
                table->data(0 , BCntheta, BCnphi, BCnt ) = BCrho; //density
                if ( BCrho < Trhomin or  BCrho > Trhomax ){ 
                        printf("density is outside of the table's range \n");
                        printf("density was set \n");
                        printf("rho = %6.40lf \n", BCrho);
                }
                table->data(1,BCntheta, BCnphi, BCnt ) = BCpress; // pressure - here I am using pressure as read from the table with pressure been previous interpolated by the eos
                if ( BCpress < Tpressmin or  BCpress > Tpressmax ){
                        printf("pressure is outside of the table's range \n");
                        printf("pressure was set \n");
                        printf("pressure = %6.40lf \n", BCpress);
                }
                table->data(2 , BCntheta, BCnphi, BCnt ) = BCvr;      // v_r
                if ( BCvr <  Tvrmin or  BCvr > Tvrmax ){
                        printf("vr is outside of the table's range \n");
                        printf("vr was set \n");
                        printf("vr = %640lf \n", BCvr);
                }
                table->data(3 , BCntheta, BCnphi, BCnt ) = BCvt;      // v_theta
                if ( BCvt <  Tvtmin or  BCvt >  Tvtmax ){
                        printf("vt is outside of the table's range \n");
                        printf("vt was set \n");
                        printf("vt = %6.40lf \n", BCvt);
                }
                table->data(4 , BCntheta, BCnphi, BCnt ) = BCvp;      // v_phi
                if ( BCvp <  Tvpmin or  BCvp >  Tvpmax ){
                        printf("vp is outside of the table's range \n");
                        printf("vp was set \n");
                        printf("vp = %6.40lf \n", BCvp);
                }
                table->data(5 , BCntheta, BCnphi, BCnt ) = BCinter;      // internal energy
                table->data(6 , BCntheta, BCnphi, BCnt ) = BCWlorentz;   // lorentz factor
                nline = nline+1;
        }
        printf("Successfuly read %d lines \n", nline); 
        printf("Filling ghost zones \n ");
        for (int kk = 0 ; kk <= BCNt ; kk++) {
                for (int ii =0 ; ii<=BCNtheta ; ii++){
                        for (int jj =0; jj <=BCNghts ; jj++){

                                // Region 1 to Region 1' (phi <0)
                                table->data(0 , BCNghts + ii, jj, kk ) = table->data(0 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //density
                                table->data(1 , BCNghts + ii, jj, kk ) = table->data(1 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //pressure
                                table->data(2 , BCNghts + ii, jj, kk ) = table->data(2 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //vr
                                table->data(3 , BCNghts + ii, jj, kk ) = table->data(3 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //vt
                                table->data(4 , BCNghts + ii, jj, kk ) = table->data(4 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //vp
                                table->data(5 , BCNghts + ii, jj, kk ) = table->data(5 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //inter
                                table->data(6 , BCNghts + ii, jj, kk ) = table->data(6 , BCNghts + ii, BCNphi + BCNghts - jj, kk ); //W

                                // Region 5 to Region 5' (phi > 2pi)

                                table->data(0 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(0 , BCNghts + ii, 2 * BCNghts - jj, kk ); //density
                                table->data(1 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(1 , BCNghts + ii, 2 * BCNghts - jj, kk ); //pressure
                                table->data(2 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(2 , BCNghts + ii, 2 * BCNghts - jj, kk ); //vr
                                table->data(3 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(3 , BCNghts + ii, 2 * BCNghts - jj, kk ); //vt
                                table->data(4 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(4 , BCNghts + ii, 2 * BCNghts - jj, kk ); //vp
                                table->data(5 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(5 , BCNghts + ii, 2 * BCNghts - jj, kk ); //inter
                                table->data(6 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(6 , BCNghts + ii, 2 * BCNghts - jj, kk ); //W

                        }
                }
                for (int ii =0 ; ii<=BCNghts ; ii++){
                        for (int jj =0; jj <=BCNphi ; jj++){
                                // Region 3 to Region 3' (theta <0)
                                table->data(0 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(0 , 2 * BCNghts - ii, BCNghts + jj, kk ); //density
                                table->data(1 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(1 , 2 * BCNghts - ii, BCNghts + jj, kk ); //pressure
                                table->data(2 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(2 , 2 * BCNghts - ii, BCNghts + jj, kk ); //vr
                                table->data(3 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) = - table->data(3 , 2 * BCNghts - ii, BCNghts + jj, kk ); //vt
                                table->data(4 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) = + table->data(4 , 2 * BCNghts - ii, BCNghts + jj, kk ); //vp 
                                table->data(5 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(5 , 2 * BCNghts - ii, BCNghts + jj, kk ); //inter  
                                table->data(6 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(6 , 2 * BCNghts - ii, BCNghts + jj, kk ); //W                               
                                // Region 7 to Region 7' (theta > pi)
                                table->data(0 , ii, BCNghts + jj, kk ) =   table->data(0 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //density
                                table->data(1 , ii, BCNghts + jj, kk ) =   table->data(1 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //pressure 
                                table->data(2 , ii, BCNghts + jj, kk ) =   table->data(2 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //vr   
                                table->data(3 , ii, BCNghts + jj, kk ) = - table->data(3 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //vt    
                                table->data(4 , ii, BCNghts + jj, kk ) = + table->data(4 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //vp     
                                table->data(5 , ii, BCNghts + jj, kk ) =   table->data(5 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //inter                
                                table->data(6 , ii, BCNghts + jj, kk ) =   table->data(6 , BCNtheta + BCNghts - ii, BCNghts + jj, kk ); //W                                             
                        }
                }

                for (int ii =0 ; ii<= BCNghts ; ii++){
                        for (int jj =0; jj <= BCNghts ; jj++){

                                // Region 2 to Region 2'
                                table->data(0 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(0 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //density         
                                table->data(1 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(1 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //pressure           
                                table->data(2 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(2 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //vr    
                                table->data(3 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(3 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //vt  
                                table->data(4 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(4 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //vp   
                                table->data(5 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(5 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //inter
				table->data(6 , BCNtheta + BCNghts + ii , jj , kk ) = table->data(6 , 2 * BCNghts - ii , BCNphi + BCNghts -jj , kk ); //W   

                                // Region 4 to Region 4'
                                table->data(0 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(0 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //density  
				table->data(1 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(1 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //pressure
                                table->data(2 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(2 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //vr       
		                table->data(3 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(3 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //vt
                                table->data(4 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(4 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //vp
                                table->data(5 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(5 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //inter
                                table->data(6 , BCNtheta + BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(6 , 2 * BCNghts - ii, 2 * BCNghts - jj, kk ); //W  

                                // Region 6 to Region 6'
                                table->data(0 , ii, BCNphi + BCNghts + jj, kk ) = table->data(0 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //density
                                table->data(1 , ii, BCNphi + BCNghts + jj, kk ) = table->data(1 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //pressure
                                table->data(2 , ii, BCNphi + BCNghts + jj, kk ) = table->data(2 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //vr
                                table->data(3 , ii, BCNphi + BCNghts + jj, kk ) = table->data(3 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //vt
				table->data(4 , ii, BCNphi + BCNghts + jj, kk ) = table->data(4 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //vp
			        table->data(5 , ii, BCNphi + BCNghts + jj, kk ) = table->data(5 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //inter
				table->data(6 , ii, BCNphi + BCNghts + jj, kk ) = table->data(6 , BCNtheta + BCNghts - ii, 2 * BCNghts - jj, kk ); //W 

                                // Region 8 to Region 8'
                                table->data(0 , ii, jj, kk ) = table->data(0 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //density
				table->data(1 , ii, jj, kk ) = table->data(1 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //pressure          
				table->data(2 , ii, jj, kk ) = table->data(2 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //vr
                                table->data(3 , ii, jj, kk ) = table->data(3 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //vt       
                                table->data(4 , ii, jj, kk ) = table->data(4 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //vp  
                                table->data(5 , ii, jj, kk ) = table->data(5 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //inter 
                                table->data(6 , ii, jj, kk ) = table->data(6 , BCNtheta + BCNghts - ii, BCNphi + BCNghts - jj, kk ); //W  
			}
		}
	}    
        printf("Done with filling ghost zones \n ");                                                                                                                                    
        printf("checking how complete is the data table \n");
        int NVAL = 0;
        int MVAL = 0;
        for (int kk = 0 ; kk <= BCNt ; kk++) {
              for (int ii =0 ; ii<= BCNtheta+2*BCNghts ; ii++){
	              for (int jj =0; jj <= BCNphi+2*BCNghts ; jj++){
                                if (table->data(0 , ii, jj, kk )==0 and table->data(1 , ii, jj, kk )==0 and table->data(2 , ii, jj, kk )==0 and table->data(3 , ii, jj, kk )==0 and table->data(4 , ii, jj, kk )==0 and table->data(5 , ii, jj, kk )==0 and table->data(6 , ii, jj, kk )==0) { 
					NVAL = NVAL+1;  
                                }
                                MVAL = MVAL +1;
                        }
                }
       }   
       printf("data table's length is = %d \n",MVAL);
       printf("data table has %d empty entries \n",NVAL);
       printf("Done with reading the data \n "); 
       return ;
}

//***********************************************************

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


Real IntEnProfile(const Real th, const Real ph, const Real tm) {
  Real IntEn;
  IntEn = table->interpolate(5, th, ph, tm);
  return IntEn ;
}

Real WProfile(const Real th, const Real ph, const Real tm) {
  Real W;
  W = table->interpolate(6, th, ph, tm);
  return W ;
}


//***********************************************************
void InnerBoundary(MeshBlock *pmb, Coordinates *pcoord,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku,
                  int ngh) {


        printf("Entering BC condition \n");

        float BCThetamaxval = BCThetamax + BCNghts*BCDtheta;
        float BCThetaminval = BCThetamin - BCNghts*BCDtheta;
        float BCPhimaxval   = BCPhimax   + BCNghts*BCDphi;
        float BCPhiminval   = BCPhimin   - BCNghts*BCDphi;
        float DeltaPhi      = BCNghts*BCDphi;
        float DeltaTheta    = BCNghts*BCDtheta;
        float BCTmaxval     = BCTmax + BCtin;

        for (int k=kl; k<=ku; ++k) {
                for (int j=jl; j<=ju; ++j){
                        for (int i=1; i<=ngh; ++i){

                                Real rval   = pcoord->x1v(il-i) ;
                                Real thval  = pcoord->x2v(j) ;
                                Real phival = pcoord->x3v(k) ;


                                float timeshift     = time   + BCtin;
                                if (timeshift < BCTmax ){
		                        if (thval > BCThetamaxval or thval < BCThetaminval){
                                                printf(" Theta outside of interpolation range  \n");
                                                printf("theta = %6.40lf \n", thval);
                                                printf("thetamin = %6.40lf \n", BCThetaminval );
                                                printf("thetamax = %6.40lf \n", BCThetamaxval );
                                                //return;
                                        }
                                        if (phival > BCPhimaxval or phival < BCPhiminval){
                                                printf(" Phi outside of interpolation range  \n");
                                                printf("phi = %6.40lf \n", phival);
                                                printf("phimin  = %6.40lf \n", BCPhiminval);
                                                printf("phimax  = %6.40lf \n", BCPhimaxval);
                                                //return;
                                        }
                                        if (timeshift > BCTmaxval or timeshift < 0 ){
                                                printf(" time outside of interpolation range  \n");
                                                printf("timeshift = %6.40lf \n", timeshift);
                                                printf("tmin = %6.40lf \n", 0);
                                                printf("tmax = %6.40lf \n", BCTmaxval);
                                                //return;
                                        }
                                        if (std::isfinite(DenProfile(thval, phival, timeshift)) == 1 ) {
                                                prim(IDN, k ,j ,il-i) = 1000000.0*DenProfile(thval, phival, timeshift);
                                                //printf("Density set to %6.40lf \n", prim(IDN, k ,j ,il-i));
                                        }
                                        else {
                                                printf("Error in setting density! \n");
//                                                return;
                                        }
                                        if (std::isfinite(PressProfile(thval, phival, timeshift)) == 1 ) {
                                                prim(IPR, k, j, il-i) = 100000000.0*PressProfile(thval, phival, timeshift);
                                        }
                                        else {
                                                printf("Error in setting pressure! \n");
//                                                return;
                                        }                                                   
                                        if (std::isfinite(IntEnProfile(thval, phival, timeshift)) == 1){
                                                prim(IEN, k, j, il-i) = DenProfile(thval, phival, timeshift)*IntEnProfile(thval, phival, timeshift);   
                                                if (RELATIVISTIC_DYNAMICS){  // this should only ever be SR with this file
                                                        prim(IEN, k, j, il-i) += DenProfile(thval, phival, timeshift);
                                                }
                                        }
                                        else{
                                                printf("Error in setting internal energy! \n");
                                        }
                                        if (std::isfinite(VRProfile(thval, phival, time + BCtin)) == 1 and std::isfinite(WProfile(thval, phival, time + BCtin))==1) {
                                                if (RELATIVISTIC_DYNAMICS){
                                                        prim(IVX, k, j, il-i) = 10.0; //WProfile(thval, phival, timeshift)*VRProfile(thval, phival, timeshift);
                                                        printf("vr rel. \n");
                                                        //printf("W = %6.40lf \n", WProfile(thval, phival, timeshift));
                                                }
                                                else {
                                                        prim(IVX, k, j, il-i) = VRProfile(thval, phival, timeshift);
                                                        printf("vr non rel. \n");
                                                        }

                                                //if (std::abs(prim(IVX, k, j, il-i)) >= 1.0 ) {
						//	printf("!!!!!!!!!!!!!!!!!!  Radial  Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                                //        printf("     r = %6.4lf \n", rval);
                                                //        printf(" theta = %6.4lf \n", thval);                                                       
						//	printf("phival = %6.4lf \n", phival);
                                                //        printf("  time = %6.4lf \n", timeshift);
                                                //        printf("    Vr = %6.40lf \n", VRProfile(thval, phival, timeshift));
                                                //        printf(" !!!!!!!!!!!!!!!!!!  Radial  Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                //        return;  
                                                //} 
                                        }
                                        else {
			                        printf("Error in setting radial velocity! \n");
                                                return;
                                        }
                                        if (std::isfinite(VTProfile(thval, phival, timeshift)) == 1 and std::isfinite(WProfile(thval, phival, time + BCtin))==1) {
                                                if (RELATIVISTIC_DYNAMICS){
	                                                prim(IVY, k, j, il-i) = 10.0 ; //WProfile(thval, phival, timeshift)*VTProfile(thval, phival, timeshift);
                                                        printf("vt rel. \n"); 
                                                 }

                                                else {
                                                        prim(IVY, k, j, il-i) = VTProfile(thval, phival, timeshift);
                                                        printf("vt non rel. \n");
                                                }  

                                                //if (std::abs(prim(IVY, k, j, il-i)) >= 1.0 ) {
                                                //        printf(" !!!!!!!!!!!!!!!!!! Theta Velocity is too big !!!!!!!!!!!!!!!!! \n "); 
                                                //        printf("     r = %6.4lf \n", rval);
                                                //        printf(" theta = %6.4lf \n", thval);
                                                //        printf("phival = %6.4lf \n", phival);
                                                //        printf("  time = %6.4lf \n", timeshift);
                                                //        printf("    Vt = %6.40lf \n", VTProfile(thval, phival, timeshift));
                                                //        printf(" !!!!!!!!!!!!!!!!!! Theta Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                //        break;
	                                       //} 
                                        }
                                        else {
                                                printf("Error in setting theta velocity! \n");
                                                break;
					}
                                        if (std::isfinite(VPProfile(thval, phival, timeshift)) == 1 and std::isfinite(WProfile(thval, phival, time + BCtin))==1) {
                                                if (RELATIVISTIC_DYNAMICS){
                                                        prim(IVZ, k, j, il-i) = 10.0; //* WProfile(thval, phival, timeshift)*VPProfile(thval, phival, timeshift);
                                                        printf("vp rel. \n");
         					}
                                                else {
                                                        prim(IVZ, k, j, il-i) = VPProfile(thval, phival, timeshift);
                                                        printf("vp non rel. \n");
                                                }
                                                //if (std::abs(prim(IVZ, k, j, il-i)) >= 1.0 ) {
                                                //        printf(" !!!!!!!!!!!!!!!!!! Phi Velocity is too big !!!!!!!!!!!!!!!!! \n "); 
                                                //        printf("     r = %6.4lf \n", rval); 
                                                //        printf(" theta = %6.4lf \n", thval);
                                                //        printf("phival = %6.4lf \n", phival); 
                                                //        printf("    Vp = %6.4lf \n", timeshift);
                                                //        printf("from interpolation = %6.40lf \n", VPProfile(thval, phival, timeshift));
                                                //        printf(" !!!!!!!!!!!!!!!!!! Phi Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                //        break;
                                                //}
                                        }
                                        else {
                                                printf("Error in setting phi velocity! \n");
                                                break;
                                        }
                                                                                                                            }
                                else{
                                        prim(IDN, k , j, il-i) = da ;
                                        prim(IPR, k , j, il-i) = pa ;
                                        prim(IVX, k , j, il-i) = 0 ;
                                        prim(IVY, k , j, il-i) = 0 ;
                                        prim(IVZ, k , j, il-i) = 0 ; 
                                        prim(IEN, k, j, il-i) = pa/(Gammaval-1);
                                        if (RELATIVISTIC_DYNAMICS){  // this should only ever be SR with this file
                                                prim(IEN, k, j, il-i) += da;
                                        }


                                                                                            
                                }
                        }
                }               
        }
        return;
}

} // namespace 
