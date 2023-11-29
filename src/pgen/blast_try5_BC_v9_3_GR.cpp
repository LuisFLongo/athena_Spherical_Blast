
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

// Comment : this version seats on top of blast_try5_BC_v9_3 but
//           including GR effects of a central mass with a given spin

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

// LFLM inclusion starts
// Configuration checking
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif
// LFLM inclusion ends


int RefinementCondition(MeshBlock *pmb);

// LFLM inclusion starts

Real threshold, method, da, pa;
float BCThetamin, BCPhimin, BCTmin, BCThetamax, BCPhimax, BCTmax, BCMtheta, BCMphi,BCMt , BCtin , BCDtheta, BCDphi, BCDt; //, BCMghts
Real Gammaval, Kval; 
float Trhomax, Tpressmax, Tvrmax, Tvtmax, Tvpmax;
float Trhomin, Tpressmin, Tvrmin, Tvtmin, Tvpmin;
int BCNtheta, BCNphi,BCNt, BCNghts;

Real m;                           // mass M of black hole
Real a;                           // spin of black hole (0 <= a < M)

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


void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                  Real *ptheta, Real *pphi);

void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3);

// LFLM inclusion ends

void Mesh::InitUserMeshData(ParameterInput *pin) {
//LFLM inclusion start
  pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  da   = pin->GetOrAddReal("problem", "damb", 1.0);
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

  //BCNghts = 0;
  BCNghts = std::ceil(std::max(std::max(Nthetadown, Nthetaup),std::max(Nphidown, Nphiup)));

  //printf("Type of BCNghts is: \n");
  std::cout << typeid(BCNghts).name() << '\n';

  //printf("Number of ghost set to = %d \n", BCNghts);
  float A1 = BCThetamax + BCNghts*BCDtheta;
  float A2 = BCThetamin - BCNghts*BCDtheta;
  //printf("New theta max = %6.40lf \n", A1 );
  //printf("New theta min = %6.40lf \n", A2 );

  float A3 = BCPhimax + BCNghts*BCDphi;
  float A4 = BCPhimin - BCNghts*BCDphi;
  //printf("New phi max = %6.40lf \n", A3 );
  //printf("New phi min = %6.40lf \n", A4 ); 

  //Gammaval = pin->GetReal("hydro","gamma");
  //Kval     = pin->GetReal("hydro","kappa");

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
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
    method = pin->GetReal("problem","mtd");
   
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

  Real bx = 0.0, by = 0.0, bz = 0.0;

  AthenaArray<Real> b, g, gi;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);
  g.NewAthenaArray(NMETRIC, ncells1);
  gi.NewAthenaArray(NMETRIC, ncells1);
  //printf("About to enter loop over coordinates\n");
  // Initialize hydro variables
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i=il; i<=iu; ++i) {

        Real rval = pcoord->x1v(i);
        Real tval = pcoord->x2v(j);
        Real pval = pcoord->x3v(k);

        Real r(0.0), theta(0.0), phi(0.0);
        GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j), pcoord->x3v(k), &r,
                                     &theta, &phi);

        Real ut , ur;
        Real u0(0.0), u1(0.0), u2(0.0), u3(0.0);
        TransformVector(ut, ur, 0.0, 0.0, rval, tval, pval, &u0, &u1, &u2, &u3);
        Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
        Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
        Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;


        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = da;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pa;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = uu1;//u1 - gi(I01,i)/gi(I00,i) * u0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uu2;//u2 - gi(I02,i)/gi(I00,i) * u0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uu3;//u3 - gi(I03,i)/gi(I00,i) * u0;

        // Calculate cell-centered magnetic fields given Minkowski values
        Real bcont = 0.0;
        Real bconx = bx;
        Real bcony = by;
        Real bconz = bz;
        Real bcon0, bcon1, bcon2, bcon3;
        TransformVector(bcont, bconx, bcony, bconz, rval, tval, pval, &bcon0, &bcon1, &bcon2,
                        &bcon3);
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);


  //printf("Exiting MeshBlock::ProblemGenerator \n");
}

// LFLM implementation of the reader for the outflow data start
// this should read the file
// import the data to an 3D AthenaArray, as in (var, theta, phi, t)


void OutflowReader(const char *filename, InterpTable3D *table) {

//first set some constants to translate the physical unities to code unities
//Some conversion factors to go between Lorene data and Athena data. Shamelessly stolen from
//the Einstein Toolkit's Mag_NS.cc.
//


    //    printf(" Entering Outflow Reader\n");

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
        // filling in the sides of the ghost zones along phi direction
        // peridiodical BC on phi
                for (int ii =0 ; ii<=BCNtheta ; ii++){
                        for (int jj =0; jj <=BCNghts ; jj++){

                                // Region 1 to Region 1' (phi <0)
                                table->data(0 , BCNghts + ii, jj, kk ) = table->data(0 , BCNghts + ii, BCNphi + jj, kk ); //density 
                                table->data(1 , BCNghts + ii, jj, kk ) = table->data(1 , BCNghts + ii, BCNphi + jj, kk ); //pressure
                                table->data(2 , BCNghts + ii, jj, kk ) = table->data(2 , BCNghts + ii, BCNphi + jj, kk ); //vr
                                table->data(3 , BCNghts + ii, jj, kk ) = table->data(3 , BCNghts + ii, BCNphi + jj, kk ); //vt
                                table->data(4 , BCNghts + ii, jj, kk ) = table->data(4 , BCNghts + ii, BCNphi + jj, kk ); //vp
                                table->data(5 , BCNghts + ii, jj, kk ) = table->data(5 , BCNghts + ii, BCNphi + jj, kk ); //inter    
                                table->data(6 , BCNghts + ii, jj, kk ) = table->data(6 , BCNghts + ii, BCNphi + jj, kk ); //W    
  

                                // Region 5 to Region 5' (phi > 2pi)
                                table->data(0 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(0 , BCNghts + ii, BCNghts + jj, kk ); //density 
                                table->data(1 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(1 , BCNghts + ii, BCNghts + jj, kk ); //pressure
                                table->data(2 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(2 , BCNghts + ii, BCNghts + jj, kk ); //vr
                                table->data(3 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(3 , BCNghts + ii, BCNghts + jj, kk ); //vt
                                table->data(4 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(4 , BCNghts + ii, BCNghts + jj, kk ); //vp
                                table->data(5 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(5 , BCNghts + ii, BCNghts + jj, kk ); //inter
                                table->data(6 , BCNghts + ii, BCNphi + BCNghts + jj, kk ) = table->data(6 , BCNghts + ii, BCNghts + jj, kk ); //W

                        }
                }
        // filling in the sides of the ghost zones along theta direction
        // parity BC on theta as in Table I of https://arxiv.org/pdf/2002.06225.pdf
                for (int ii =0 ; ii<=BCNghts ; ii++){
                        for (int jj =0; jj <=BCNphi ; jj++){
                                // Region 3 to Region 3' (theta <0)
				table->data(0 , ii, BCNghts + jj, kk ) =   table->data(0 , BCNtheta + ii, BCNghts + jj, kk ); //density
                                table->data(1 , ii, BCNghts + jj, kk ) =   table->data(1 , BCNtheta + ii, BCNghts + jj, kk ); //pressure
                                table->data(2 , ii, BCNghts + jj, kk ) =   table->data(2 , BCNtheta + ii, BCNghts + jj, kk ); //vr
                                table->data(3 , ii, BCNghts + jj, kk ) = - table->data(3 , BCNtheta + ii, BCNghts + jj, kk ); //vt
                                table->data(4 , ii, BCNghts + jj, kk ) = + table->data(4 , BCNtheta + ii, BCNghts + jj, kk ); //vp 
                                table->data(5 , ii, BCNghts + jj, kk ) =   table->data(5 , BCNtheta + ii, BCNghts + jj, kk ); //inter  
                                table->data(6 , ii, BCNghts + jj, kk ) =   table->data(6 , BCNtheta + ii, BCNghts + jj, kk ); //W  

                                // Region 7 to Region 7' (theta > pi)
                                table->data(0 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(0 , BCNghts + ii, BCNghts + jj, kk ); //density
                                table->data(1 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(1 , BCNghts + ii, BCNghts + jj, kk ); //pressure
                                table->data(2 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(2 , BCNghts + ii, BCNghts + jj, kk ); //vr
                                table->data(3 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) = - table->data(3 , BCNghts + ii, BCNghts + jj, kk ); //vt
                                table->data(4 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) = + table->data(4 , BCNghts + ii, BCNghts + jj, kk ); //vp 
                                table->data(5 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(5 , BCNghts + ii, BCNghts + jj, kk ); //inter
                                table->data(6 , BCNtheta + BCNghts + ii, BCNghts + jj, kk ) =   table->data(6 , BCNghts + ii, BCNghts + jj, kk ); //W 
                        }
                }
        // filling in the corners of the ghost zones
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
                                        //printf("*********************** \n");
                                        //printf("BCntheta = %d \n", ii);
                                        //printf("BCnphi   = %d \n", jj);
                                        //printf("BCnt     = %d \n", kk);
                                }
                                MVAL = MVAL +1;
                        }
                }
        }
//        printf("data table's length is = %d \n",MVAL);
//        printf("data table has %d empty entries \n",NVAL);


//	printf("Done with reading the data \n ");
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



// LFLM inplementation of the intepolation ends
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






// LFLM implementation of the inner boundary condition start

void InnerBoundary(MeshBlock *pmb, Coordinates *pcoord,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku,
                  int ngh) {

  //      printf("Entering BC condition \n");

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

				//printf("r     = %6.40lf \n",rval);
				//printf("theta = %6.40lf \n",thval);
				//printf("phi   = %6.40lf \n",phival);

                                float timeshift     = time   + BCtin;
				if (timeshift < BCTmax ){
	                                if (thval > BCThetamaxval or thval < BCThetaminval){
	                                        printf(" Theta outside of interpolation range  \n");
        	                                printf("theta = %6.40lf \n", thval);
                	                        printf("thetamin = %6.40lf \n", BCThetaminval );
                        	                printf("thetamax = %6.40lf \n", BCThetamaxval );
                                	        return;
	                                }
        	                        if (phival > BCPhimaxval or phival < BCPhiminval){
                	                        printf(" Phi outside of interpolation range  \n");
                        	                printf("phi = %6.40lf \n", phival);
                        	                printf("phimin  = %6.40lf \n", BCPhiminval);
                        	                printf("phimax  = %6.40lf \n", BCPhimaxval);
                        	                return;
	                                }
	                                if (timeshift > BCTmaxval or timeshift < 0 ){
	                                        printf(" time outside of interpolation range  \n");
	                                        printf("timeshift = %6.40lf \n", timeshift);
	       	                                printf("tmin = %6.40lf \n", 0);
	                                        printf("tmax = %6.40lf \n", BCTmaxval);
	                                        return;
	                                }
        	                        if (std::isfinite(DenProfile(thval, phival, timeshift)) == 1 ) {
	                                        prim(IDN, k ,j ,il-i) = DenProfile(thval, phival, timeshift);
						//printf("Density set to %6.40lf \n", prim(IDN, k ,j ,il-i));
                	                }
                        	        else {
                                	        printf("Error in setting density! \n");
                                        	return;
	                            	}
	                                if (std::isfinite(PressProfile(thval, phival, timeshift)) == 1 ) {
	                                        prim(IPR, k, j, il-i) = PressProfile(thval, phival, timeshift);
        	                        }
	                                else {
        	                                printf("Error in setting pressure! \n");
                	                        return;
                        	        }

//					if (std::isfinite(IntEnProfile(thval, phival, timeshift)) == 1){
//						prim(IEN, k, j, il-i) = DenProfile(thval, phival, timeshift)*IntEnProfile(thval, phival, timeshift);
//                                        	if (RELATIVISTIC_DYNAMICS){  // this should only ever be SR with this file
//                                                	prim(IEN, k, j, il-i) += DenProfile(thval, phival, timeshift);
//                                        	}
//					}
//					else{
//						printf("Error in setting internal energy! \n");
//					}
                                	if (std::isfinite(VRProfile(thval, phival, time + BCtin)) == 1 and std::isfinite(WProfile(thval, phival, time + BCtin))==1) {
        					if (RELATIVISTIC_DYNAMICS){
							prim(IVX, k, j, il-i) = WProfile(thval, phival, timeshift)*VRProfile(thval, phival, timeshift); 
                                                        //printf("vr rel. \n");
							//printf("W = %6.40lf \n", WProfile(thval, phival, timeshift));
						}
						else {	
							prim(IVX, k, j, il-i) = VRProfile(thval, phival, timeshift);
                                                        //printf("vr non rel. \n");
						}
                                        	if (std::abs(prim(IVX, k, j, il-i)) >= 1.0 ) {
                                                	printf("!!!!!!!!!!!!!!!!!!  Radial  Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                                	printf("     r = %6.4lf \n", rval);
                                                	printf(" theta = %6.4lf \n", thval);
                                                	printf("phival = %6.4lf \n", phival);
                                                	printf("  time = %6.4lf \n", timeshift);
                                                	printf("    Vr = %6.40lf \n", VRProfile(thval, phival, timeshift));    
                                                	printf(" !!!!!!!!!!!!!!!!!!  Radial  Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                	return;
                                        	}
                                	}
                                	else {
                                        	printf("Error in setting radial velocity! \n");
                                        	return;
                                	}
                                	if (std::isfinite(VTProfile(thval, phival, timeshift)) == 1 and std::isfinite(WProfile(thval, phival, time + BCtin))==1) {
                                        	
                                                if (RELATIVISTIC_DYNAMICS){
                                                        prim(IVY, k, j, il-i) = WProfile(thval, phival, timeshift)*VTProfile(thval, phival, timeshift); 
                                                	//printf("vt rel. \n");
						}
                                                else {
                                                        prim(IVY, k, j, il-i) = VTProfile(thval, phival, timeshift);
                                                        //printf("vt non rel. \n");

                                                }
	                                  	if (std::abs(prim(IVY, k, j, il-i)) >= 1.0 ) {
                                                	printf(" !!!!!!!!!!!!!!!!!! Theta Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                                	printf("     r = %6.4lf \n", rval);
                                                	printf(" theta = %6.4lf \n", thval);
                                                	printf("phival = %6.4lf \n", phival);
                                                	printf("  time = %6.4lf \n", timeshift);
                                                	printf("    Vt = %6.40lf \n", VTProfile(thval, phival, timeshift));
                                                	printf(" !!!!!!!!!!!!!!!!!! Theta Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                	break;
                                        	}
                                	}
                                	else {
                                        	printf("Error in setting theta velocity! \n");
                                        	break;
                                	}
                                	if (std::isfinite(VPProfile(thval, phival, timeshift)) == 1 and std::isfinite(WProfile(thval, phival, time + BCtin))==1) {
                                                if (RELATIVISTIC_DYNAMICS){
                                                        prim(IVZ, k, j, il-i) = WProfile(thval, phival, timeshift)*VPProfile(thval, phival, timeshift); 
                                                        //printf("vp rel. \n");

						}
                                                else {
                                                        prim(IVZ, k, j, il-i) = VPProfile(thval, phival, timeshift);
                                                        //printf("vp non rel. \n");
						}
                                        	if (std::abs(prim(IVZ, k, j, il-i)) >= 1.0 ) {
                                                	printf(" !!!!!!!!!!!!!!!!!! Phi Velocity is too big !!!!!!!!!!!!!!!!! \n ");
                                                	printf("     r = %6.4lf \n", rval);
                                                	printf(" theta = %6.4lf \n", thval);
                                                	printf("phival = %6.4lf \n", phival);
                                                	printf("    Vp = %6.4lf \n", timeshift);
                                                	printf("from interpolation = %6.40lf \n", VPProfile(thval, phival, timeshift));
                                                	printf(" !!!!!!!!!!!!!!!!!! Phi Velocity is too big !!!!!!!!!!!!!!!!! \n");
                                                	break;
                                        	}	

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
//	       				prim(IEN, k, j, il-i) = pa/(Gammaval-1);
//        				if (RELATIVISTIC_DYNAMICS){  // this should only ever be SR with this file
//          					prim(IEN, k, j, il-i) += da;
//        				}


				}	
                        }
                }
        }
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


