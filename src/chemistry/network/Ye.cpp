//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file H2.cpp
//! \brief implementation of functions in class ChemNetwork, using the simple
//! network for H2 formation and destruction.


// this class header
#include "Ye.hpp"

// C headers

// C++ header
#include <iostream>   // endl
#include <limits>     // inf
#include <sstream>    // stringstream
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>

// Athena++ header
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../../units/units.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../utils/thermo.hpp"
#include "network.hpp"
#include "../../utils/interp_table3D.hpp" //LFLM inclusion


// constants
//const Real ChemNetwork::kgr_ = 3e-17;

// species names
const std::array<std::string, NSPECIES> ChemNetwork::species_names = {"YE", "YE2"};

const int ChemNetwork::iYE_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "YE");

const int ChemNetwork::iYE2_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "YE2");


InterpTable3D *tableLATER;
InterpTable3D *tableEARLY;

// flag for Cv
//static bool is_const_Cv;


//----------------------------------------------------------------------------------------
//! \brief ChemNetwork constructor

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
  // number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
  pmy_mb_ = pmb;


  std::string fnlater;
  std::string fnearly;
  fnlater = pin->GetString("chemistry","eps_later_table");
  fnearly = pin->GetString("chemistry","eps_early_table");
//  // set the parameters from input file
  ReationRate_ = pin->GetOrAddReal("chemistry", "reation_rate", 2e-16);
//  kcr_ = xi_cr_ * 3.;
//  // set Cv: constant or H2 abundance dependent
//  is_const_Cv = pin->GetOrAddBoolean("problem", "is_const_Cv", true);
}

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork destructor

ChemNetwork::~ChemNetwork() {
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::InitializeNextStep(const int k, const int j, const int i)
//! \brief Set the rates of chemical reactions, eg. through density and radiation field.
//!
//! k, j, i are the corresponding index of the grid

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, nb, mb , mbMeV, vr, taucgs;
  // density
  rho = pmy_mb_->phydro->u(IDN, k, j, i);
  Ye_  = pmy_mb_->pscalars->r(0, k, j, i);

  vr = pmy_mb_->phydro->u(IVX, k, j, i); 

  taucgs = 1e-3/vr ; // expansion time scale in seconds 
	 	     // first approximationg of the homologous expansion given at MNRAS 478, 3298-3334 (2018) 
	             // private comunication with Fabio Magistrelly 
 		     // should I implement the full expression?
 		     // notice also thar vr should be in units of c, which is assumed when passing the outflow datas to the BC

  Real const c_light  = 299792458.0; // Speed of light [m/s]
  Real const G_grav   = 6.67428e-11; // Gravitational constant [m^3/kg/s^2]
  Real const M_sun    = 1.98892e+30; // Solar mass [kg]
  Real const e_charge = 1.60218e-19; // electron charge [C] (retrived from https://pdg.lbl.gov/2023/reviews/contents_sports.html)

  S_   = 50.0 ; //fix entropy per baryon for now I will test the code and then elaborate on the entropy
	       // questions : 1) which kind of process? adiabatic? isothermal?
	       // 	      2) which eos? polytopic? 
	       // 	      3) I athena already calculating entropy internaly that I can pass?
	       // 	      4) If not, should I implement a given analytic expression or a table?

  Tau_ = taucgs * std::pow(c_light,3) / ( M_sun  * G_grav ) ; // tau to code units

  mbMeV = 939.56542052 ; //MeV/c^2 (retrived from https://pdg.lbl.gov/2023/reviews/contents_sports.html)
  mb    = mbMeV * 1e6 * e_charge /(c_light*c_light*M_sun) ;// baryon mass in code unities (I am using the neutron mass for now)
  nb  = rho / mb ; //baryon number density 

  ne_ = nb * Ye_ ;  //electron number density

  // apply density floor
  //rho_floor = pmy_mb_->peos->GetDensityFloor();
  //rho = (rho > rho_floor) ?  rho : rho_floor;
  // hydrogen atom number density
  //nH_ =  rho;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::RHS(const Real t, const Real *y, const Real ED,
//!                       Real *ydot)
//! \brief RHS: right-hand-side of ODE.
//!
//! dy/dt = ydot(t, y). Here y are the abundance
//! of species. details see CVODE package documentation.
//! all input/output variables are in code units

void ChemNetwork::RHS(const Real t, const Real *y, const Real ED, Real *ydot) {
  const Real rate = ReationRate_ * y[iYE_];
  ydot[iYE_] = rate;
  for (int i=0; i<NSPECIES; i++) {
    // return in code units
    ydot[i] *= pmy_mb_->pmy_mesh->punit->code_time_cgs; //not sure what is done here
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED)
//! \brief energy equation dED/dt
//!
//! all input/output variables are in code units (ED is the energy density)

Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED) {
  // isothermal


  Real heating = EPSsmooth( Tau_, S_, Ye_, t);
  Real dEDdt  = ne_ * heating;


  return dEDdt;
}


//LFLM inclusion starts

//LFLM development of reading the knec tables
//
//reading eatly time table

void EARLYepsdatafitReader(const char *filenameearly , InterpTable3D *tableEARLY){

// Reading to feed to the fit of A(tau, s, ye). t0(tau, s, ye) and sigma(tau,s, ye)
// for the construction of the fit for the heating rate
// the is heating rate = A * (0.5-arctan((t-t_0)/sigma)/pi)^alpha
// https://arxiv.org/pdf/2111.06870.pdf


	std::ifstream finearly;
        finearly.open("filenameearly");
        const int MAX_CHARS_PER_LINE  = 512;
        const int MAX_TOKENS_PER_LINE = 5;
        const char* const DELIMITER = "  ";
        int nline = 0;

        if (!finearly.good()){
                printf("Outflow file not found \n");
                return;
        }

        int Ntau = 15;
        int Ns   = 21;
        int Nye  = 24;

        tableEARLY->SetSize( 4 , Ntau , Ns  , Nye  );

	Real taumin = 1.360; //ms 
	Real taumax = 100.0; //ms

	Real smin   = 1.820; //k_B/baryon
	Real smax   = 100.0; //k_B/baryon

	Real yemin  = 0.02 ; //dimensionless
	Real yemax  = 0.48 ; //dimensionless

        tableEARLY -> SetY3lim( std::log10(taumin) , std::log10(taumax) ); // log10 of tau, as the grid is evenly spaced in log scale
        tableEARLY -> SetY2lim(  std::log10(smin)  ,  std::log10(smax)  ); // log10 of entropy, as the grid is evenly spaced in log scale
        tableEARLY -> SetY1lim(  yemin ,  yemax );

	float Dtau = 1.2/9.0 ; // in log scale
	float Ds   = 2.0/23.0; // in log scale
	float Dye  = 0.02    ; // in linear scale

        float Token[MAX_TOKENS_PER_LINE] ;

        while (!finearly.eof()){
                char buf[MAX_CHARS_PER_LINE];
                finearly.getline(buf, MAX_CHARS_PER_LINE);
                int n = 0; // a for-loop index
                const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
                token[0] = strtok(buf, DELIMITER); // first token
                if (token[0]){ // zero if line is blank
                        for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
                                token[n] = strtok(0, DELIMITER); // subsequent tokens
                                Token[n] = std::stod(token[n]);
                        }
                }

		float TAU   = std::log10(Token[0]); // log10 of tau, as the grid is evenly spaced in log scale
		float ENT   = std::log10(Token[1]); // log10 of entropy, as the grid is evenly spaced in log scale
        	float YE    = Token[2];             // ye, as the grid is evenly spaced in linear scale
        	float A     = Token[3];             //
        	float T0    = Token[4];             //
        	float SIGMA = Token[5];             //06-07

	        int Ntau = round( (TAU - std::log10(taumin) ) / Dtau );
        	int Ns   = round( (ENT - std::log10(smin)   ) / Ds   );
        	int Nye     = round( (YE  - yemin )           / Dye  );


		tableEARLY->data(0 , Ntau , Ns , Nye ) = A;     //
        	tableEARLY->data(1 , Ntau , Ns , Nye ) = T0;    //
        	tableEARLY->data(2 , Ntau , Ns , Nye ) = SIGMA; //
	}



	return;
}



void LATERepsdatafitReader(const char *filenamelater , InterpTable3D *tableLATER){

// Reading to feed to the fit of A(tau, s, ye) and alpha(tau,s, ye)
// for the construction of the fit for the heating rate
// the is heating rate = A * (t )^(-alpha)
// https://arxiv.org/pdf/2111.06870.pdf

        std::ifstream finlater;
        finlater.open(filenamelater);
        const int MAX_CHARS_PER_LINE  = 512;
        const int MAX_TOKENS_PER_LINE = 5;
        const char* const DELIMITER = "  ";
        int nline = 0;

        if (!finlater.good()){
                printf("Outflow file not found \n");
                return;
        }

        int Ntau = 15;
        int Ns   = 21;
        int Nye  = 24;

        tableLATER->SetSize( 3 , Ntau , Ns  , Nye  );

        Real taumin = 1.360; //ms 
        Real taumax = 100.0; //ms

        Real smin   = 1.820; //k_B/baryon
        Real smax   = 100.0; //k_B/baryon

        Real yemin  = 0.02 ; //dimensionless
        Real yemax  = 0.48 ; //dimensionless

        tableLATER -> SetY3lim( std::log10(taumin) , std::log10(taumax) ); // log10 of tau, as the grid is evenly spaced in log scale
        tableLATER -> SetY2lim(  std::log10(smin)  ,  std::log10(smax)  ); // log10 of entropy, as the grid is evenly spaced in log scale
        tableLATER -> SetY1lim(  yemin ,  yemax );

        float Dtau = 1.2/9.0 ; // in log scale
        float Ds   = 2.0/23.0; // in log scale
        float Dye  = 0.02    ; // in linear scale


        float Token[MAX_TOKENS_PER_LINE] ;

        while (!finlater.eof()){
                char buf[MAX_CHARS_PER_LINE];
                finlater.getline(buf, MAX_CHARS_PER_LINE);
                int n = 0; // a for-loop index
                const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
                token[0] = strtok(buf, DELIMITER); // first token
                if (token[0]){ // zero if line is blank
                        for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
                                token[n] = strtok(0, DELIMITER); // subsequent tokens
                                Token[n] = std::stod(token[n]);
                        }
                }

                float TAU   = std::log10(Token[0]); // log10 of tau, as the grid is evenly spaced in log scale
                float ENT   = std::log10(Token[1]); // log10 of entropy, as the grid is evenly spaced in log scale
                float YE    = Token[2];             // ye, as the grid is evenly spaced in linear scale
                float A     = Token[3];             //
                float ALPHA = Token[4];             //

                int Ntau = round( (TAU - std::log10(taumin) ) / Dtau );
                int Ns   = round( (ENT - std::log10(smin)   ) / Ds   );
                int Nye  = round( (YE  - yemin )              / Dye  );


                tableLATER->data(0 , Ntau , Ns , Nye ) = A;     //
                tableLATER->data(1 , Ntau , Ns , Nye ) = ALPHA;    //

        }



        return;
}




Real EPSearly(const Real tauval, const Real sval, const Real yeval, const Real T) {
                                
// Fit for the heating rate at early times                                
//  A * (0.5-arctan((t-t_0)/sigma)/pi)^alpha                    
// https://arxiv.org/pdf/2111.06870.pdf (eq 2)


  Real T0val;
  Real Aval;
  Real SIGMAval;

  Aval     = tableEARLY->interpolate(0, std::log10(tauval), std::log10(sval), yeval); // remember that the table is evenly spaced on log space for 
  T0val    = tableEARLY->interpolate(1, std::log10(tauval), std::log10(sval), yeval); // tau and sval , I still have to check the units of the intepolation (cgs)
  SIGMAval = tableEARLY->interpolate(2, std::log10(tauval), std::log10(sval), yeval); // and the code 
  
  Real EPSval = Aval * ( 0.5 - std::atan( (T - T0val) / SIGMAval )/M_PI );


  return EPSval ;
}


Real EPSlater(const Real tauval, const Real sval, const Real yeval, const Real T) {

// Fit for the heating rate at late times
// A * t ** (-alpha)
// https://arxiv.org/pdf/2111.06870.pdf (eq 3)

  Real Aval;
  Real ALPHAval;

  Aval     = tableLATER->interpolate(0, std::log10(tauval), std::log10(sval), yeval); // remember that the table is evenly spaced on log space for
  ALPHAval = tableLATER->interpolate(1, std::log10(tauval), std::log10(sval), yeval); // tau and sval , I still have to check the units of the intepolation (cgs)
                                                                                      // and the code

  Real EPSval = Aval * std::pow( T , - ALPHAval );


  return EPSval ;
}


Real EPSsmooth(const Real tauval, const Real sval, const Real yeval, const Real T) {
// Every other function of the fit are assuming cgs units
// Here I have to convert all the input values from code to cgs
// the final answer I should convert from cgs to code units


        Real const c_light = 29979245800.0; // Speed of light [cm/s]
        Real const G_grav = 6.67428e-8; // Gravitational constant [cm^3/g/s^2]
        Real const M_sun = 1.98892e+33; // Solar mass [g]

        Real const ergperspergtocode = ( M_sun * G_grav )/ std::pow(c_light,5) ;
        Real const stocode           = std::pow(c_light,3) / ( M_sun * G_grav ) ;

        Real Tcgs = T / stocode ;

        Real T1 = 1.0*std::pow(10,3); // in cgs  see below equation (3) in https://arxiv.org/pdf/2111.06870.pdf
        Real T2 = 4.0*std::pow(10,4); // in cgs

	Real taucgs = tauval / stocode;

	Real EPSval ;
	if (T < T1){
		EPSval = EPSearly( taucgs, sval, yeval, T) ;
	}
	if (T > T2 ){
		EPSval = EPSlater( taucgs, sval, yeval, T) ;
	}
	else {

        	Real x  = (T-T2)/(T1-T2);
		Real EPSearlyval = EPSearly( taucgs, sval, yeval, T) ;
		Real EPSlaterval = EPSlater( taucgs, sval, yeval, T) ;
		
		EPSval = std::pow( 10 , (1-x)*std::log10( EPSearlyval ) + x*std::log10( EPSlaterval ) );
 
	}
	Real EPScode = EPSval / ergperspergtocode ;
	return EPScode ; 

}



//LFLM inclusion ends
