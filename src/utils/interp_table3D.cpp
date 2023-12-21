//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file interp_table.cpp
//  \brief implements functions in class InterpTable2D an intpolated lookup table

// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../coordinates/coordinates.hpp" // Coordinates
#include "interp_table3D.hpp"

//LFLM inclusion starts
//
//3D
InterpTable3D::InterpTable3D(const int mvar, const int my3, const int my2, const int my1) {
  SetSize(mvar, my3, my2, my1);
}


void InterpTable3D::SetSize(const int mvar, const int my3, const int my2, const int my1) {
  //printf("Entering  InterpTable3D SetSize \n");
  mvar_ = mvar; // number of variables/tables
  //printf("mvar = %d", mvar);
  my3_ = my3;
  //printf("my3 = %d", my3);
  my2_ = my2; // slower indexing dimension
  //printf("my2 = %d", my2);
  my1_ = my1; // faster indexing dimension
  //printf("mx1 = %d", my1);
  //printf("Read all the variables \n");
  data.NewAthenaArray(mvar, my3, my2, my1);
  //printf("Exiting  InterpTable3D SetSize \n");
}

void InterpTable3D::SetY1lim(Real y1min, Real y1max) {
  y1min_ = y1min;
  y1max_ = y1max;
  y1norm_ = (my1_ - 1) / (y1max - y1min);
}

void InterpTable3D::SetY2lim(Real y2min, Real y2max) {
  y2min_ = y2min;
  y2max_ = y2max;
  y2norm_ = (my2_ - 1) / (y2max - y2min);
}

void InterpTable3D::SetY3lim(Real y3min, Real y3max) {
  y3min_ = y3min;
  y3max_ = y3max;
  y3norm_ = (my3_ - 1) / (y3max - y3min);
}


void InterpTable3D::GetY1lim(Real &y1min, Real &y1max) {
  y1min = y1min_;
  y1max = y1max_;
}


void InterpTable3D::GetY2lim(Real &y2min, Real &y2max) {
  y2min = y2min_;
  y2max = y2max_;
}

void InterpTable3D::GetY3lim(Real &y3min, Real &y3max) {
  y3min = y3min_;
  y3max = y3max_;
}


void InterpTable3D::GetSize(int &mvar, int &my3, int &my2, int &my1) {
  mvar = mvar_;
  my3 = my3_;
  my2 = my2_;
  my1 = my1_;
}


//Trilinear interpolation
Real InterpTable3D::interpolate(int var, Real y3, Real y2, Real y1) {
  Real x, y, z, xrl, yrl, zrl, out;
  z = (y3 - y3min_) * y3norm_;
  y = (y2 - y2min_) * y2norm_;
  x = (y1 - y1min_) * y1norm_;
  int zil = static_cast<int>(z);
  int yil = static_cast<int>(y); // lower y index
  int xil = static_cast<int>(x); // lower x index
  int mz = my3_;
  int my = my2_;
  int mx = my1_;
  //if off table, do linear extrapolation
  if (xil < 0) { // below xmin
    xil = 0;
  } else if (xil >= mx - 1) { // above xmax
    xil = mx - 2;
  }
  xrl = 1 + xil - x;  // x residual

  if (yil < 0) { // below ymin
    yil = 0;
  } else if (yil >= my - 1) { // above ymax
    yil = my - 2;
  }
  yrl = 1 + yil - y;  // y residual

  if (zil < 0) { // below zmin
    zil = 0;
  } else if (zil >= mz - 1) { // above zmax
    zil = mz - 2;
  }
  zrl = 1 + zil - z;  // z residual
  //if ( data(var, zil  , yil  , xil  ) == 0 or data(var, zil  , yil  , xil+1) == 0 or data(var, zil  , yil+1, xil  ) == 0 or data(var, zil  , yil+1, xil+1) == 0 or data(var, zil+1, yil  , xil  ) == 0 or data(var, zil+1, yil  , xil+1) == 0 or data(var, zil+1, yil+1, xil  ) == 0 or  data(var, zil+1, yil+1, xil+1) == 0){
//	printf("Interpolation table is possibly not set correctly: there are exact zeros on the tables \n");
	//printf(" zval = %6.40lf \n", zil);
	//printf(" yval = %6.40lf \n", yil);
	//printf(" xval = %6.40lf \n", xil);
//  }
//  else if ( std::isfinite(data(var, zil  , yil  , xil  )) !=1  or std::isfinite(data(var, zil  , yil  , xil+1))!=1 or std::isfinite(data(var, zil  , yil+1, xil  )) !=1 or std::isfinite(data(var, zil  , yil+1, xil+1)) !=1 or std::isfinite(data(var, zil+1, yil  , xil  )) !=1 or std::isfinite(data(var, zil+1, yil  , xil+1)) !=1 or std::isfinite(data(var, zil+1, yil+1, xil  )) != 1 or  std::isfinite(data(var, zil+1, yil+1, xil+1)) !=1){
//        printf("Interpolation table is possibly not set correctly: there are non numerical tables  on the tables \n");
	//printf(" zval = %6.40lf \n", zil);
	//printf(" yval = %6.40lf \n", yil);
	//printf(" xval = %6.40lf \n", xil);
//  }
//  else {printf("Interpolation if fine in this point \n");}

  //Sample from the 8 nearest data points and weight appropriately
  //data is an attribute of the eos class
   out =   xrl  *  yrl  *  zrl  * data(var, zil  , yil  , xil  )
        +(1-xrl)*  yrl  *  zrl  * data(var, zil  , yil  , xil+1)
        +  xrl  *(1-yrl)*  zrl  * data(var, zil  , yil+1, xil  )
        +(1-xrl)*(1-yrl)*  zrl  * data(var, zil  , yil+1, xil+1)
        +  xrl  *  yrl  *(1-zrl)* data(var, zil+1, yil  , xil  )
        +(1-xrl)*  yrl  *(1-zrl)* data(var, zil+1, yil  , xil+1)
        +  xrl  *(1-yrl)*(1-zrl)* data(var, zil+1, yil+1, xil  )
        +(1-xrl)*(1-yrl)*(1-zrl)* data(var, zil+1, yil+1, xil+1);
       
 return out;
}

