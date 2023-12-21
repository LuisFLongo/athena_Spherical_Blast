#ifndef UTILS_INTERP_TABLE3D_HPP_
#define UTILS_INTERP_TABLE3D_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file interp_table.hpp
//  \brief defines class InterpTable2D
//  Contains functions that implement an intpolated lookup table

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray


//3D
class InterpTable3D {
 public:
  InterpTable3D() = default; 
  InterpTable3D(const int mvar, const int my3, const int my2, const int my1);

  void SetSize(const int mvar, const int my3, const int my2, const int my1);
  Real interpolate(int mvar, Real y3, Real y2, Real y1);
  int nvar();
  AthenaArray<Real> data;
  void SetY1lim(Real Y1min, Real Y1max);
  void SetY2lim(Real Y2min, Real Y2max);
  void SetY3lim(Real Y3min, Real y3max);
  
  void GetY1lim(Real &y1min, Real &y1max);
  void GetY2lim(Real &y2min, Real &y2max);
  void GetY3lim(Real &y3min, Real &y3max);
  
  void GetSize(int &mvar, int &my3, int &my2, int &my1);

 private:
  int mvar_;
  int my1_;
  int my2_;
  int my3_;
  Real y1min_;
  Real y1max_;
  Real y1norm_;
  Real y2min_;
  Real y2max_;
  Real y2norm_;
  Real y3min_;
  Real y3max_;
  Real y3norm_;
};



#endif //UTILS_INTERP_TABLE3D_HPP_
