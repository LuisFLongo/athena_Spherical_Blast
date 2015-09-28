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
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file integrators.cpp
//  \brief 
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../fluid.hpp"
#include "../../mesh.hpp"
#include "../../parameter_input.hpp" 

// this class header
#include "fluid_integrator.hpp"

// constructor

HydroIntegrator::HydroIntegrator(Hydro *pf, ParameterInput *pin)
{
  pmy_fluid = pf;

// Allocate memory for scratch vectors

  int nthreads = pf->pmy_block->pmy_mesh->GetNumMeshThreads();
  int ncells1 = pf->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pf->pmy_block->block_size.nx2 + 2*(NGHOST);

  wl_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  wr_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  flx_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  jflx_j_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  kflx_k_.NewAthenaArray(nthreads,(NWAVE),ncells2,ncells1);
  face_area_.NewAthenaArray(nthreads,ncells1);
  face_area_p1_.NewAthenaArray(nthreads,ncells1);
  cell_volume_.NewAthenaArray(nthreads,ncells1);
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS)  // only used in (SR/GR)MHD
  {
    bb_normal_.NewAthenaArray(ncells1);
    lambdas_p_l_.NewAthenaArray(ncells1);
    lambdas_m_l_.NewAthenaArray(ncells1);
    lambdas_p_r_.NewAthenaArray(ncells1);
    lambdas_m_r_.NewAthenaArray(ncells1);
  }
  if (GENERAL_RELATIVITY)  // only used in GR
  {
    g_.NewAthenaArray(NMETRIC,ncells1);
    gi_.NewAthenaArray(NMETRIC,ncells1);
    cons_.NewAthenaArray(NWAVE,ncells1);
  }
}

// destructor

HydroIntegrator::~HydroIntegrator()
{
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
  jflx_j_.DeleteAthenaArray();
  kflx_k_.DeleteAthenaArray();
  face_area_.DeleteAthenaArray();
  face_area_p1_.DeleteAthenaArray();
  cell_volume_.DeleteAthenaArray();
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS)  // only used in (SR/GR)MHD
  {
    bb_normal_.DeleteAthenaArray();
    lambdas_p_l_.DeleteAthenaArray();
    lambdas_m_l_.DeleteAthenaArray();
    lambdas_p_r_.DeleteAthenaArray();
    lambdas_m_r_.DeleteAthenaArray();
  }
  if (GENERAL_RELATIVITY)
  {
    g_.DeleteAthenaArray();
    gi_.DeleteAthenaArray();
    cons_.DeleteAthenaArray();
  }
}
