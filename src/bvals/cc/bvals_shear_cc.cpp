//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_cc.cpp
//  \brief functions that apply shearing box BCs for cell-centered variables
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals_interfaces.hpp"
#include "../bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//--------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadShearing(Real *buf, int nb)
//  \brief Load shearing box hydro boundary buffers

void CellCenteredBoundaryVariable::LoadShearing(Real *buf, int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int nx2 = pmb->block_size.nx2 - NGHOST;

  si = pmb->is - NGHOST; ei = pmb->is - 1;
  sk = pmb->ks;        ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1)  ek += NGHOST, sk -= NGHOST;
  // nb=0-3 for inner boundary; nb=4-7 for outer boundary
  switch (nb) {
    case 0:
      sj = pmb->je - joverlap_ - (NGHOST - 1); ej = pmb->je;
      if (joverlap_ > nx2) sj = pmb->js;
      break;
    case 1:
      sj = pmb->js; ej = pmb->je - joverlap_ + NGHOST;
      if (joverlap_ < NGHOST) ej = pmb->je;
      break;
    case 2:
      sj = pmb->je - (NGHOST - 1); ej = pmb->je;
      if (joverlap_ > nx2) sj = pmb->je - (joverlap_ - nx2) + 1;
      break;
    case 3:
      sj = pmb->js; ej = pmb->js + (NGHOST - 1);
      if (joverlap_ < NGHOST) ej = pmb->js + (NGHOST - joverlap_) - 1;
      break;
    case 4:
      sj = pmb->js; ej = pmb->js + joverlap_ + NGHOST - 1;
      if (joverlap_ > nx2) ej = pmb->je;
      break;
    case 5:
      sj = pmb->js + joverlap_ - NGHOST; ej = pmb->je;
      if (joverlap_ < NGHOST) sj = pmb->js;
      break;
    case 6:
      sj = pmb->js; ej = pmb->js + (NGHOST - 1);
      if (joverlap_ > nx2) ej = pmb->js + (joverlap_ - nx2) - 1;
      break;
    case 7:
      sj = pmb->je - (NGHOST - 1); ej = pmb->je;
      if (joverlap_ < NGHOST) sj = pmb->je - (NGHOST - joverlap_) + 1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in CellCenteredBoundaryVariable:LoadShearing "
          << std::endl << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }
  int p = 0;
  BufferUtility::PackData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffersForInit()
//  \brief Send shearing box boundary buffers for hydro variables

void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffersForInit() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;

  // KGF: hidden assumption that 2D?
  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  Real qomL = qshear_*Omega_0_*x1size_;

  if (shbb_.inner == true) {
    // step 1. -- add shear to the inner periodic boundary values
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          // add shear to conservative
          shboxvar_inner_hydro_(IM2,k,j,i) = src(IM2,k,j,i) + qomL*src(IDN,k,j,i);
          if (NON_BAROTROPIC_EOS) {
            src(IEN,k,j,i) += (0.5/src(IDN,k,j,i))*(SQR(shboxvar_inner_hydro_(IM2,k,j,i))
                                                    - SQR(src(IM2,k,j,i)));
          } // update energy
          src(IM2,k,j,i) = shboxvar_inner_hydro_(IM2,k,j,i);// update IM2
        }
      }
    }
  }

  if (shbb_.outer == true) {
    int ib = pmb->ie + 1;
    int ii;
    // step 2. -- add shear to the outer periodic boundary values
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib + i;
          // add shear to conservative
          shboxvar_outer_hydro_(IM2,k,j,i) = src(IM2,k,j,ii) - qomL*src(IDN,k,j,ii);
          if (NON_BAROTROPIC_EOS) {
            src(IEN,k,j,ii) += (0.5/src(IDN,k,j,ii))
                               *(SQR(shboxvar_outer_hydro_(IM2,k,j,i))
                                 - SQR(src(IM2,k,j,ii)));
          } // update energy
          src(IM2,k,j,ii) = shboxvar_outer_hydro_(IM2,k,j,i);// update IM2
        }
      }
    }
  }
  return;
}
// --------------------------------------------------------------------------------------
// ! \fn void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers()
//  \brief Send shearing box boundary buffers for hydro variables

void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;

  // KGF: hidden assumption that 2D?
  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  int js = pmb->js;
  int je = pmb->je;

  Real qomL = qshear_*Omega_0_*x1size_;
  int ssize = ssize_*NHYDRO;

  if (shbb_.inner == true) {
    int ib = pmb->is - NGHOST;
    int ii;
    // step 1. -- load shboxvar_hydro_
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_inner_hydro_(IDN,k,j,i) = src(IDN,k,j,ii);
          shboxvar_inner_hydro_(IM1,k,j,i) = src(IM1,k,j,ii);
          shboxvar_inner_hydro_(IM2,k,j,i) = src(IM2,k,j,ii)
                                             + qomL*src(IDN,k,j,ii);
          shboxvar_inner_hydro_(IM3,k,j,i) = src(IM3,k,j,ii);
          if (NON_BAROTROPIC_EOS) {
            shboxvar_inner_hydro_(IEN,k,j,i) = src(IEN,k,j,ii)
                                               + (0.5/src(IDN,k,j,ii))
                                               *(SQR(shboxvar_inner_hydro_(IM2,k,j,i))
                                                 - SQR(src(IM2,k,j,ii)));
          }
        }
      }
    }

    // step 2. -- conservative remaping
    for (int n=0; n<NHYDRO; n++) {
      for (int k=kl; k<=ku; k++) {
        for (int i=0; i<NGHOST; i++) {
          RemapFlux(n, k, js, je+2, i, eps_, shboxvar_inner_hydro_, flx_inner_hydro_);
          for (int j=js; j<=je+1; j++) {
            shboxvar_inner_hydro_(n,k,j,i) -= flx_inner_hydro_(j+1)
                                              - flx_inner_hydro_(j);
          }
        }
      }
    }

    // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, else post MPI_Isend
    for (int n=0; n<4; n++) {
      if (send_inner_rank_[n] != -1) {
        LoadShearing(shboxvar_inner_hydro_, send_innerbuf_hydro_[n], n);
        if (send_inner_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[n]);
          std::memcpy(pbl->pbval->recv_innerbuf_hydro_[n],send_innerbuf_hydro_[n],
                      send_innersize_hydro_[n]*ssize*sizeof(Real));
          pbl->pbval->shbox_inner_hydro_flag_[n] = BoundaryStatus::arrived;
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag=CreateBvalsMPITag(send_inner_lid_[n], n, AthenaTagMPI::shbox_hydro);
          MPI_Isend(send_innerbuf_hydro_[n], send_innersize_hydro_[n]*ssize,
                    MPI_ATHENA_REAL, send_inner_rank_[n], tag, MPI_COMM_WORLD,
                    &rq_innersend_hydro_[n]);
#endif
        }
      }
    }
  } // inner boundaries

  if (shbb_.outer == true) {
    int  ib = pmb->ie + 1;
    qomL = -qomL;
    int ii;
    // step 1. -- load shboxvar_hydro_
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_outer_hydro_(IDN,k,j,i) = src(IDN,k,j,ii);
          shboxvar_outer_hydro_(IM1,k,j,i) = src(IM1,k,j,ii);
          shboxvar_outer_hydro_(IM2,k,j,i) = src(IM2,k,j,ii) + qomL*src(IDN,k,j,ii);
          shboxvar_outer_hydro_(IM3,k,j,i) = src(IM3,k,j,ii);
          if (NON_BAROTROPIC_EOS) {
            shboxvar_outer_hydro_(IEN,k,j,i) = src(IEN,k,j,ii)
                                               + (0.5/src(IDN,k,j,ii))
                                               *(SQR(shboxvar_outer_hydro_(IM2,k,j,i))
                                                 - SQR(src(IM2,k,j,ii)));
          }
        }
      }
    }

    // step 2. -- conservative remaping
    for (int n=0; n<NHYDRO; n++) {
      for (int k=kl; k<=ku; k++) {
        for (int i=0; i<NGHOST; i++) {
          RemapFlux(n, k, js-1, je+1, i, -eps_, shboxvar_outer_hydro_, flx_outer_hydro_);
          for (int j=js-1; j<=je; j++) {
            shboxvar_outer_hydro_(n,k,j,i) -= flx_outer_hydro_(j+1) - flx_outer_hydro_(j);
          }
        }
      }
    }

    // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, post
    // MPI_Isend otherwise
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (send_outer_rank_[n] != -1) {
        LoadShearing(shboxvar_outer_hydro_, send_outerbuf_hydro_[n], n+offset);
        if (send_outer_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[n]);
          std::memcpy(pbl->pbval->recv_outerbuf_hydro_[n],
                      send_outerbuf_hydro_[n],
                      send_outersize_hydro_[n]*ssize*sizeof(Real));
          pbl->pbval->shbox_outer_hydro_flag_[n] = BoundaryStatus::arrived;
        } else { // MPI
#ifdef MPI_PARALLEL
          // bufid for outer(inner): 2(0) and 3(1)
          int tag = CreateBvalsMPITag(send_outer_lid_[n],
                                      n+offset, AthenaTagMPI::shbox_hydro);
          MPI_Isend(send_outerbuf_hydro_[n], send_outersize_hydro_[n]*ssize,
                    MPI_ATHENA_REAL, send_outer_rank_[n], tag, MPI_COMM_WORLD,
                    &rq_outersend_hydro_[n]);
#endif
        }
      }
    }
  } // outer boundaries
  return;
}

// --------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(Real *buf, const int nb)
//  \brief Set hydro shearing box boundary received from a block on the same level

void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  int nxo = pmb->block_size.nx2 - joverlap_;

  sk = pmb->ks; ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) ek += NGHOST, sk -= NGHOST;
  // nb=0-3 for inner boundary; 4-7 for outer boundary.
  switch (nb) {
    case 0:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js - NGHOST; ej = pmb->js + (joverlap_ - 1);
      if (joverlap_ > nx2) sj = pmb->js - nxo;
      break;
    case 1:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js + joverlap_; ej = pmb->je + NGHOST;
      if (joverlap_ < NGHOST) ej = pmb->je + joverlap_;
      break;
    case 2:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js - NGHOST; ej = pmb->js - 1;
      if (joverlap_ > nx2) ej = pmb->js - nxo - 1;
      break;
    case 3:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->je + joverlap_ + 1; ej = pmb->je + NGHOST;
      break;
    case 4:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->je - (joverlap_ - 1); ej = pmb->je + NGHOST;
      if (joverlap_ > nx2) ej = pmb->je + nxo;
      break;
    case 5:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->js - NGHOST; ej = pmb->je - joverlap_;
      if (joverlap_ < NGHOST)   sj = pmb->js - joverlap_;
      break;
    case 6:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->je + 1; ej = pmb->je + NGHOST;
      if (joverlap_ > nx2) sj = pmb->je + nxo + 1;
      break;
    case 7:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->js - NGHOST; ej = pmb->js - joverlap_-1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in CellCenteredBoundaryVariable:SetShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }

  // set [sj:ej] of current meshblock
  int p = 0;
  BufferUtility::UnpackData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}


// --------------------------------------------------------------------------------------
// ! \fn bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers()
//  \brief receive shearing box boundary data for hydro variables

bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  bool flagi = true, flago = true;

  if (shbb_.inner == true) { // check inner boundaries
    for (int n=0; n<4; n++) {
      if (shbox_inner_hydro_flag_[n] == BoundaryStatus::completed) continue;
      if (shbox_inner_hydro_flag_[n] == BoundaryStatus::waiting) {
        if (recv_inner_rank_[n] == Globals::my_rank) {// on the same process
          flagi = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&rq_innerrecv_hydro_[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flagi = false;
            continue;
          }
          shbox_inner_hydro_flag_[n] = BoundaryStatus::arrived;
#endif
        }
      }
      // set dst if boundary arrived
      SetShearingBoxBoundarySameLevel(dst, recv_innerbuf_hydro_[n], n);
      shbox_inner_hydro_flag_[n] = BoundaryStatus::completed; // completed
    } // loop over recv[0] to recv[3]
  } // inner boundary

  if (shbb_.outer == true) { // check outer boundaries
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (shbox_outer_hydro_flag_[n] == BoundaryStatus::completed) continue;
      if (shbox_outer_hydro_flag_[n] == BoundaryStatus::waiting) {
        if (recv_outer_rank_[n] == Globals::my_rank) {// on the same process
          flago = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&rq_outerrecv_hydro_[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flago = false;
            continue;
          }
          shbox_outer_hydro_flag_[n] = BoundaryStatus::arrived;
#endif
        }
      }
      SetShearingBoxBoundarySameLevel(dst, recv_outerbuf_hydro_[n], n+offset);
      shbox_outer_hydro_flag_[n] = BoundaryStatus::completed; // completed
    }
  } // outer boundary
  return (flagi && flago);
}

//--------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::FindShearBlock(const Real time)
//  \brief Calculate the following quantities:
//  send_gid recv_gid send_lid recv_lid send_rank recv_rank,
//  send_size_hydro  recv_size_hydro: for MPI_Irecv
//  eps_,joverlap_: for update the conservative

void CellCenteredBoundaryVariable::FindShearBlock(const Real time) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  Mesh *pmesh = pmb->pmy_mesh;

  int js = pmb->js; int je = pmb->je;

  int level = pmb->loc.level - pmesh->root_level;
  std::int64_t nrbx2 = pmesh->nrbx2*(1L<<level);
  int nx2   = pmb->block_size.nx2; // # of cells per meshblock
  int nx3   = pmb->block_size.nx3; // # of cells per meshblock
  // KGF: for symmetry reasons, how can ncells3 but not ncells2 be used in this fn?
  // int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
  int ncells3 = pmb->block_size.nx3;
  if (pmesh->mesh_size.nx3 > 1) ncells3 += 2*NGHOST;

  Real qomL = qshear_*Omega_0_*x1size_;
  Real yshear = qomL*time;
  Real deltay = fmod(yshear,x2size_);
  int joffset = static_cast<int>(deltay/pco->dx2v(js)); // assumes uniform grid in azimuth
  int Ngrids  = static_cast<int>(joffset/nx2);
  joverlap_   = joffset - Ngrids*nx2;
  eps_ = (std::fmod(deltay, pco->dx2v(js)))/pco->dx2v(js);

  if (shbb_.inner == true) { // if inner block
    for (int n=0; n<4; n++) {
      send_inner_gid_[n]  = -1;
      send_inner_rank_[n] = -1;
      send_inner_lid_[n]  = -1;
      recv_inner_gid_[n]  = -1;
      recv_inner_rank_[n] = -1;
      recv_inner_lid_[n]  = -1;
      send_innersize_hydro_[n] = 0;
      recv_innersize_hydro_[n] = 0;
      shbox_inner_hydro_flag_[n] = BoundaryStatus::completed;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[n] = 0;
        recv_innersize_field_[n] = 0;
        shbox_inner_field_flag_[n] = BoundaryStatus::completed;
        send_innersize_emf_[n] = 0;
        recv_innersize_emf_[n] = 0;
        shbox_inner_emf_flag_[n] = BoundaryStatus::completed;
      }
    }
    int jblock = 0;
    for (int j=0; j<nrbx2; j++) {
      // index of current meshblock on the shearingboundary block list
      if (shbb_.igidlist[j] == pmb->gid)  jblock = j;
    }
    // send [js:je-joverlap] of the meshblock to other
    // attach [je-joverlap+1:MIN(je-joverlap+(NGHOST),je-js+1)]
    // to its right end.
    std::int64_t jtmp = jblock + Ngrids;
    if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    send_inner_gid_[1]  = shbb_.igidlist[jtmp];
    send_inner_rank_[1] = shbb_.irnklist[jtmp];
    send_inner_lid_[1]  = shbb_.ilidlist[jtmp];
    send_innersize_hydro_[1] = std::min(je-js-joverlap_+1+NGHOST, je-js+1);
    // recv [js+joverlap:je] from other
    // attach [je+1:MIN(je+NGHOST,je+joverlap)] to its right end.
    jtmp = jblock - Ngrids;
    if (jtmp < 0) jtmp += nrbx2;
    recv_inner_gid_[1]  = shbb_.igidlist[jtmp];
    recv_inner_rank_[1] = shbb_.irnklist[jtmp];
    recv_inner_lid_[1]  = shbb_.ilidlist[jtmp];
    recv_innersize_hydro_[1] = send_innersize_hydro_[1];
    shbox_inner_hydro_flag_[1] = BoundaryStatus::waiting;
    if (MAGNETIC_FIELDS_ENABLED) {
      send_innersize_field_[1] = send_innersize_hydro_[1]
                                 *NGHOST*(NFIELD*ncells3+1)
                                 +NGHOST*ncells3;
      recv_innersize_field_[1] = send_innersize_field_[1];
      shbox_inner_field_flag_[1] = BoundaryStatus::waiting;
      send_innersize_emf_[1] = send_innersize_hydro_[1]*(2*nx3 + 1) + nx3;
      recv_innersize_emf_[1] = send_innersize_emf_[1];
      shbox_inner_emf_flag_[1] = BoundaryStatus::waiting;
    }


    // if there is overlap to next blocks
    if (joverlap_ != 0) {
      // send to the right
      jtmp = jblock + (Ngrids + 1);
      if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[0]  = shbb_.igidlist[jtmp];
      send_inner_rank_[0] = shbb_.irnklist[jtmp];
      send_inner_lid_[0]  = shbb_.ilidlist[jtmp];
      send_innersize_hydro_[0] = std::min(joverlap_+NGHOST, je-js+1);
      // receive from its left
      jtmp = jblock - (Ngrids + 1);
      if (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[0]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[0] = shbb_.irnklist[jtmp];
      recv_inner_lid_[0]  = shbb_.ilidlist[jtmp];
      recv_innersize_hydro_[0] = send_innersize_hydro_[0];
      shbox_inner_hydro_flag_[0] = BoundaryStatus::waiting;// switch on if overlap
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[0] = send_innersize_hydro_[0]
                                   *NGHOST*(NFIELD*ncells3+1)
                                   +NGHOST*ncells3;
        recv_innersize_field_[0] = send_innersize_field_[0];
        shbox_inner_field_flag_[0] = BoundaryStatus::waiting;
        send_innersize_emf_[0] = send_innersize_hydro_[0]*(2*nx3+1)+nx3;
        recv_innersize_emf_[0] = send_innersize_emf_[0];
        shbox_inner_emf_flag_[0] = BoundaryStatus::waiting;
      }
      // deal the left boundary cells with send[2]
      if (joverlap_ > (nx2 - NGHOST)) {
        // send to Right
        jtmp = jblock + (Ngrids + 2);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        send_inner_gid_[2]  = shbb_.igidlist[jtmp];
        send_inner_rank_[2] = shbb_.irnklist[jtmp];
        send_inner_lid_[2]  = shbb_.ilidlist[jtmp];
        send_innersize_hydro_[2] = joverlap_-(nx2-NGHOST);
        // recv from Left
        jtmp = jblock - (Ngrids+2);
        while (jtmp < 0) jtmp += nrbx2;
        recv_inner_gid_[2]  = shbb_.igidlist[jtmp];
        recv_inner_rank_[2] = shbb_.irnklist[jtmp];
        recv_inner_lid_[2]  = shbb_.ilidlist[jtmp];
        recv_innersize_hydro_[2] = send_innersize_hydro_[2];
        shbox_inner_hydro_flag_[2] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_field_[2] = send_innersize_hydro_[2]*NGHOST*(NFIELD*ncells3+1);
          recv_innersize_field_[2] = send_innersize_field_[2];
          shbox_inner_field_flag_[2] = BoundaryStatus::waiting;
          send_innersize_emf_[2] = send_innersize_hydro_[2]*(2*nx3+1);
          recv_innersize_emf_[2] = send_innersize_emf_[2];
          shbox_inner_emf_flag_[2] = BoundaryStatus::waiting;
        }
      }
      // deal with the right boundary cells with send[3]
      if (joverlap_ < NGHOST) {
        jtmp = jblock + (Ngrids - 1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        send_inner_gid_[3]  = shbb_.igidlist[jtmp];
        send_inner_rank_[3] = shbb_.irnklist[jtmp];
        send_inner_lid_[3]  = shbb_.ilidlist[jtmp];
        send_innersize_hydro_[3] = NGHOST-joverlap_;
        jtmp = jblock - (Ngrids-1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        recv_inner_gid_[3]  = shbb_.igidlist[jtmp];
        recv_inner_rank_[3] = shbb_.irnklist[jtmp];
        recv_inner_lid_[3]  = shbb_.ilidlist[jtmp];
        recv_innersize_hydro_[3] = send_innersize_hydro_[3];
        shbox_inner_hydro_flag_[3] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_field_[3] = send_innersize_hydro_[3]*NGHOST
                                     *(NFIELD*ncells3+1);
          recv_innersize_field_[3] = send_innersize_field_[3];
          shbox_inner_field_flag_[3] = BoundaryStatus::waiting;
          send_innersize_emf_[3] = send_innersize_hydro_[3]*(2*nx3+1);
          recv_innersize_emf_[3] = send_innersize_emf_[3];
          shbox_inner_emf_flag_[3] = BoundaryStatus::waiting;
        }
      }
    } else { // joverlap_ == 0
      // send [je-(NGHOST-1):je] to Right
      jtmp = jblock + (Ngrids+1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[2]  = shbb_.igidlist[jtmp];
      send_inner_rank_[2] = shbb_.irnklist[jtmp];
      send_inner_lid_[2]  = shbb_.ilidlist[jtmp];
      send_innersize_hydro_[2] = NGHOST;
      // recv [js-NGHOST:js-1] from Left
      jtmp = jblock - (Ngrids+1);
      while (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[2]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[2] = shbb_.irnklist[jtmp];
      recv_inner_lid_[2]  = shbb_.ilidlist[jtmp];
      recv_innersize_hydro_[2] = send_innersize_hydro_[2];
      shbox_inner_hydro_flag_[2] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[2] = send_innersize_hydro_[2]*NGHOST
                                   *(NFIELD*ncells3+1);
        recv_innersize_field_[2] = send_innersize_field_[2];
        shbox_inner_field_flag_[2] = BoundaryStatus::waiting;
        send_innersize_emf_[2] = send_innersize_hydro_[2]*(2*nx3+1);
        recv_innersize_emf_[2] = send_innersize_emf_[2];
        shbox_inner_emf_flag_[2] = BoundaryStatus::waiting;
      }

      // send [js:js+(NGHOST-1)] to Left
      jtmp = jblock + (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      send_inner_gid_[3]  = shbb_.igidlist[jtmp];
      send_inner_rank_[3] = shbb_.irnklist[jtmp];
      send_inner_lid_[3]  = shbb_.ilidlist[jtmp];
      send_innersize_hydro_[3] = NGHOST;
      // recv [je+1:je+(NGHOST-1)] from Right
      jtmp = jblock - (Ngrids-1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[3]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[3] = shbb_.irnklist[jtmp];
      recv_inner_lid_[3]  = shbb_.ilidlist[jtmp];
      recv_innersize_hydro_[3] = send_innersize_hydro_[3];
      shbox_inner_hydro_flag_[3] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[3] = send_innersize_hydro_[3]*NGHOST
                                   *(NFIELD*ncells3+1);
        recv_innersize_field_[3] = send_innersize_field_[3];
        shbox_inner_field_flag_[3] = BoundaryStatus::waiting;
        send_innersize_emf_[3] = send_innersize_hydro_[3]*(2*nx3+1);
        recv_innersize_emf_[3] = send_innersize_emf_[3];
        shbox_inner_emf_flag_[3] = BoundaryStatus::waiting;
      }
    }
  } // inner bc

  if (shbb_.outer == true) { // if outer block
    for (int n=0; n<4; n++) {
      send_outer_gid_[n]  = -1;
      send_outer_rank_[n] = -1;
      send_outer_lid_[n]  = -1;
      recv_outer_gid_[n]  = -1;
      recv_outer_rank_[n] = -1;
      recv_outer_lid_[n]  = -1;
      send_outersize_hydro_[n] = 0;
      recv_outersize_hydro_[n] = 0;
      shbox_outer_hydro_flag_[n] = BoundaryStatus::completed;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[n] = 0;
        recv_outersize_field_[n] = 0;
        shbox_outer_field_flag_[n] = BoundaryStatus::completed;
        send_outersize_emf_[n] = 0;
        recv_outersize_emf_[n] = 0;
        shbox_outer_emf_flag_[n] = BoundaryStatus::completed;
      }
    }
    int jblock = 0;
    for (int j=0; j<nrbx2; j++) {
      // index of current meshblock on the shearingboundary block list
      if (shbb_.ogidlist[j] == pmb->gid) jblock = j;
    }
    // recv [js-NGHOST:je-joverlap] of the meshblock from other
    std::int64_t jtmp = jblock + Ngrids;
    if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    recv_outer_gid_[1]  = shbb_.ogidlist[jtmp];
    recv_outer_rank_[1] = shbb_.ornklist[jtmp];
    recv_outer_lid_[1]  = shbb_.olidlist[jtmp];
    recv_outersize_hydro_[1] = std::min(je-js-joverlap_+1+NGHOST,je-js+1);
    // send [js+joverlap-NGHOST:je] of the meshblock to other
    jtmp = jblock - Ngrids;
    if (jtmp < 0) jtmp += nrbx2;
    send_outer_gid_[1]  = shbb_.ogidlist[jtmp];
    send_outer_rank_[1] = shbb_.ornklist[jtmp];
    send_outer_lid_[1]  = shbb_.olidlist[jtmp];
    send_outersize_hydro_[1] = recv_outersize_hydro_[1];
    shbox_outer_hydro_flag_[1] = BoundaryStatus::waiting;
    if (MAGNETIC_FIELDS_ENABLED) {
      send_outersize_field_[1] = send_outersize_hydro_[1]*NGHOST*(NFIELD*ncells3+1) +
                                 NGHOST*ncells3;
      recv_outersize_field_[1] = send_outersize_field_[1];
      shbox_outer_field_flag_[1] = BoundaryStatus::waiting;
      send_outersize_emf_[1] = send_outersize_hydro_[1]*(2*nx3+1)+nx3;
      recv_outersize_emf_[1] = send_outersize_emf_[1];
      shbox_outer_emf_flag_[1] = BoundaryStatus::waiting;
    }

    // if there is overlap to next blocks
    if (joverlap_ != 0) {
      // recv from right
      jtmp = jblock + (Ngrids + 1);
      if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[0]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[0] = shbb_.ornklist[jtmp];
      recv_outer_lid_[0]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[0] = std::min(joverlap_+NGHOST,je-js+1);
      // send to left
      jtmp = jblock - (Ngrids + 1);
      if (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[0]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[0] = shbb_.ornklist[jtmp];
      send_outer_lid_[0]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[0] = recv_outersize_hydro_[0];
      shbox_outer_hydro_flag_[0] = BoundaryStatus::waiting; // switch on if overlap
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[0] = send_outersize_hydro_[0]
                                   *NGHOST*(NFIELD*ncells3+1)
                                   +NGHOST*ncells3;
        recv_outersize_field_[0] = send_outersize_field_[0];
        shbox_outer_field_flag_[0] = BoundaryStatus::waiting;
        send_outersize_emf_[0] = send_outersize_hydro_[0]*(2*nx3+1)+nx3;
        recv_outersize_emf_[0] = send_outersize_emf_[0];
        shbox_outer_emf_flag_[0] = BoundaryStatus::waiting;
      }
      // deal the left boundary cells with send[2]
      if (joverlap_ > (nx2 - NGHOST)) {
        // recv from left
        jtmp = jblock + (Ngrids + 2);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        recv_outer_gid_[2]  = shbb_.ogidlist[jtmp];
        recv_outer_rank_[2] = shbb_.ornklist[jtmp];
        recv_outer_lid_[2]  = shbb_.olidlist[jtmp];
        recv_outersize_hydro_[2] = joverlap_-(nx2-NGHOST);
        // send to right
        jtmp = jblock - (Ngrids+2);
        while (jtmp < 0) jtmp += nrbx2;
        send_outer_gid_[2]  = shbb_.ogidlist[jtmp];
        send_outer_rank_[2] = shbb_.ornklist[jtmp];
        send_outer_lid_[2]  = shbb_.olidlist[jtmp];
        send_outersize_hydro_[2] = recv_outersize_hydro_[2];
        shbox_outer_hydro_flag_[2] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_field_[2] = send_outersize_hydro_[2]*NGHOST*(NFIELD*ncells3+1);
          recv_outersize_field_[2] = send_outersize_field_[2];
          shbox_outer_field_flag_[2] = BoundaryStatus::waiting;
          send_outersize_emf_[2] = send_outersize_hydro_[2]*(2*nx3+1);
          recv_outersize_emf_[2] = send_outersize_emf_[2];
          shbox_outer_emf_flag_[2] = BoundaryStatus::waiting;
        }
      }
      // deal the right boundary cells with send[3]
      if (joverlap_ < NGHOST) {
        jtmp = jblock + (Ngrids - 1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        recv_outer_gid_[3]  = shbb_.ogidlist[jtmp];
        recv_outer_rank_[3] = shbb_.ornklist[jtmp];
        recv_outer_lid_[3]  = shbb_.olidlist[jtmp];
        recv_outersize_hydro_[3] = NGHOST-joverlap_;
        jtmp = jblock - (Ngrids-1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        send_outer_gid_[3]  = shbb_.ogidlist[jtmp];
        send_outer_rank_[3] = shbb_.ornklist[jtmp];
        send_outer_lid_[3]  = shbb_.olidlist[jtmp];
        send_outersize_hydro_[3] = recv_outersize_hydro_[3];
        shbox_outer_hydro_flag_[3] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_field_[3] = send_outersize_hydro_[3]*NGHOST
                                     *(NFIELD*ncells3+1);
          recv_outersize_field_[3] = send_outersize_field_[3];
          shbox_outer_field_flag_[3] = BoundaryStatus::waiting;
          send_outersize_emf_[3] = send_outersize_hydro_[3]*(2*nx3+1);
          recv_outersize_emf_[3] = send_outersize_emf_[3];
          shbox_outer_emf_flag_[3] = BoundaryStatus::waiting;
        }
      }
    } else { // joverlap_ == 0
      // recv [je+1:je+NGHOST] from Left
      jtmp = jblock + (Ngrids + 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[2]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[2] = shbb_.ornklist[jtmp];
      recv_outer_lid_[2]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[2] = NGHOST;
      // send [js:js+NGHOST-1] to Right
      jtmp = jblock - (Ngrids+1);
      while (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[2]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[2] = shbb_.ornklist[jtmp];
      send_outer_lid_[2]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[2] = NGHOST;
      shbox_outer_hydro_flag_[2] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[2] = send_outersize_hydro_[2]
                                   *NGHOST*(NFIELD*ncells3+1);
        recv_outersize_field_[2] = send_outersize_field_[2];
        shbox_outer_field_flag_[2] = BoundaryStatus::waiting;
        send_outersize_emf_[2] = send_outersize_hydro_[2]*(2*nx3+1);
        recv_outersize_emf_[2] = send_outersize_emf_[2];
        shbox_outer_emf_flag_[2] = BoundaryStatus::waiting;
      }

      // recv [js-NGHOST:js-1] from Left
      jtmp = jblock + (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      recv_outer_gid_[3]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[3] = shbb_.ornklist[jtmp];
      recv_outer_lid_[3]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[3] = NGHOST;
      // send [je-(NGHOST-1):je] to Right
      jtmp = jblock - (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[3]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[3] = shbb_.ornklist[jtmp];
      send_outer_lid_[3]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[3] = NGHOST;
      shbox_outer_hydro_flag_[3] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[3] = send_outersize_hydro_[3]
                                   *NGHOST*(NFIELD*ncells3+1);
        recv_outersize_field_[3] = send_outersize_field_[3];
        shbox_outer_field_flag_[3] = BoundaryStatus::waiting;
        send_outersize_emf_[3] = send_outersize_hydro_[3]*(2*nx3+1);
        recv_outersize_emf_[3] = send_outersize_emf_[3];
        shbox_outer_emf_flag_[3] = BoundaryStatus::waiting;
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

//  \brief compute the flux along j indices for remapping adopted from 2nd
//  order RemapFlux of Athena 4.2 (C-version)

void CellCenteredBoundaryVariable::RemapFlux(const int n, const int k, const int jinner,
                               const int jouter, const int i, const Real eps,
                               const AthenaArray<Real> &var, AthenaArray<Real> &flux) {
  int j, jl, ju;
  Real dUc, dUl, dUr, dUm, lim_slope;

  // jinner, jouter are index range over which flux must be returned.  Set loop
  // limits depending on direction of upwind differences
  if (eps > 0.0) { // eps always > 0 for inner i boundary
    jl = jinner - 1;
    ju = jouter - 1;
  } else {         // eps always < 0 for outer i boundary
    jl = jinner;
    ju = jouter;
  }

  // TODO(felker): do not reimplement PLM here; use plm.cpp.
  // TODO(felker): relax assumption that 2nd order reconstruction must be used
  for (j=jl; j<=ju; j++) {
    dUc = var(n,k,j+1,i) - var(n,k,j-1,i);
    dUl = var(n,k,j,  i) - var(n,k,j-1,i);
    dUr = var(n,k,j+1,i) - var(n,k,j,  i);

    dUm = 0.0;
    if (dUl*dUr > 0.0) {
      lim_slope = std::min(std::fabs(dUl), std::fabs(dUr));
      dUm = SIGN(dUc)*std::min(0.5*std::fabs(dUc), 2.0*lim_slope);
    }

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      flux(j+1) = eps*(var(n,k,j,i) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      flux(j  ) = eps*(var(n,k,j,i) - 0.5*(1.0 + eps)*dUm);
    }
  }
  return;
}
