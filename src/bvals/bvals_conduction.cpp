//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_conduction.cpp
//  \brief functions that apply BCs for Landau fluid conduction operator

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "bvals_conduction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class FFTBlock;
class FFTDriver;

//----------------------------------------------------------------------------------------
//! \fn ConductionBoundaryValues::ConductionBoundaryValues(MeshbBlock *pmb,
//                                                   enum BoundaryFlag *input_bcs)
//  \brief Constructor of the ConductionBoundaryValues class
//    This sets the boundaries after the fft calls on d_prl(T) etc.
//    Tprl, Tprp, and B all have the same boundary conditions (periodic at this stage)

ConductionBoundaryValues::ConductionBoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs)
: BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs) {
  pmy_block_=pmb;
  for (int i=0; i<6; i++) {
    if (block_bcs[i] == PERIODIC_BNDRY || block_bcs[i]==BLOCK_BNDRY)
      ConductionBoundaryFunction_[i]=NULL;
    // else
  }
  
  SearchAndSetNeighbors(pmy_mesh_->tree, pmy_mesh_->ranklist, pmy_mesh_->nslist);
  
  num_fields_ = 3;
  InitBoundaryData(bd_conduction_);
}


//----------------------------------------------------------------------------------------
//! \fn ConductionBoundaryValues::~ConductionBoundaryValues()
//  \brief Destructor of the ConductionBoundaryValues class

ConductionBoundaryValues::~ConductionBoundaryValues() {
  DestroyBoundaryData(bd_conduction_);
}

//----------------------------------------------------------------------------------------
//! \fn void ConductionBoundaryValues::InitBoundaryData(ConductionBoundaryData &bd)
//  \brief Initialize ConductionBoundaryData structure
void ConductionBoundaryValues::InitBoundaryData(ConductionBoundaryData &bd) {
  MeshBlock *pmb=pmy_block_;
  int size;
  bd.nbmax=maxneighbor_;
  for (int n=0;n<bd.nbmax;n++) {
    // Clear flags and requests
    bd.flag[n]=BNDRY_WAITING;
    bd.sflag[n]=BNDRY_WAITING;
    bd.send[n]=NULL;
    bd.recv[n]=NULL;
#ifdef MPI_PARALLEL
    bd.req_send[n]=MPI_REQUEST_NULL;
    bd.req_recv[n]=MPI_REQUEST_NULL;
#endif
    
    // Allocate buffers
    // calculate the buffer size
    size=((BoundaryValues::ni[n].ox1==0)?pmb->block_size.nx1:NGHOST)
        *((BoundaryValues::ni[n].ox2==0)?pmb->block_size.nx2:NGHOST)
        *((BoundaryValues::ni[n].ox3==0)?pmb->block_size.nx3:NGHOST);
    
    size *= num_fields_;
    
    bd.send[n] = new Real[size];
    bd.recv[n] = new Real[size];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ConductionBoundaryValues::DestroyBoundaryData(ConductionBoundaryData &bd)
//  \brief Destroy ConductionBoundaryData structure
void ConductionBoundaryValues::DestroyBoundaryData(ConductionBoundaryData &bd) {
  for (int n=0;n<bd.nbmax;n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
#ifdef MPI_PARALLEL
    if (bd.req_send[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_send[n]);
    if (bd.req_recv[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_recv[n]);
#endif
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ConductionBoundaryValues::ApplyPhysicalBoundaries(void)
//  \brief Apply physical boundary conditions to the parallel derivatives of T and B

void ConductionBoundaryValues::ApplyPhysicalBoundaries(void) {
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> &dst=pmb->phydro->phdif->dprl_cond;
  int bis=pmb->is-NGHOST, bie=pmb->ie+NGHOST, bjs=pmb->js, bje=pmb->je,
  bks=pmb->ks, bke=pmb->ke;
  Real time=pmy_mesh_->time;
  Real dt=pmy_mesh_->dt;
  if (ConductionBoundaryFunction_[INNER_X2]==NULL
      && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
  if (ConductionBoundaryFunction_[OUTER_X2]==NULL
      && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
  if (ConductionBoundaryFunction_[INNER_X3]==NULL
      && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
  if (ConductionBoundaryFunction_[OUTER_X3]==NULL
      && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;
  
  // Apply boundary function on inner-x1
  if (ConductionBoundaryFunction_[INNER_X1] != NULL)
    ConductionBoundaryFunction_[INNER_X1](pmb, pco, dst, time, dt,
                                       pmb->is, pmb->ie, bjs, bje, bks, bke);
  // Apply boundary function on outer-x1
  if (ConductionBoundaryFunction_[OUTER_X1] != NULL)
    ConductionBoundaryFunction_[OUTER_X1](pmb, pco, dst, time, dt,
                                       pmb->is, pmb->ie, bjs, bje, bks, bke);
  
  if (pmb->block_size.nx2>1) { // 2D or 3D
    
    // Apply boundary function on inner-x2
    if (ConductionBoundaryFunction_[INNER_X2] != NULL)
      ConductionBoundaryFunction_[INNER_X2](pmb, pco, dst, time, dt,
                                         bis, bie, pmb->js, pmb->je, bks, bke);
    // Apply boundary function on outer-x2
    if (ConductionBoundaryFunction_[OUTER_X2] != NULL)
      ConductionBoundaryFunction_[OUTER_X2](pmb, pco, dst, time, dt,
                                         bis, bie, pmb->js, pmb->je, bks, bke);
  }
  
  if (pmb->block_size.nx3>1) { // 3D
    
    // Apply boundary function on inner-x3
    if (ConductionBoundaryFunction_[INNER_X3] != NULL)
      ConductionBoundaryFunction_[INNER_X3](pmb, pco, dst, time, dt,
                                         bis, bie, bjs, bje, pmb->ks, pmb->ke);
    // Apply boundary function on outer-x3
    if (ConductionBoundaryFunction_[OUTER_X3] != NULL)
      ConductionBoundaryFunction_[OUTER_X3](pmb, pco, dst, time, dt,
                                         bis, bie, bjs, bje, pmb->ks, pmb->ke);
  }
  
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ConductionBoundaryValues::StartReceivingConduction(void)
//  \brief initiate MPI_Irecv for Conduction

void ConductionBoundaryValues::StartReceivingConduction(void) {
  MeshBlock *pmb=pmy_block_;
  int tag;
  ConductionBoundaryData *pbd;
  
  pbd=&bd_conduction_;
  
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank) {
      int size;
      size=((nb.ox1==0)?pmb->block_size.nx1:NGHOST)
        *((nb.ox2==0)?pmb->block_size.nx2:NGHOST)
        *((nb.ox3==0)?pmb->block_size.nx3:NGHOST);
      size *= num_fields_;
      tag=CreateBvalsMPITag(pmb->lid, TAG_CONDUCTION, nb.bufid);
      MPI_Irecv(pbd->recv[nb.bufid], size, MPI_ATHENA_REAL,
                nb.rank, tag, MPI_COMM_WORLD, &(pbd->req_recv[nb.bufid]));
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ConductionBoundaryValues::ClearBoundaryConduction(void)
//  \brief clean up the boundary flags after each loop

void ConductionBoundaryValues::ClearBoundaryConduction(void) {
  ConductionBoundaryData *pbd;
  
  pbd=&bd_conduction_;
  
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    pbd->flag[nb.bufid] = BNDRY_WAITING;
    pbd->sflag[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank)
      MPI_Wait(&(pbd->req_send[nb.bufid]),MPI_STATUS_IGNORE); // Wait for Isend
#endif
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int ConductionBoundaryValues::LoadConductionBoundaryBufferSameLevel(AthenaArray<Real>
//                                 &src, int ns, int ne, Real *buf, const NeighborBlock& nb)
//  \brief Set Conduction boundary buffers for sending to a block on the same level

int ConductionBoundaryValues::LoadConductionBoundaryBufferSameLevel(AthenaArray<Real> &src,
                                                                    int ns, int ne,
                                                      Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  
  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;
  int p=0;
  BufferUtility::Pack4DData(src, buf, ns, ne, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void ConductionBoundaryValues::SendConductionBoundaryBuffers(AthenaArray<Real> &src)
//  \brief Send boundary buffers

bool ConductionBoundaryValues::SendConductionBoundaryBuffers(AthenaArray<Real> &src) {
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  ConductionBoundaryData *pbd, *ptarget;
  bool bflag=true;
  int tag;
  int ns=0, ne=num_fields_-1;
  
  
  pbd=&bd_conduction_;
  
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (pbd->sflag[nb.bufid]==BNDRY_COMPLETED) continue;
    int ssize;
    
    if (nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      ptarget=&(pbl->pcondbval->bd_conduction_);
      if (ptarget->flag[nb.targetid] != BNDRY_WAITING) {
        bflag=false;
        continue;
      }
    }
    if (nb.level==mylevel)
      ssize=LoadConductionBoundaryBufferSameLevel(src, ns, ne, pbd->send[nb.bufid],nb);
    
    if (nb.rank == Globals::my_rank) { // on the same process
      std::memcpy(ptarget->recv[nb.targetid], pbd->send[nb.bufid], ssize*sizeof(Real));
      ptarget->flag[nb.targetid]=BNDRY_ARRIVED;
#ifdef MPI_PARALLEL
    } else { // MPI
      tag=CreateBvalsMPITag(nb.lid, TAG_CONDUCTION, nb.targetid);
      MPI_Isend(pbd->send[nb.bufid], ssize, MPI_ATHENA_REAL, nb.rank, tag,
                MPI_COMM_WORLD, &(pbd->req_send[nb.bufid]));
    }
#else
  }
#endif
  pbd->sflag[nb.bufid] = BNDRY_COMPLETED;
}

return bflag;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetConductionBoundarySameLevel(AthenaArray<Real> &dst,
//                                           Real *buf, const NeighborBlock& nb)
//  \brief Set Conduction boundary received from a block on the same level

void ConductionBoundaryValues::SetConductionBoundarySameLevel(AthenaArray<Real> &dst,
                                                        int ns, int ne, Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  
  if (nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if (nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if (nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if (nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if (nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  
  int p=0;
  // Conduction only works with Cartesian coordinate
  BufferUtility::Unpack4DData(buf, dst, ns, ne, si, ei, sj, ej, sk, ek, p);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool ConductionBoundaryValues::ReceiveConductionBoundaryBuffers(AthenaArray<Real> &dst)
//  \brief receive the boundary data

bool ConductionBoundaryValues::ReceiveConductionBoundaryBuffers(AthenaArray<Real> &dst) {
  MeshBlock *pmb=pmy_block_;
  bool flag=true;
  ConductionBoundaryData *pbd;
  int ns=0, ne=num_fields_-1;
  
  pbd=&bd_conduction_;
  
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (pbd->flag[nb.bufid]==BNDRY_COMPLETED) continue;
    if (pbd->flag[nb.bufid]==BNDRY_WAITING) {
      if (nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
#ifdef MPI_PARALLEL
      } else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&(pbd->req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
        if (static_cast<bool>(test)==false) {
          flag=false;
          continue;
        }
        pbd->flag[nb.bufid] = BNDRY_ARRIVED;
      }
#else
    }
#endif
  }
  if (nb.level==pmb->loc.level)
    SetConductionBoundarySameLevel(dst, ns, ne, pbd->recv[nb.bufid], nb);
  //    else
  // error message
  pbd->flag[nb.bufid] = BNDRY_COMPLETED; // completed
}

return flag;
}
