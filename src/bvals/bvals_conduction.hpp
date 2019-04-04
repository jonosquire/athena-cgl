#ifndef BVALS_BVALS_CONDUCTION_HPP_
#define BVALS_BVALS_CONDUCTION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_conduction.hpp
//  \brief defines ConductionBoundaryValues class

// C++ headers
#include <string>   // string

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "bvals.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class FFTBlock;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;

//! \struct ConductionBoundaryData
//  \brief structure storing heat conduction boundary information
typedef struct ConductionBoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56], sflag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} ConductionBoundaryData;

//----------------------------------------------------------------------------------------
//! \class ConductionBoundaryValues
//  \brief BVals data and functions

class ConductionBoundaryValues : public BoundaryBase {
public:
  ConductionBoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs);
  ~ConductionBoundaryValues();
  
  void InitBoundaryData(ConductionBoundaryData &bd);
  void DestroyBoundaryData(ConductionBoundaryData &bd);
  
  void ApplyPhysicalBoundaries(void);
  void StartReceivingConduction(void);
  void ClearBoundaryConduction(void);
  int LoadConductionBoundaryBufferSameLevel(AthenaArray<Real> &src, int ns, int ne,
                                            Real *buf, const NeighborBlock& nb);
  bool SendConductionBoundaryBuffers(AthenaArray<Real> &src);
  void SetConductionBoundarySameLevel(AthenaArray<Real> &dst, int ns, int ne,
                                      Real *buf, const NeighborBlock& nb);
  bool ReceiveConductionBoundaryBuffers(AthenaArray<Real> &dst);
  
private:
  MeshBlock *pmy_block_;
  ConductionBoundaryFunc_t ConductionBoundaryFunction_[6];
  ConductionBoundaryData bd_conduction_;
  int num_fields_; // Number of fields across which to apply BCs
};

#endif // BVALS_BVALS_CONDUCTION_HPP_
