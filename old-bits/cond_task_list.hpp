#ifndef TASK_LIST_COND_TASK_LIST_HPP_
#define TASK_LIST_COND_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file cond_task_list.hpp
//    \brief Handles boundary send/receives for application of LF heat conduction operator

#include <stdint.h>

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//----------------------------------------------------------------------------------------
//! \class ConductionSolverTaskList
//  \brief data and function definitions for ConductionSolverTaskList derived class

class ConductionSolverTaskList : public TaskList {
public:
  ConductionSolverTaskList(ParameterInput *pin, Mesh *pm);
  ~ConductionSolverTaskList() {}
  
  void AddConductionSolverTask(uint64_t id, uint64_t dep);
  
  // functions
  enum TaskStatus StartConductionReceive(MeshBlock *pmb, int stage);
  enum TaskStatus ClearConductionBoundary(MeshBlock *pmb, int stage);
  enum TaskStatus SendConductionBoundary(MeshBlock *pmb, int stage);
  enum TaskStatus ReceiveConductionBoundary(MeshBlock *pmb, int stage);
  enum TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);
};


//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.

namespace ConductionSolverTaskNames {
  const uint64_t NONE=0;
  const uint64_t START_COND_RECV=1LL<<0;
  const uint64_t CLEAR_COND=1LL<<1;
  
  const uint64_t SEND_COND_BND=1LL<<2;
  const uint64_t RECV_COND_BND=1LL<<3;
  const uint64_t COND_PHYS_BND=1LL<<4;
};

#endif // TASK_LIST_COND_TASK_LIST_HPP_
