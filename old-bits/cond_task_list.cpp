//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cond_task_list.cpp
//  \brief

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "task_list.hpp"
#include "cond_task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../bvals/bvals_conduction.hpp"

//----------------------------------------------------------------------------------------
//  ConductionSolverTaskList constructor

ConductionSolverTaskList::ConductionSolverTaskList(ParameterInput *pin, Mesh *pm)
: TaskList(pm) {
  
  // Now assemble list of tasks for each stage of time integrator
  {using namespace ConductionSolverTaskNames; // NOLINT (build/namespace)
    AddConductionSolverTask(START_COND_RECV,NONE);
    
    // Send parallel gradient boundaries
    AddConductionSolverTask(SEND_COND_BND,START_COND_RECV);
    AddConductionSolverTask(RECV_COND_BND,START_COND_RECV);
    AddConductionSolverTask(COND_PHYS_BND,SEND_COND_BND|RECV_COND_BND);
    AddConductionSolverTask(CLEAR_COND, COND_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void ConductionSolverTaskList::AddConductionSolverTask(uint64_t id, uint64_t dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void ConductionSolverTaskList::AddConductionSolverTask(uint64_t id, uint64_t dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;
  
  using namespace ConductionSolverTaskNames; // NOLINT (build/namespace)
  switch((id)) {
    case (START_COND_RECV):
      task_list_[ntasks].TaskFunc=
      static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
      (&ConductionSolverTaskList::StartConductionReceive);
      break;
    case (CLEAR_COND):
      task_list_[ntasks].TaskFunc=
      static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
      (&ConductionSolverTaskList::ClearConductionBoundary);
      break;
    case (SEND_COND_BND):
      task_list_[ntasks].TaskFunc=
      static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
      (&ConductionSolverTaskList::SendConductionBoundary);
      break;
    case (RECV_COND_BND):
      task_list_[ntasks].TaskFunc=
      static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
      (&ConductionSolverTaskList::ReceiveConductionBoundary);
      break;
    case (COND_PHYS_BND):
      task_list_[ntasks].TaskFunc=
      static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
      (&ConductionSolverTaskList::PhysicalBoundary);
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddConductionSolverTask" << std::endl
      << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus ConductionSolverTaskList::StartConductionReceive(MeshBlock *pmb, int stage) {
  pmb->pcondbval->StartReceivingConduction();
  std::cout << "At StartConductionReceive \n";
  return TASK_SUCCESS;
}

enum TaskStatus ConductionSolverTaskList::ClearConductionBoundary(MeshBlock *pmb, int stage) {
  pmb->pcondbval->ClearBoundaryConduction();
  std::cout << "At ClearConductionBoundary \n";
  return TASK_SUCCESS;
}

enum TaskStatus ConductionSolverTaskList::SendConductionBoundary(MeshBlock *pmb, int stage) {
  if (pmb->pcondbval->SendConductionBoundaryBuffers(pmb->phydro->phdif->dprl_tprl)==false)
    return TASK_FAIL;
  std::cout << "At SendConductionBoundary \n";
  return TASK_SUCCESS;
}

enum TaskStatus ConductionSolverTaskList::ReceiveConductionBoundary(MeshBlock *pmb, int stage) {
  if (pmb->pcondbval->ReceiveConductionBoundaryBuffers(pmb->phydro->phdif->dprl_tprl)==false)
    return TASK_FAIL;
  std::cout << "At ReceiveConductionBoundary \n";
  return TASK_SUCCESS;
}

enum TaskStatus ConductionSolverTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  pmb->pcondbval->ApplyPhysicalBoundaries();
  std::cout << "At PhysicalBoundary \n";
  return TASK_NEXT;
}
