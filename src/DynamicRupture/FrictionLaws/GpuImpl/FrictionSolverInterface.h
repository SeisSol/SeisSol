#ifndef SEISSOL_FRICTION_SOLVER_INTERFACE_H
#define SEISSOL_FRICTION_SOLVER_INTERFACE_H

#include "DynamicRupture/Parameters.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

// A sycl-independent interface is required for interacting with the wp solver
// which, in its turn, is not supposed to know anything about SYCL
namespace seissol::dr::friction_law::gpu {
class FrictionSolverInterface : public seissol::dr::friction_law::FrictionSolver {
  public:
  explicit FrictionSolverInterface(dr::DRParameters* drParameters)
      : seissol::dr::friction_law::FrictionSolver(drParameters) {}
  ~FrictionSolverInterface() override{};

  virtual void initSyclQueue() = 0;
  void setMaxClusterSize(size_t size) { maxClusterSize = size; }
  virtual void allocateAuxiliaryMemory() = 0;
  virtual void copyStaticDataToDevice() = 0;

  protected:
  size_t maxClusterSize{};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_FRICTION_SOLVER_INTERFACE_H
