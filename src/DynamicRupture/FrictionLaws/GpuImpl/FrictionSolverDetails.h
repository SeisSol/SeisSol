#ifndef SEISSOL_FRICTION_SOLVER_DETAILS_H
#define SEISSOL_FRICTION_SOLVER_DETAILS_H

#include <yaml-cpp/yaml.h>
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include <CL/sycl.hpp>
#include <Parallel/Graph/SyclGraphCapturing.hpp>
#include <Parallel/Graph/GraphCapturing.hpp>

#ifndef __DPCPP_COMPILER
namespace sycl = cl::sycl;
#endif

namespace seissol::dr::friction_law::gpu {
class FrictionSolverDetails : public FrictionSolverInterface {
  public:
  explicit FrictionSolverDetails(seissol::initializer::parameters::DRParameters* drParameters);
  ~FrictionSolverDetails() override;

  void initSyclQueue() override;
  void allocateAuxiliaryMemory() override;
  void copyStaticDataToDevice() override;

  virtual void
      copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                     seissol::initializer::DynamicRupture const* const dynRup,
                                     real fullUpdateTime) = 0;

  protected:
  size_t currLayerSize{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::numPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::numPaddedPoints]{nullptr};
  real* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  real* devSpaceWeights{nullptr};

  sycl::device device;
  sycl::queue queue;

  std::unordered_map<real, seissol::parallel::ComputeGraph<seissol::parallel::SyclGraph>> graphs;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_FRICTION_SOLVER_DETAILS_H
