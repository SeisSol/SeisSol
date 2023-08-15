#ifndef SEISSOL_FRICTION_SOLVER_DETAILS_H
#define SEISSOL_FRICTION_SOLVER_DETAILS_H

#include <yaml-cpp/yaml.h>
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include <CL/sycl.hpp>

#ifndef __DPCPP_COMPILER
namespace sycl = cl::sycl;
#endif

namespace seissol::dr::friction_law::gpu {
template <typename Config>
class FrictionSolverDetails : public FrictionSolverInterface {
  public:
  using RealT = typename Config::RealT;
  explicit FrictionSolverDetails(dr::DRParameters* drParameters);
  ~FrictionSolverDetails() override;

  void initSyclQueue() override;
  void allocateAuxiliaryMemory() override;
  void copyStaticDataToDevice() override;

  virtual void copySpecificLtsDataTreeToLocal(
      seissol::initializers::Layer& layerData,
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      RealT fullUpdateTime) = 0;

  protected:
  size_t currLayerSize{};

  FaultStresses<Config>* faultStresses{nullptr};
  TractionResults<Config>* tractionResults{nullptr};
  RealT (*stateVariableBuffer)[misc::numPaddedPoints<Config>]{nullptr};
  RealT (*strengthBuffer)[misc::numPaddedPoints<Config>]{nullptr};
  RealT* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  RealT* devSpaceWeights{nullptr};

  sycl::device device;
  sycl::queue queue;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_FRICTION_SOLVER_DETAILS_H
