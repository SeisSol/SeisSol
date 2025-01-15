// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include <CL/sycl.hpp>
#include <yaml-cpp/yaml.h>

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
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) = 0;

  protected:
  size_t currLayerSize{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::NumPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::NumPaddedPoints]{nullptr};
  real* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  real* devSpaceWeights{nullptr};

  sycl::device device;
  sycl::queue queue;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
