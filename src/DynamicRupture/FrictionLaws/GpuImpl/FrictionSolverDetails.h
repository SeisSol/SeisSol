// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"

#include <yaml-cpp/yaml.h>

namespace seissol::dr::friction_law::gpu {

class FrictionSolverDetails : public FrictionSolverInterface {
  public:
  explicit FrictionSolverDetails(const FrictionLawParameters& drParameters)
      : FrictionSolverInterface(drParameters) {}

  ~FrictionSolverDetails() override = default;

  void allocateAuxiliaryMemory(GlobalData* globalData) override {
    {
      data_ =
          seissol::memory::allocTyped<FrictionLawData>(1, 1, memory::Memkind::DeviceGlobalMemory);
    }

    resampleMatrix_ = globalData->resampleMatrix;
    devSpaceWeights_ = globalData->spaceWeights;
    devTpInverseFourierCoefficients_ = globalData->tpInverseFourierCoefficients;
    devHeatSource_ = globalData->heatSource;
    devTpGridPoints_ = globalData->tpGridPoints;
  }

  protected:
  size_t currLayerSize_{};

  real* resampleMatrix_{nullptr};
  real* devSpaceWeights_{nullptr};
  real* devTpInverseFourierCoefficients_{nullptr};
  real* devTpGridPoints_{nullptr};
  real* devHeatSource_{nullptr};

  FrictionLawData* data_{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
