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
  explicit FrictionSolverDetails(seissol::initializer::parameters::DRParameters* drParameters);
  ~FrictionSolverDetails() override;

  void allocateAuxiliaryMemory(GlobalData* globalData) override;

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
