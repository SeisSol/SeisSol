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

#include <Memory/GlobalData.h>

namespace seissol::dr::friction_law::gpu {

template <typename Cfg>
class FrictionSolverDetails : public FrictionSolverInterface<Cfg> {
  public:
  explicit FrictionSolverDetails(seissol::initializer::parameters::DRParameters* drParameters);
  ~FrictionSolverDetails() override;

  void allocateAuxiliaryMemory(const GlobalData& globalData) override;

  protected:
  size_t currLayerSize{};

  using real = Real<Cfg>;

  real* resampleMatrix{nullptr};
  real* devSpaceWeights{nullptr};
  real* devTpInverseFourierCoefficients{nullptr};
  real* devTpGridPoints{nullptr};
  real* devHeatSource{nullptr};

  FrictionLawData<Cfg>* data{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
