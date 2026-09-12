// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/TPCommon.h"
#include "DynamicRupture/Misc.h"

#include <yaml-cpp/yaml.h>

namespace seissol::dr::friction_law::gpu {

class FrictionSolverDetails : public FrictionSolverInterface {
  public:
  explicit FrictionSolverDetails(const FrictionLawParameters& drParameters)
      : FrictionSolverInterface(drParameters) {}

  ~FrictionSolverDetails() override = default;

  void allocateAuxiliaryMemory(GlobalData* globalData) override {
    // call the device module directly here
    {
#ifdef ACL_DEVICE
      data_ = reinterpret_cast<FrictionLawData*>(
          device::DeviceInstance::getInstance().api->allocGlobMem(sizeof(FrictionLawData)));
#endif
    }

    resampleMatrix_ = globalData->resample;
    devSpaceWeights_ = globalData->quadweights;

#ifdef ACL_DEVICE
    // The thermal-pressurization tables are functions of the grid alone, and
    // only the device path reads them -- the CPU friction law holds its own
    // copies. So they are built and uploaded here, alongside this solver's
    // other device memory, rather than travelling through the global
    // matrices. Per solver rather than per process: they live and die with
    // the device allocation they sit next to, and are rebuilt whenever it is.
    const auto upload = [](const auto& source) {
      auto& device = device::DeviceInstance::getInstance();
      const std::size_t bytes = source.data().size() * sizeof(real);
      auto* target = reinterpret_cast<real*>(device.api->allocGlobMem(bytes));
      device.api->copyTo(target, source.data().data(), bytes);
      return target;
    };
    devTpGridPoints_ = upload(tp::GridPoints<misc::NumTpGridPoints>());
    devTpInverseFourierCoefficients_ =
        upload(tp::InverseFourierCoefficients<misc::NumTpGridPoints>());
    devHeatSource_ = upload(tp::GaussianHeatSource<misc::NumTpGridPoints>());
#endif
  }

  size_t currLayerSize_{};

  const real* resampleMatrix_{nullptr};
  const real* devSpaceWeights_{nullptr};
  real* devTpInverseFourierCoefficients_{nullptr};
  real* devTpGridPoints_{nullptr};
  real* devHeatSource_{nullptr};

  FrictionLawData* data_{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
