// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_

#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/TPCommon.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"

#include <yaml-cpp/yaml.h>

namespace seissol::dr::friction_law::gpu {

template <typename Cfg>
class FrictionSolverDetails : public FrictionSolverInterface<Cfg> {
  public:
  using real = Real<Cfg>;

  explicit FrictionSolverDetails(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverInterface<Cfg>(drParameters) {}

  ~FrictionSolverDetails() override = default;

  void allocateAuxiliaryMemory(GlobalData<Cfg>* globalData) override {
    // call the device module directly here
    {
#ifdef ACL_DEVICE
      data_ = reinterpret_cast<FrictionLawData<Cfg>*>(
          device::DeviceInstance::instance().api().allocGlobMem(sizeof(FrictionLawData<Cfg>)));
#endif
    }

    resampleMatrix_ = globalData->*init::resample<Cfg>::PoolMember;
    devSpaceWeights_ = globalData->*init::quadweights<Cfg>::PoolMember;

#ifdef ACL_DEVICE
    // The thermal-pressurization tables are functions of the grid alone, and
    // only the device path reads them -- the CPU friction law holds its own
    // copies. So they are built and uploaded here, alongside this solver's
    // other device memory, rather than travelling through the global
    // matrices. Per solver rather than per process: they live and die with
    // the device allocation they sit next to, and are rebuilt whenever it is.
    const auto upload = [](const auto& source) {
      auto& device = device::DeviceInstance::instance();
      const std::size_t bytes = source.data().size() * sizeof(real);
      auto* target = reinterpret_cast<real*>(device.api().allocGlobMem(bytes));
      device.api().copyTo(target, source.data().data(), bytes);
      return target;
    };
    devTpGridPoints_ = upload(tp::GridPoints<misc::NumTpGridPoints, real>());
    devTpInverseFourierCoefficients_ =
        upload(tp::InverseFourierCoefficients<misc::NumTpGridPoints, real>());
    devHeatSource_ = upload(tp::GaussianHeatSource<misc::NumTpGridPoints, real>());
#endif
  }

  protected:
  size_t currLayerSize_{};

  const real* resampleMatrix_{nullptr};
  const real* devSpaceWeights_{nullptr};
  real* devTpInverseFourierCoefficients_{nullptr};
  real* devTpGridPoints_{nullptr};
  real* devHeatSource_{nullptr};

  FrictionLawData<Cfg>* data_{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERDETAILS_H_
