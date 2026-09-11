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

    resampleMatrix_ = globalData->resampleMatrix;
    devSpaceWeights_ = globalData->spaceWeights;

    const auto& tables = thermalPressurizationTables();
    devTpInverseFourierCoefficients_ = tables.inverseFourierCoefficients;
    devHeatSource_ = tables.heatSource;
    devTpGridPoints_ = tables.gridPoints;
  }

  protected:
  //! Device copies of the thermal-pressurization tables.
  struct TpTables {
    real* gridPoints{nullptr};
    real* inverseFourierCoefficients{nullptr};
    real* heatSource{nullptr};
  };

  /**
   * The tables are functions of the grid alone, and only the device path reads
   * them -- the CPU friction law holds its own copies. So they are built and
   * uploaded here instead of travelling through the global matrices, once for
   * the process rather than once per solver.
   */
  static const TpTables& thermalPressurizationTables() {
    static const TpTables tables = [] {
      TpTables result;
#ifdef ACL_DEVICE
      auto& device = device::DeviceInstance::getInstance();
      const auto upload = [&device](const auto& source) {
        const std::size_t bytes = source.data().size() * sizeof(real);
        auto* target = reinterpret_cast<real*>(device.api->allocGlobMem(bytes));
        device.api->copyTo(target, source.data().data(), bytes);
        return target;
      };
      result.gridPoints = upload(tp::GridPoints<misc::NumTpGridPoints>());
      result.inverseFourierCoefficients =
          upload(tp::InverseFourierCoefficients<misc::NumTpGridPoints>());
      result.heatSource = upload(tp::GaussianHeatSource<misc::NumTpGridPoints>());
#endif
      return result;
    }();
    return tables;
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
