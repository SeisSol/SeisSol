// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_FRICTIONSOLVERHOST_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_FRICTIONSOLVERHOST_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Parallel/Runtime/HostRuntime.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Parallel/Trace.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Physics/TPCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawDriver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "GeneratedCode/init.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <algorithm>
#include <cstddef>
#include <utils/logger.h>

namespace seissol::dr::friction_law::cpu {

/**
 * The friction laws that are written once for the host and the device (see FrictionLawDriver), on
 * the host: the law `Derived` evaluates the faces of a layer one after the other, all points of a
 * face at once.
 */
template <typename Cfg, typename Derived>
class FrictionSolverHost : public FrictionSolverImpl<Cfg> {
  private:
  // only the CRTP subclass `Derived` constructs this base
  explicit FrictionSolverHost(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverImpl<Cfg>(drParameters), tpGridPoints_(drParameters.tpGridPoints),
        tpInverseFourierCoefficients_(drParameters.tpGridPoints),
        tpHeatSource_(drParameters.tpGridPoints) {}
  friend Derived;

  public:
  using real = Real<Cfg>;

  void allocateAuxiliaryMemory(GlobalData<Cfg>* globalData) override {
    FrictionSolverImpl<Cfg>::allocateAuxiliaryMemory(globalData);
    resampleMatrix_ = globalData->*init::resample<Cfg>::PoolMember;
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& /*runtime*/) override {
    layerSize_ = layerData.size();
    copyCommonStorageToLocal(&data_, layerData, seissol::initializer::AllocationPlace::Host);
    Derived::copySpecificStorageDataToLocal(
        &data_, layerData, seissol::initializer::AllocationPlace::Host);
    data_.drParameters = this->drParameters_;
  }

  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& /*runtime*/) override {
    if (layerSize_ == 0) {
      return;
    }

    if constexpr (model::MaterialOf<Cfg>::SupportsDR) {
      FrictionLawArgs<Cfg> args{};
      args.data = &data_;
      args.spaceWeights = this->spaceWeights_;
      args.resampleMatrix = resampleMatrix_;
      args.tpGridPoints = tpGridPoints_.data().data();
      args.tpInverseFourierCoefficients = tpInverseFourierCoefficients_.data().data();
      args.heatSource = tpHeatSource_.data().data();
      std::copy_n(timeWeights, misc::TimeSteps<Cfg>, args.timeWeights);
      std::copy_n(frictionTime.deltaT.data(), misc::TimeSteps<Cfg>, args.deltaT);
      args.fullUpdateTime = fullUpdateTime;

      // loop over all dynamic rupture faces, in this LTS layer
      static const auto TraceName = monitoring::Trace::name("friction law");
      const parallel::runtime::HostRuntime::LoopName loopName(TraceName);
      parallel::runtime::HostRuntime::forEach(layerSize_, [&](std::size_t ltsFace) {
        FrictionLawContext<Cfg, common::RangeType::CPU> ctx{};
        ctx.ltsFace = ltsFace;
        ctx.data = &data_;
        ctx.args = &args;
        ctx.fullUpdateTime = args.fullUpdateTime;
        FrictionLawDriver<Cfg, Derived>::evaluate(ctx);
      });
    } else {
      logError() << "The material" << model::MaterialOf<Cfg>::Text
                 << "does not support DR friction law computations.";
    }
  }

  private:
  std::size_t layerSize_{};
  const real* resampleMatrix_{nullptr};
  FrictionLawData<Cfg> data_{};
  /// the wavenumber grid of the thermal pressurization and the coefficients on it
  tp::GridPoints<real> tpGridPoints_;
  tp::InverseFourierCoefficients<real> tpInverseFourierCoefficients_;
  tp::GaussianHeatSource<real> tpHeatSource_;
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_FRICTIONSOLVERHOST_H_
