// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SCRIPTEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SCRIPTEDSLIPRATES_H_

#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/ImposedSlipRates.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/SourceTimeFunction.h"
#include "DynamicRupture/FrictionLaws/SlipRateScript.h"
#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Runtime/Stream.h"
#include "Reader/Scripting/DataTable.h"

#include <memory>
#include <utility>

namespace seissol::dr::friction_law::gpu {

/**
 * The imposed slip rates of a script (FL 36) on a device. Before the friction kernel of a step,
 * a kernel of the script writes the slip rates of all sub-steps at all points of the layer, one
 * launch per sub-step on the stream of the friction kernel, at the times of the sub-steps from the
 * clock of the cluster on the device: as the friction kernel, so that both go into the graph of
 * the step. Where the build has no runtime compiler for its device, the script runs on the host
 * and its slip rates are copied over.
 */
template <typename Cfg>
class ScriptedSlipRates : public ImposedSlipRates<Cfg, ScriptedSTF<Cfg>> {
  public:
  using Base = ImposedSlipRates<Cfg, ScriptedSTF<Cfg>>;

  ScriptedSlipRates(const FrictionLawParameters<Real<Cfg>>& drParameters,
                    std::shared_ptr<const SlipRateScript> script)
      : Base(drParameters), evaluator_(std::move(script),
                                       reader::scripting::DataTypeTraits<Real<Cfg>>::Type,
                                       misc::TimeSteps<Cfg>) {}

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& runtime) override {
    Base::setupLayer(layerData, runtime);
    layer_ = &layerData;
    const auto dataAt = [&](seissol::initializer::AllocationPlace place) {
      SlipRateEvaluator::LayerData data;
      data.parameters = reinterpret_cast<const double*>(
          layerData.var<LTSImposedSlipRatesScript::ScriptParameters>(Cfg(), place));
      data.slipRates = layerData.var<LTSImposedSlipRatesScript::ScriptSlipRates>(Cfg(), place);
      return data;
    };
    evaluator_.setLayer(layerData.size() * misc::NumPaddedPoints<Cfg>,
                        dataAt(seissol::initializer::AllocationPlace::Host),
                        dataAt(seissol::initializer::AllocationPlace::Device));
    // before the time stepping, which records the work of a step
    evaluator_.prepareDevice(runtime.stream());
  }

  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (!evaluator_.evaluate(
            fullUpdateTime, frictionTime.deltaT, runtime.stream(), true, this->clock_)) {
      layer_->varSynchronizeTo<LTSImposedSlipRatesScript::ScriptSlipRates>(
          seissol::initializer::AllocationPlace::Device, runtime.stream());
    }
    Base::evaluate(fullUpdateTime, frictionTime, timeWeights, runtime);
  }

  private:
  SlipRateEvaluator evaluator_;
  DynamicRupture::Layer* layer_{nullptr};
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SCRIPTEDSLIPRATES_H_
