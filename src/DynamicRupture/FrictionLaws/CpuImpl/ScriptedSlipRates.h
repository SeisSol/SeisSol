// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SCRIPTEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SCRIPTEDSLIPRATES_H_

#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/CpuImpl/ImposedSlipRates.h"
#include "DynamicRupture/FrictionLaws/CpuImpl/SourceTimeFunction.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/SlipRateScript.h"
#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Parallel/Runtime/Stream.h"
#include "Reader/Scripting/DataTable.h"

#include <memory>
#include <utility>

namespace seissol::dr::friction_law::cpu {

/**
 * The imposed slip rates of a script (FL 36). Before the law evaluates a step of its layer, the
 * script gives the slip rates of all sub-steps at all points of the layer; the law imposes them
 * as FL 33 imposes those of its time function.
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
    SlipRateEvaluator::LayerData data;
    data.parameters = reinterpret_cast<const double*>(
        layerData.var<LTSImposedSlipRatesScript::ScriptParameters>(Cfg()));
    data.slipRates = layerData.var<LTSImposedSlipRatesScript::ScriptSlipRates>(Cfg());
    evaluator_.setLayer(layerData.size() * misc::NumPaddedPoints<Cfg>, data, data);
  }

  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    evaluator_.evaluate(fullUpdateTime, frictionTime.deltaT, nullptr, false);
    Base::evaluate(fullUpdateTime, frictionTime, timeWeights, runtime);
  }

  private:
  SlipRateEvaluator evaluator_;
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SCRIPTEDSLIPRATES_H_
