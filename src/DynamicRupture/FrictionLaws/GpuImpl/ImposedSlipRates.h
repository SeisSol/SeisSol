#ifndef SEISSOL_GPU_IMPOSEDSLIPRATES_H
#define SEISSOL_GPU_IMPOSEDSLIPRATES_H

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
/**
 * Slip rates are set fixed values
 */
template <typename STF>
class ImposedSlipRates : public BaseFrictionSolver<ImposedSlipRates<STF>> {
  public:
  using BaseFrictionSolver<ImposedSlipRates>::BaseFrictionSolver;

  void copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                      seissol::initializer::DynamicRupture const* const dynRup,
                                      real fullUpdateTime) override {
    auto* concreteLts =
        dynamic_cast<seissol::initializer::LTSImposedSlipRates const* const>(dynRup);
    imposedSlipDirection1 = layerData.var(concreteLts->imposedSlipDirection1);
    imposedSlipDirection2 = layerData.var(concreteLts->imposedSlipDirection2);
    stf.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void allocateAuxiliaryMemory() override { FrictionSolverDetails::allocateAuxiliaryMemory(); }

  void updateFrictionAndSlip(unsigned timeIndex) {

    this->updateFrictionAndSlipImposedSR(this->faultStresses, this->tractionResults, timeIndex);
  }

  void updateFrictionAndSlipImposedSR(FaultStresses* devFaultStresses,
                                      TractionResults* devTractionResults,
                                      unsigned timeIndex) {

    auto timeIncrement{this->deltaT[timeIndex]};
    auto devDeltaT{this->deltaT};
    auto currentTime{this->mFullUpdateTime};
    auto* devImpAndEta{this->impAndEta};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devSlipRate1{this->slipRate1};
    auto* devSlipRate2{this->slipRate2};
    auto* devTraction1{this->traction1};
    auto* devTraction2{this->traction2};
    auto* devSlip1{this->slip1};
    auto* devSlip2{this->slip2};
    auto* devImposedSlipDirection1{this->imposedSlipDirection1};
    auto* devImposedSlipDirection2{this->imposedSlipDirection2};
    auto* devAccumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto devStf{this->stf};

    for (unsigned i = 0; i <= timeIndex; i++) {
      currentTime += devDeltaT[i];
    }

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = devFaultStresses[ltsFace];
        auto& tractionResults = devTractionResults[ltsFace];

        const real stfEvaluated = devStf.evaluate(currentTime, timeIncrement, ltsFace, pointIndex);

        devSlipRate1[ltsFace][pointIndex] =
            devImposedSlipDirection1[ltsFace][pointIndex] * stfEvaluated;
        devSlipRate2[ltsFace][pointIndex] =
            devImposedSlipDirection2[ltsFace][pointIndex] * stfEvaluated;
        devSlipRateMagnitude[ltsFace][pointIndex] =
            misc::magnitude(devSlipRate1[ltsFace][pointIndex], devSlipRate2[ltsFace][pointIndex]);

        // calculate traction
        devTraction1[timeIndex][pointIndex] =
            faultStresses.traction1[timeIndex][pointIndex] -
            devImpAndEta[ltsFace].etaS * devSlipRate1[ltsFace][pointIndex];
        devTraction2[timeIndex][pointIndex] =
            faultStresses.traction2[timeIndex][pointIndex] -
            devImpAndEta[ltsFace].etaS * devSlipRate2[ltsFace][pointIndex];

        // Update slip
        devSlip1[ltsFace][pointIndex] += devSlipRate1[ltsFace][pointIndex] * timeIncrement;
        devSlip2[ltsFace][pointIndex] += devSlipRate2[ltsFace][pointIndex] * timeIncrement;
        devAccumulatedSlipMagnitude[ltsFace][pointIndex] +=
            devSlipRateMagnitude[ltsFace][pointIndex] * timeIncrement;

        tractionResults.traction1[timeIndex][pointIndex] = devTraction1[ltsFace][pointIndex];
        tractionResults.traction2[timeIndex][pointIndex] = devTraction2[ltsFace][pointIndex];
      });
    });
  }

  void preHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {}
  void postHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {}
  void saveDynamicStressOutput() {}

  protected:
  real (*imposedSlipDirection1)[misc::numPaddedPoints];
  real (*imposedSlipDirection2)[misc::numPaddedPoints];
  STF stf{};
};

} // namespace seissol::dr::friction_law::gpu
#endif // SEISSOL_IMPOSEDSLIPRATES_H
