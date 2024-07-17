#ifndef SEISSOL_GPU_DRCODEGEN_H
#define SEISSOL_GPU_DRCODEGEN_H

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"

#include <cstddef>

#include "init.h"
#include "kernel.h"
#include "tensor.h"
#include <device.h>

namespace seissol::dr::friction_law::gpu {
class FrictionSolverAdapter : public FrictionSolverInterface {
  public:
  explicit FrictionSolverAdapter(seissol::initializer::parameters::DRParameters* drParameters, int lawId)
    : FrictionSolverInterface(drParameters), lawId(lawId)
  {

  }
  ~FrictionSolverAdapter() override = default;

  void initSyclQueue() override {}
  void allocateAuxiliaryMemory() override {
    spaceWeights = reinterpret_cast<real*>(device::DeviceInstance::getInstance().api->allocGlobMem(sizeof(real) * tensor::quadweights::size()));
    device::DeviceInstance::getInstance().api->copyTo(spaceWeights, init::quadweights::Values, sizeof(real) * tensor::quadweights::size());
  }
  void copyStaticDataToDevice() override {

  }

  virtual void
      copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
}

void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const double timeWeights[CONVERGENCE_ORDER],
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (lawId <= 0) {
        // noFL
        return;
    }

    const auto place = allocationPlace();

    FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.getNumberOfCells();

    seissol::dynamicRupture::kernel::gpu_frictionLaw frictionLawKernel;

    // frictionLawKernel.resampleMatrix = resampleMatrix;
    frictionLawKernel.spaceWeights = spaceWeights;

    double sumDt = 0;
    for (int i = 0; i < CONVERGENCE_ORDER; ++i) {
        frictionLawKernel.timeWeights(i) = timeWeights[i];
        frictionLawKernel.dt(i) = deltaT[i];
        sumDt += deltaT[i];
    }
    frictionLawKernel.t0 = this->drParameters->t0;
    frictionLawKernel.sumDt = sumDt;
    frictionLawKernel.fullUpdateTime = fullUpdateTime;

    frictionLawKernel.impAndEta = reinterpret_cast<real*>(layerData.var(dynRup->impAndEta, place));
    // not yet: only needed for poroelasticity
    // frictionLawKernel.impedanceMatrices = layerData.var(dynRup->impedanceMatrices, place);
    frictionLawKernel.initialStressInFaultCS = reinterpret_cast<real*>(layerData.var(dynRup->initialStressInFaultCS, place));
    frictionLawKernel.nucleationStressInFaultCS = reinterpret_cast<real*>(layerData.var(dynRup->nucleationStressInFaultCS, place));
    // not yet: only for LSW, R+S
    // frictionLawKernel.mu = reinterpret_cast<real*>(layerData.var(dynRup->mu, place));
    frictionLawKernel.accumulatedSlipMagnitude = reinterpret_cast<real*>(layerData.var(dynRup->accumulatedSlipMagnitude, place));
    frictionLawKernel.slip1 = reinterpret_cast<real*>(layerData.var(dynRup->slip1, place));
    frictionLawKernel.slip2 = reinterpret_cast<real*>(layerData.var(dynRup->slip2, place));
    frictionLawKernel.slipRateMagnitude = reinterpret_cast<real*>(layerData.var(dynRup->slipRateMagnitude, place));
    frictionLawKernel.slipRate1 = reinterpret_cast<real*>(layerData.var(dynRup->slipRate1, place));
    frictionLawKernel.slipRate2 = reinterpret_cast<real*>(layerData.var(dynRup->slipRate2, place));
    frictionLawKernel.ruptureTime = reinterpret_cast<real*>(layerData.var(dynRup->ruptureTime, place));
    frictionLawKernel.ruptureTimePending = reinterpret_cast<bool*>(layerData.var(dynRup->ruptureTimePending, place));
    frictionLawKernel.peakSlipRate = reinterpret_cast<real*>(layerData.var(dynRup->peakSlipRate, place));
    frictionLawKernel.traction1 = reinterpret_cast<real*>(layerData.var(dynRup->traction1, place));
    frictionLawKernel.traction2 = reinterpret_cast<real*>(layerData.var(dynRup->traction2, place));
    frictionLawKernel.imposedStatePlus = reinterpret_cast<real*>(layerData.var(dynRup->imposedStatePlus, place));
    frictionLawKernel.imposedStateMinus = reinterpret_cast<real*>(layerData.var(dynRup->imposedStateMinus, place));
    frictionLawKernel.accumulatedSlip = reinterpret_cast<real*>(layerData.var(dynRup->drEnergyOutput, place));
    frictionLawKernel.frictionalEnergy = reinterpret_cast<real*>(layerData.var(dynRup->drEnergyOutput, place));
    frictionLawKernel.energyDataSlip = reinterpret_cast<real*>(layerData.var(dynRup->drEnergyOutput, place));
    // frictionLawKernel.timeSinceSlipRateBelowThreshold = reinterpret_cast<real*>(layerData.var(dynRup->drEnergyOutput, place));
    // frictionLawKernel.dynStressTime = reinterpret_cast<real*>(layerData.var(dynRup->dynStressTime, place));
    // frictionLawKernel.dynStressTimePending = reinterpret_cast<bool*>(layerData.var(dynRup->dynStressTimePending, place));
    frictionLawKernel.qInterpolatedPlus = reinterpret_cast<real*>(layerData.var(dynRup->qInterpolatedPlus, place));
    frictionLawKernel.qInterpolatedMinus = reinterpret_cast<real*>(layerData.var(dynRup->qInterpolatedMinus, place));
    frictionLawKernel.initialPressure = reinterpret_cast<real*>(layerData.var(dynRup->initialPressure, place));
    frictionLawKernel.nucleationPressure = reinterpret_cast<real*>(layerData.var(dynRup->nucleationPressure, place));

    frictionLawKernel.doubledSurfaceArea = reinterpret_cast<double*>(layerData.var(dynRup->godunovData, place));

    // TODO: make more OOP-y

    /*if ((dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup);
        frictionLawKernel.dC = layerData.var(specializedDynRup->dC, place);
        frictionLawKernel.muS = layerData.var(specializedDynRup->muS, place);
        frictionLawKernel.muD = layerData.var(specializedDynRup->muD, place);
        frictionLawKernel.cohesion = layerData.var(specializedDynRup->cohesion, place);
        frictionLawKernel.forcedRuptureTime = layerData.var(specializedDynRup->forcedRuptureTime, place);
    }
    if ((dynamic_cast<const seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup);
        frictionLawKernel.regularisedStrength = layerData.var(specializedDynRup->regularisedStrength, place);
    }*/
    if ((dynamic_cast<const seissol::initializer::LTSImposedSlipRates*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSImposedSlipRates*>(dynRup);
        frictionLawKernel.imposedSlipDirection1 = reinterpret_cast<real*>(layerData.var(specializedDynRup->imposedSlipDirection1, place));
        frictionLawKernel.imposedSlipDirection2 = reinterpret_cast<real*>(layerData.var(specializedDynRup->imposedSlipDirection2, place));
        frictionLawKernel.onsetTime = reinterpret_cast<real*>(layerData.var(specializedDynRup->onsetTime, place));
    }
    if ((dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe*>(dynRup);
        frictionLawKernel.tauS = reinterpret_cast<real*>(layerData.var(specializedDynRup->tauS, place));
        frictionLawKernel.tauR = reinterpret_cast<real*>(layerData.var(specializedDynRup->tauR, place));
    }
    if ((dynamic_cast<const seissol::initializer::LTSImposedSlipRatesGaussian*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSImposedSlipRatesGaussian*>(dynRup);
        frictionLawKernel.riseTime = reinterpret_cast<real*>(layerData.var(specializedDynRup->riseTime, place));
    }
    /*if ((dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);
        frictionLawKernel.rsA = layerData.var(specializedDynRup->rsA, place);
        frictionLawKernel.rsSl0 = layerData.var(specializedDynRup->rsSl0, place);
        frictionLawKernel.stateVariable = layerData.var(specializedDynRup->stateVariable, place);
    }
    if ((dynamic_cast<const seissol::initializer::LTSRateAndStateFastVelocityWeakening*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSRateAndStateFastVelocityWeakening*>(dynRup);
        frictionLawKernel.rsSrW = layerData.var(specializedDynRup->rsSrW, place);
    }
    if ((dynamic_cast<const seissol::initializer::LTSRateAndStateThermalPressurization*>(dynRup)) != nullptr) {
        auto specializedDynRup = dynamic_cast<const seissol::initializer::LTSRateAndStateThermalPressurization*>(dynRup);
        frictionLawKernel.temperature = layerData.var(specializedDynRup->temperature, place);
        frictionLawKernel.pressure = layerData.var(specializedDynRup->pressure, place);
        frictionLawKernel.theta = layerData.var(specializedDynRup->theta, place);
        frictionLawKernel.sigma = layerData.var(specializedDynRup->sigma, place);
        frictionLawKernel.thetaTmpBuffer = layerData.var(specializedDynRup->thetaTmpBuffer, place);
        frictionLawKernel.sigmaTmpBuffer = layerData.var(specializedDynRup->sigmaTmpBuffer, place);
        frictionLawKernel.faultStrength = layerData.var(specializedDynRup->faultStrength, place);
        frictionLawKernel.halfWidthShearZone = layerData.var(specializedDynRup->halfWidthShearZone, place);
        frictionLawKernel.hydraulicDiffusivity = layerData.var(specializedDynRup->hydraulicDiffusivity, place);
    }*/

    frictionLawKernel.numElements = currLayerSize;
    frictionLawKernel.streamPtr = runtime.stream();
    frictionLawKernel.execute(lawId);
}

  protected:
  size_t currLayerSize{};

  real* resampleMatrix{nullptr};
  real* spaceWeights{nullptr};

  double* devTimeWeights{nullptr};
  int lawId;
};
} // namespace seissol::dr::friction_law::gpu


#endif
