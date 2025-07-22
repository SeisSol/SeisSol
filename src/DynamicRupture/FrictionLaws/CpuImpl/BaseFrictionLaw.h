// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Monitoring/Instrumentation.h"

namespace seissol::dr::friction_law::cpu {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  public:
  explicit BaseFrictionLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolver(drParameters) {}

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const double timeWeights[ConvergenceOrder],
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (layerData.size() == 0) {
      return;
    }

    SCOREP_USER_REGION_DEFINE(myRegionHandle)
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.size(); ++ltsFace) {
      alignas(Alignment) FaultStresses<Executor::Host> faultStresses{};
      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePrecomputeStress");
      const auto etaPDamp =
          drParameters->etaStop > this->mFullUpdateTime ? drParameters->etaHack : 1.0;
      
      #ifdef USE_DAMAGE
      // convert from strain to stress here
      //TODO: extract this into a separate function
        alignas(PagesizeStack) real qStressInterpolatedPlus[ConvergenceOrder][seissol::tensor::QInterpolated::size()] = {{0.0}};
        alignas(PagesizeStack) real qStressInterpolatedMinus[ConvergenceOrder][seissol::tensor::QInterpolated::size()] = {{0.0}};
        using QInterpolatedShapeT = const real(*)[seissol::dr::misc::NumQuantities][seissol::dr::misc::NumPaddedPoints];
        using QStressInterpolatedShapeT = real(*)[seissol::dr::misc::NumQuantities][seissol::dr::misc::NumPaddedPoints];
        auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(this->qInterpolatedPlus[ltsFace]);
        auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(this->qInterpolatedMinus[ltsFace]);
        auto* qStressIPlus = reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedPlus);
        auto* qStressIMinus = reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedMinus);

        using namespace seissol::dr::misc::quantity_indices;
        const unsigned DAM = 9;
        const unsigned BRE = 10;

        const real epsInitxx = damageParameters->epsInitxx;
        const real epsInityy = damageParameters->epsInityy;
        const real epsInitzz = damageParameters->epsInitzz;
        const real epsInitxy = damageParameters->epsInitxy;
        const real epsInityz = damageParameters->epsInityz;
        const real epsInitxz = damageParameters->epsInitxz;
        const real aB0 = damageParameters->aB0;
        const real aB1 = damageParameters->aB1;
        const real aB2 = damageParameters->aB2;
        const real aB3 = damageParameters->aB3;
        const real lamda0P = impAndEta[ltsFace].lambda0P;
        const real mu0P = impAndEta[ltsFace].mu0P;
        const real rho0P = impAndEta[ltsFace].rho0P;
        const real lambda0M = impAndEta[ltsFace].lambda0M;
        const real mu0M = impAndEta[ltsFace].mu0M;
        const real rho0M = impAndEta[ltsFace].rho0M;

        for(std::size_t o = 0; o < ConvergenceOrder ; o++){
          for(std::size_t i = 0; i < seissol::dr::misc::NumPaddedPoints; i++){
            real EspIp = (qIPlus[o][XX][i]+epsInitxx) + (qIPlus[o][YY][i]+epsInityy) + (qIPlus[o][ZZ][i]+epsInitzz);
            
            real EspIIp = (qIPlus[o][XX][i]+epsInitxx)*(qIPlus[o][XX][i]+epsInitxx)
            + (qIPlus[o][YY][i]+epsInityy)*(qIPlus[o][YY][i]+epsInityy)
            + (qIPlus[o][ZZ][i]+epsInitzz)*(qIPlus[o][ZZ][i]+epsInitzz)
            + 2*(qIPlus[o][XY][i]+epsInitxy)*(qIPlus[o][XY][i]+epsInitxy)
            + 2*(qIPlus[o][YZ][i]+epsInityz)*(qIPlus[o][YZ][i]+epsInityz)
            + 2*(qIPlus[o][XZ][i]+epsInitxz)*(qIPlus[o][XZ][i]+epsInitxz);
          real alphap = qIPlus[o][DAM][i];
          real xip;
          if (EspIIp > 1e-30){
            xip = EspIp / std::sqrt(EspIIp);
          } else{
            xip = 0.0;
          }
          }
        }

      #endif
      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                impedanceMatrices[ltsFace],
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace],
                                                etaPDamp);
      LIKWID_MARKER_STOP("computeDynamicRupturePrecomputeStress");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePreHook");
      // define some temporary variables
      std::array<real, misc::NumPaddedPoints> stateVariableBuffer{0};
      std::array<real, misc::NumPaddedPoints> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
      LIKWID_MARKER_STOP("computeDynamicRupturePreHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRuptureUpdateFrictionAndSlip",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRuptureUpdateFrictionAndSlip");
      TractionResults<Executor::Host> tractionResults = {};

      // loop over sub time steps (i.e. quadrature points in time)
      real updateTime = this->mFullUpdateTime;
      for (std::size_t timeIndex = 0; timeIndex < ConvergenceOrder; timeIndex++) {
        updateTime += this->deltaT[timeIndex];
        for (unsigned i = 0; i < this->drParameters->nucleationCount; ++i) {
          common::adjustInitialStress(initialStressInFaultCS[ltsFace],
                                      nucleationStressInFaultCS[i][ltsFace],
                                      initialPressure[ltsFace],
                                      nucleationPressure[i][ltsFace],
                                      updateTime,
                                      this->drParameters->t0[i],
                                      this->drParameters->s0[i],
                                      this->deltaT[timeIndex]);
        }

        static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses,
                                                           tractionResults,
                                                           stateVariableBuffer,
                                                           strengthBuffer,
                                                           ltsFace,
                                                           timeIndex);
      }
      LIKWID_MARKER_STOP("computeDynamicRuptureUpdateFrictionAndSlip");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePostHook", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePostHook");
      static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

      common::saveRuptureFrontOutput(ruptureTimePending[ltsFace],
                                     ruptureTime[ltsFace],
                                     slipRateMagnitude[ltsFace],
                                     mFullUpdateTime);

      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);
      LIKWID_MARKER_STOP("computeDynamicRupturePostHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRupturePostcomputeImposedState",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePostcomputeImposedState");
      common::postcomputeImposedStateFromNewStress(faultStresses,
                                                   tractionResults,
                                                   impAndEta[ltsFace],
                                                   impedanceMatrices[ltsFace],
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qInterpolatedPlus[ltsFace],
                                                   qInterpolatedMinus[ltsFace],
                                                   timeWeights);
      LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
      SCOREP_USER_REGION_END(myRegionHandle)

      if (this->drParameters->isFrictionEnergyRequired) {

        if (this->drParameters->isCheckAbortCriteraEnabled) {
          common::updateTimeSinceSlipRateBelowThreshold(
              slipRateMagnitude[ltsFace],
              ruptureTimePending[ltsFace],
              energyData[ltsFace],
              this->sumDt,
              this->drParameters->terminatorSlipRateThreshold);
        }
        common::computeFrictionEnergy(energyData[ltsFace],
                                      qInterpolatedPlus[ltsFace],
                                      qInterpolatedMinus[ltsFace],
                                      impAndEta[ltsFace],
                                      timeWeights,
                                      spaceWeights,
                                      godunovData[ltsFace],
                                      slipRateMagnitude[ltsFace],
                                      this->drParameters->energiesFromAcrossFaultVelocities);
      }
    }
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
