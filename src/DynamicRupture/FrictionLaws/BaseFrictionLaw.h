#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <Initializer/Parameters/ModelParameters.h>
#include <Kernels/Time.h>
#include <Kernels/precision.hpp>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "FrictionSolver.h"
#include "FrictionSolverCommon.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Monitoring/instrumentation.hpp"
#include "SeisSol.h"

namespace seissol::dr::friction_law {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  public:
  explicit BaseFrictionLaw(seissol::SeisSol& seissolInstance)
      : FrictionSolver(&seissolInstance.getSeisSolParameters().drParameters,
                       &seissolInstance.getSeisSolParameters().model.damagedElasticParameters){};

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializer::Layer& layerData,
                seissol::initializer::DynamicRupture const* const dynRup,
                real fullUpdateTime,
                const double timeWeights[CONVERGENCE_ORDER]) override {
    SCOREP_USER_REGION_DEFINE(myRegionHandle)
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      alignas(ALIGNMENT) FaultStresses faultStresses{};
      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePrecomputeStress");

      alignas(PAGESIZE_STACK)
          real qStressInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] =
              {{0.0}};
      alignas(PAGESIZE_STACK)
          real qStressInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] =
              {{0.0}};

#ifdef USE_DAMAGEDELASTIC
      // TODO: convert from strain to stress

      using QInterpolatedShapeT =
          const real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];
      using QStressInterpolatedShapeT =
          real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];

      auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus[ltsFace]));
      auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus[ltsFace]));
      auto* qStressIPlus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedPlus));
      auto* qStressIMinus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedMinus));
      seissol::kernels::Time m_timeKernel;
      m_timeKernel.setDamagedElasticParameters(damagedElasticParameters);

      m_timeKernel.computeNonLinearBaseFrictionLaw(
          impAndEta, ltsFace, *qIPlus[0], *qStressIPlus[0], *qIMinus[0], *qStressIMinus[0]);
#else
      for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
        for (unsigned i = 0; i < seissol::tensor::QInterpolated::size(); i++) {
          qStressInterpolatedPlus[o][i] = 0.0;
          qStressInterpolatedMinus[o][i] = 0.0;
        }
      }
#endif

      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                impedanceMatrices[ltsFace],
                                                qStressInterpolatedPlus,
                                                qStressInterpolatedMinus,
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace],
                                                damagedElasticParameters);
      LIKWID_MARKER_STOP("computeDynamicRupturePrecomputeStress");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePreHook");
      // define some temporary variables
      std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
      std::array<real, misc::numPaddedPoints> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
      LIKWID_MARKER_STOP("computeDynamicRupturePreHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRuptureUpdateFrictionAndSlip",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRuptureUpdateFrictionAndSlip");
      TractionResults tractionResults = {};

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        common::adjustInitialStress(initialStressInFaultCS[ltsFace],
                                    nucleationStressInFaultCS[ltsFace],
                                    initialPressure[ltsFace],
                                    nucleationPressure[ltsFace],
                                    this->mFullUpdateTime,
                                    this->drParameters->t0,
                                    this->deltaT[timeIndex]);

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
                                                   // TODO(NONLINEAR): unify
                                                   qStressInterpolatedPlus,
                                                   qStressInterpolatedMinus,
                                                   timeWeights);
      LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
      SCOREP_USER_REGION_END(myRegionHandle)

      if (this->drParameters->isFrictionEnergyRequired) {
        common::computeFrictionEnergy(energyData[ltsFace],
                                      qStressInterpolatedPlus,
                                      qStressInterpolatedMinus,
                                      impAndEta[ltsFace],
                                      timeWeights,
                                      spaceWeights,
                                      godunovData[ltsFace]);
      }
    }
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H
