#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <Initializer/Parameters/ModelParameters.h>
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

      using namespace seissol::dr::misc::quantity_indices;
      auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus[ltsFace]));
      auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus[ltsFace]));
      auto* qStressIPlus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedPlus));
      auto* qStressIMinus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedMinus));

      real lambda0P = impAndEta[ltsFace].lambda0P;
      real mu0P = impAndEta[ltsFace].mu0P;
      real rho0P = impAndEta[ltsFace].rho0P;
      real lambda0M = impAndEta[ltsFace].lambda0M;
      real mu0M = impAndEta[ltsFace].mu0M;
      real rho0M = impAndEta[ltsFace].rho0M;

      // TODO(NONLINEAR) What are these values?
      const real aB0 = damagedElasticParameters->aB0;
      const real aB1 = damagedElasticParameters->aB1;
      const real aB2 = damagedElasticParameters->aB2;
      const real aB3 = damagedElasticParameters->aB3;

      for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
        for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {

          real EspIp = (qIPlus[o][XX][i]) + (qIPlus[o][YY][i]) + (qIPlus[o][ZZ][i]);
          real EspIIp = (qIPlus[o][XX][i]) * (qIPlus[o][XX][i]) +
                        (qIPlus[o][YY][i]) * (qIPlus[o][YY][i]) +
                        (qIPlus[o][ZZ][i]) * (qIPlus[o][ZZ][i]) +
                        2 * (qIPlus[o][XY][i]) * (qIPlus[o][XY][i]) +
                        2 * (qIPlus[o][YZ][i]) * (qIPlus[o][YZ][i]) +
                        2 * (qIPlus[o][XZ][i]) * (qIPlus[o][XZ][i]);
          real alphap = qIPlus[o][DAM][i];
          real xip;
          if (EspIIp > 1e-30) {
            xip = EspIp / std::sqrt(EspIIp);
          } else {
            xip = 0.0;
          }

          // damage stress impAndEtaGet->gammaRP, mu0P
          real mu_eff = mu0P - alphap * impAndEta[ltsFace].gammaRP * impAndEta[ltsFace].xi0P -
                        0.5 * alphap * impAndEta[ltsFace].gammaRP * xip;
          real sxx_sp = lambda0P * EspIp - alphap * impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp) +
                        2 * mu_eff * (qIPlus[o][XX][i]);
          real syy_sp = lambda0P * EspIp - alphap * impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp) +
                        2 * mu_eff * (qIPlus[o][YY][i]);
          real szz_sp = lambda0P * EspIp - alphap * impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp) +
                        2 * mu_eff * (qIPlus[o][ZZ][i]);

          const real sxy_sp = 2 * mu_eff * (qIPlus[o][XY][i]);
          const real syz_sp = 2 * mu_eff * (qIPlus[o][YZ][i]);
          const real szx_sp = 2 * mu_eff * (qIPlus[o][XZ][i]);

          // breakage stress
          const real sxx_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                              (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o][XX][i]);
          const real syy_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                              (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o][YY][i]);
          const real szz_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                              (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o][ZZ][i]);

          const real sxy_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o][XY][i]);
          const real syz_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o][YZ][i]);
          const real szx_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o][XZ][i]);

          qStressIPlus[o][XX][i] = (1 - qIPlus[o][BRE][i]) * sxx_sp + qIPlus[o][BRE][i] * sxx_bp;

          qStressIPlus[o][YY][i] = (1 - qIPlus[o][BRE][i]) * syy_sp + qIPlus[o][BRE][i] * syy_bp;

          qStressIPlus[o][ZZ][i] = (1 - qIPlus[o][BRE][i]) * szz_sp + qIPlus[o][BRE][i] * szz_bp;

          qStressIPlus[o][XY][i] = (1 - qIPlus[o][BRE][i]) * sxy_sp + qIPlus[o][BRE][i] * sxy_bp;

          qStressIPlus[o][YZ][i] = (1 - qIPlus[o][BRE][i]) * syz_sp + qIPlus[o][BRE][i] * syz_bp;

          qStressIPlus[o][XZ][i] = (1 - qIPlus[o][BRE][i]) * szx_sp + qIPlus[o][BRE][i] * szx_bp;

          const real EspIm = (qIMinus[o][XX][i]) + (qIMinus[o][YY][i]) + (qIMinus[o][ZZ][i]);
          const real EspIIm = (qIMinus[o][XX][i]) * (qIMinus[o][XX][i]) +
                              (qIMinus[o][YY][i]) * (qIMinus[o][YY][i]) +
                              (qIMinus[o][ZZ][i]) * (qIMinus[o][ZZ][i]) +
                              2 * (qIMinus[o][XY][i]) * (qIMinus[o][XY][i]) +
                              2 * (qIMinus[o][YZ][i]) * (qIMinus[o][YZ][i]) +
                              2 * (qIMinus[o][XZ][i]) * (qIMinus[o][XZ][i]);
          real alpham = qIMinus[o][DAM][i];
          real xim;
          if (EspIIm > 1e-30) {
            xim = EspIm / std::sqrt(EspIIm);
          } else {
            xim = 0.0;
          }

          // damage stress minus
          mu_eff = mu0M - alpham * impAndEta[ltsFace].gammaRM * impAndEta[ltsFace].xi0M -
                   0.5 * alpham * impAndEta[ltsFace].gammaRM * xim;
          const real sxx_sm = lambda0M * EspIm -
                              alpham * impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm) +
                              2 * mu_eff * (qIMinus[o][XX][i]);
          const real syy_sm = lambda0M * EspIm -
                              alpham * impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm) +
                              2 * mu_eff * (qIMinus[o][YY][i]);
          const real szz_sm = lambda0M * EspIm -
                              alpham * impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm) +
                              2 * mu_eff * (qIMinus[o][ZZ][i]);

          const real sxy_sm = 2 * mu_eff * (qIMinus[o][XY][i]);
          const real syz_sm = 2 * mu_eff * (qIMinus[o][YZ][i]);
          const real szx_sm = 2 * mu_eff * (qIMinus[o][XZ][i]);

          // breakage stress
          const real sxx_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                              (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o][XX][i]);
          const real syy_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                              (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o][YY][i]);
          const real szz_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                              (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o][ZZ][i]);

          const real sxy_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o][XY][i]);
          const real syz_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o][YZ][i]);
          const real szx_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o][XZ][i]);

          qStressIMinus[o][XX][i] = (1 - qIMinus[o][BRE][i]) * sxx_sm + qIMinus[o][BRE][i] * sxx_bm;

          qStressIMinus[o][YY][i] = (1 - qIMinus[o][BRE][i]) * syy_sm + qIMinus[o][BRE][i] * syy_bm;

          qStressIMinus[o][ZZ][i] = (1 - qIMinus[o][BRE][i]) * szz_sm + qIMinus[o][BRE][i] * szz_bm;

          qStressIMinus[o][XY][i] = (1 - qIMinus[o][BRE][i]) * sxy_sm + qIMinus[o][BRE][i] * sxy_bm;

          qStressIMinus[o][YZ][i] = (1 - qIMinus[o][BRE][i]) * syz_sm + qIMinus[o][BRE][i] * syz_bm;

          qStressIMinus[o][XZ][i] = (1 - qIMinus[o][BRE][i]) * szx_sm + qIMinus[o][BRE][i] * szx_bm;

          qStressIPlus[o][U][i] = qIPlus[o][U][i];
          qStressIPlus[o][V][i] = qIPlus[o][V][i];
          qStressIPlus[o][W][i] = qIPlus[o][W][i];
          qStressIPlus[o][DAM][i] = qIPlus[o][DAM][i];
          qStressIPlus[o][BRE][i] = qIPlus[o][BRE][i];

          qStressIMinus[o][U][i] = qIMinus[o][U][i];
          qStressIMinus[o][V][i] = qIMinus[o][V][i];
          qStressIMinus[o][W][i] = qIMinus[o][W][i];
          qStressIMinus[o][DAM][i] = qIMinus[o][DAM][i];
          qStressIMinus[o][BRE][i] = qIMinus[o][BRE][i];
        }
      } // time integration loop
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
