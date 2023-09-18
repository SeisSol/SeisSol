#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "FrictionSolver.h"
#include "FrictionSolverCommon.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol::dr::friction_law {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  public:
  explicit BaseFrictionLaw(dr::DRParameters* drParameters) : FrictionSolver(drParameters){};

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture const* const dynRup,
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
      /// TODO: convert from strain to stress
      alignas(PAGESIZE_STACK) real qStressInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};
      alignas(PAGESIZE_STACK) real qStressInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};

      using QInterpolatedShapeT = const real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];
      using QStressInterpolatedShapeT = real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];
      // std::cout << seissol::dr::misc::numQuantities << std::endl;
      auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus[ltsFace]));
      auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus[ltsFace]));
      auto* qStressIPlus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedPlus));
      auto* qStressIMinus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedMinus));

      using namespace seissol::dr::misc::quantity_indices;
      unsigned DAM = 9;

      real epsInitxx = -0e-2; // eps_xx0
      real epsInityy = -0e-1; // eps_yy0
      real epsInitzz = -0e-1; // eps_zz0
      real lambda0P = impAndEta[ltsFace].lambda0P;
      real mu0P = impAndEta[ltsFace].mu0P;
      real rho0P = impAndEta[ltsFace].rho0P;
      real lambda0M = impAndEta[ltsFace].lambda0M;
      real mu0M = impAndEta[ltsFace].mu0M;
      real rho0M = impAndEta[ltsFace].rho0M;

      for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
        for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints;
            i ++) {

          real EspIp = (qIPlus[o][XX][i]+epsInitxx) + (qIPlus[o][YY][i]+epsInityy) + (qIPlus[o][ZZ][i]+epsInitzz);
          real EspIIp = (qIPlus[o][XX][i]+epsInitxx)*(qIPlus[o][XX][i]+epsInitxx)
            + (qIPlus[o][YY][i]+epsInityy)*(qIPlus[o][YY][i]+epsInityy)
            + (qIPlus[o][ZZ][i]+epsInitzz)*(qIPlus[o][ZZ][i]+epsInitzz)
            + 2*qIPlus[o][XY][i]*qIPlus[o][XY][i]
            + 2*qIPlus[o][YZ][i]*qIPlus[o][YZ][i]
            + 2*qIPlus[o][XZ][i]*qIPlus[o][XZ][i];
          real alphap = qIPlus[o][DAM][i];
          real xip;
          if (EspIIp > 1e-30){
            xip = EspIp / std::sqrt(EspIIp);
          } else{
            xip = 0.0;
          }

          qStressIPlus[o][XX][i] = (lambda0P*EspIp - alphap*impAndEta[ltsFace].gammaRP*std::sqrt(EspIIp))
                + (2*(mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P)
                    - alphap*impAndEta[ltsFace].gammaRP*xip)
                  *(qIPlus[o][XX][i]+epsInitxx);

          qStressIPlus[o][YY][i] = (lambda0P*EspIp - alphap*impAndEta[ltsFace].gammaRP*std::sqrt(EspIIp))
                + (2*(mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P)
                    - alphap*impAndEta[ltsFace].gammaRP*xip)
                  *(qIPlus[o][YY][i]+epsInityy);

          qStressIPlus[o][ZZ][i] = (lambda0P*EspIp - alphap*impAndEta[ltsFace].gammaRP*std::sqrt(EspIIp))
                + (2*(mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P)
                    - alphap*impAndEta[ltsFace].gammaRP*xip)
                  *(qIPlus[o][ZZ][i]+epsInitzz);

          qStressIPlus[o][XY][i] = 0
                + (2*(mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P)
                    - alphap*impAndEta[ltsFace].gammaRP*xip)
                  *qIPlus[o][XY][i];

          qStressIPlus[o][YZ][i] = 0
                + (2*(mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P)
                    - alphap*impAndEta[ltsFace].gammaRP*xip)
                  *qIPlus[o][YZ][i];

          qStressIPlus[o][XZ][i] = 0
                + (2*(mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P)
                    - alphap*impAndEta[ltsFace].gammaRP*xip)
                  *qIPlus[o][XZ][i];

          real EspIm = (qIMinus[o][XX][i]+epsInitxx) + (qIMinus[o][YY][i]+epsInityy) + (qIMinus[o][ZZ][i]+epsInitzz);
          real EspIIm = (qIMinus[o][XX][i]+epsInitxx)*(qIMinus[o][XX][i]+epsInitxx)
            + (qIMinus[o][YY][i]+epsInityy)*(qIMinus[o][YY][i]+epsInityy)
            + (qIMinus[o][ZZ][i]+epsInitzz)*(qIMinus[o][ZZ][i]+epsInitzz)
            + 2*qIMinus[o][XY][i]*qIMinus[o][XY][i]
            + 2*qIMinus[o][YZ][i]*qIMinus[o][YZ][i]
            + 2*qIMinus[o][XZ][i]*qIMinus[o][XZ][i];
          real alpham = qIMinus[o][DAM][i];
          real xim;
          if (EspIIm > 1e-30){
            xim = EspIm / std::sqrt(EspIIm);
          } else{
            xim = 0.0;
          }

          qStressIMinus[o][XX][i] = (lambda0M*EspIm - alpham*impAndEta[ltsFace].gammaRM*std::sqrt(EspIIm))
                + (2*(mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M)
                    - alpham*impAndEta[ltsFace].gammaRM*xim)
                  *(qIMinus[o][XX][i]+epsInitxx);

          qStressIMinus[o][YY][i] = (lambda0M*EspIm - alpham*impAndEta[ltsFace].gammaRM*std::sqrt(EspIIm))
                + (2*(mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M)
                    - alpham*impAndEta[ltsFace].gammaRM*xim)
                  *(qIMinus[o][YY][i]+epsInityy);

          qStressIMinus[o][ZZ][i] = (lambda0M*EspIm - alpham*impAndEta[ltsFace].gammaRM*std::sqrt(EspIIm))
                + (2*(mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M)
                    - alpham*impAndEta[ltsFace].gammaRM*xim)
                  *(qIMinus[o][ZZ][i]+epsInitzz);

          qStressIMinus[o][XY][i] = 0
                + (2*(mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M)
                    - alpham*impAndEta[ltsFace].gammaRM*xim)
                  *qIMinus[o][XY][i];

          qStressIMinus[o][YZ][i] = 0
                + (2*(mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M)
                    - alpham*impAndEta[ltsFace].gammaRM*xim)
                  *qIMinus[o][YZ][i];

          qStressIMinus[o][XZ][i] = 0
                + (2*(mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M)
                    - alpham*impAndEta[ltsFace].gammaRM*xim)
                  *qIMinus[o][XZ][i];

          qStressIPlus[o][U][i] = qIPlus[o][U][i];
          qStressIPlus[o][V][i] = qIPlus[o][V][i];
          qStressIPlus[o][W][i] = qIPlus[o][W][i];
          qStressIPlus[o][DAM][i] = qIPlus[o][DAM][i];

          qStressIMinus[o][U][i] = qIMinus[o][U][i];
          qStressIMinus[o][V][i] = qIMinus[o][V][i];
          qStressIMinus[o][W][i] = qIMinus[o][W][i];
          qStressIMinus[o][DAM][i] = qIMinus[o][DAM][i];
        }
        // std::cout << qInterpolatedPlus[ltsFace][o][3*seissol::dr::misc::numPaddedPoints+0]
        // << ", " << qStressInterpolatedPlus[o][3*seissol::dr::misc::numPaddedPoints+0]
        // << ", " << qStressInterpolatedPlus[o][3*seissol::dr::misc::numPaddedPoints+0]/qInterpolatedPlus[ltsFace][o][3*seissol::dr::misc::numPaddedPoints+0]
        // << std::endl;
      } // time integration loop

      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                qStressInterpolatedPlus,
                                                qStressInterpolatedMinus);
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
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qStressInterpolatedPlus,
                                                   qStressInterpolatedMinus,
                                                   timeWeights);
      LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
      SCOREP_USER_REGION_END(myRegionHandle)

      common::computeFrictionEnergy(energyData[ltsFace],
                                    qStressInterpolatedPlus,
                                    qStressInterpolatedMinus,
                                    impAndEta[ltsFace],
                                    timeWeights,
                                    spaceWeights,
                                    godunovData[ltsFace]);
    }
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H
