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
      auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus[ltsFace]));
      auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus[ltsFace]));
      auto* qStressIPlus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedPlus));
      auto* qStressIMinus = (reinterpret_cast<QStressInterpolatedShapeT>(qStressInterpolatedMinus));

      using namespace seissol::dr::misc::quantity_indices;
      unsigned DAM = 9;
      unsigned BRE = 10;

      real epsInitxx = 0.0e-4; // eps_xx0
      real epsInityy = 0.0e-3; // eps_yy0
      real epsInitzz = 0.0e-4; // eps_zz0
      real epsInitxy = 0.0e-3; // eps_xy0
      real epsInityz = -0e-1; // eps_yz0
      real epsInitzx = -0e-1; // eps_zx0
      real lambda0P = impAndEta[ltsFace].lambda0P;
      real mu0P = impAndEta[ltsFace].mu0P;
      real rho0P = impAndEta[ltsFace].rho0P;
      real lambda0M = impAndEta[ltsFace].lambda0M;
      real mu0M = impAndEta[ltsFace].mu0M;
      real rho0M = impAndEta[ltsFace].rho0M;

      //TODO(NONLINEAR) What are these values?
      real aB0 = 7.43e9;
      real aB1 = -12.14e9;
      real aB2 = 18.93e9;
      real aB3 = -5.067e9;

      for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
        for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints;
            i ++) {

          real EspIp = (qIPlus[o][XX][i]+epsInitxx) + (qIPlus[o][YY][i]+epsInityy) + (qIPlus[o][ZZ][i]+epsInitzz);
          real EspIIp = (qIPlus[o][XX][i]+epsInitxx)*(qIPlus[o][XX][i]+epsInitxx)
            + (qIPlus[o][YY][i]+epsInityy)*(qIPlus[o][YY][i]+epsInityy)
            + (qIPlus[o][ZZ][i]+epsInitzz)*(qIPlus[o][ZZ][i]+epsInitzz)
            + 2*(qIPlus[o][XY][i]+epsInitxy)*(qIPlus[o][XY][i]+epsInitxy)
            + 2*(qIPlus[o][YZ][i]+epsInityz)*(qIPlus[o][YZ][i]+epsInityz)
            + 2*(qIPlus[o][XZ][i]+epsInitzx)*(qIPlus[o][XZ][i]+epsInitzx);
          real alphap = qIPlus[o][DAM][i];
          real xip;
          if (EspIIp > 1e-30){
            xip = EspIp / std::sqrt(EspIIp);
          } else{
            xip = 0.0;
          }

          // damage stress impAndEtaGet->gammaRP, mu0P
          real mu_eff = mu0P - alphap*impAndEta[ltsFace].gammaRP*impAndEta[ltsFace].xi0P
              - 0.5*alphap*impAndEta[ltsFace].gammaRP*xip;
          real sxx_sp = lambda0P*EspIp
                        - alphap*impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp)
                        + 2*mu_eff*(qIPlus[o][XX][i]+epsInitxx);
          real syy_sp = lambda0P*EspIp
                        - alphap*impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp)
                        + 2*mu_eff*(qIPlus[o][YY][i]+epsInityy);
          real szz_sp = lambda0P*EspIp
                        - alphap*impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp)
                        + 2*mu_eff*(qIPlus[o][ZZ][i]+epsInitzz);

          real sxy_sp = 2*mu_eff*(qIPlus[o][XY][i]+epsInitxy);
          real syz_sp = 2*mu_eff*(qIPlus[o][YZ][i]+epsInityz);
          real szx_sp = 2*mu_eff*(qIPlus[o][XZ][i]+epsInitzx);

          // breakage stress
          real sxx_bp = (2.0*aB2 + 3.0*xip*aB3)*EspIp
                        + aB1 * std::sqrt(EspIIp)
                        + (2.0*aB0 + aB1*xip - aB3*xip*xip*xip)*(qIPlus[o][XX][i]+epsInitxx);
          real syy_bp = (2.0*aB2 + 3.0*xip*aB3)*EspIp
                        + aB1 * std::sqrt(EspIIp)
                        + (2.0*aB0 + aB1*xip - aB3*xip*xip*xip)*(qIPlus[o][YY][i]+epsInityy);
          real szz_bp = (2.0*aB2 + 3.0*xip*aB3)*EspIp
                        + aB1 * std::sqrt(EspIIp)
                        + (2.0*aB0 + aB1*xip - aB3*xip*xip*xip)*(qIPlus[o][ZZ][i]+epsInitzz);

          real sxy_bp = (2.0*aB0 + aB1*xip - aB3*xip*xip*xip)*(qIPlus[o][XY][i]+epsInitxy);
          real syz_bp = (2.0*aB0 + aB1*xip - aB3*xip*xip*xip)*(qIPlus[o][YZ][i]+epsInityz);
          real szx_bp = (2.0*aB0 + aB1*xip - aB3*xip*xip*xip)*(qIPlus[o][XZ][i]+epsInitzx);

          qStressIPlus[o][XX][i] =
            (1-qIPlus[o][BRE][i]) * sxx_sp
            + qIPlus[o][BRE][i] * sxx_bp;

          qStressIPlus[o][YY][i] = 
            (1-qIPlus[o][BRE][i]) * syy_sp
            + qIPlus[o][BRE][i] * syy_bp;

          qStressIPlus[o][ZZ][i] = 
            (1-qIPlus[o][BRE][i]) * szz_sp
            + qIPlus[o][BRE][i] * szz_bp;

          qStressIPlus[o][XY][i] =
            (1-qIPlus[o][BRE][i]) * sxy_sp
            + qIPlus[o][BRE][i] * sxy_bp;

          qStressIPlus[o][YZ][i] =
            (1-qIPlus[o][BRE][i]) * syz_sp
            + qIPlus[o][BRE][i] * syz_bp;

          qStressIPlus[o][XZ][i] =
            (1-qIPlus[o][BRE][i]) * szx_sp
            + qIPlus[o][BRE][i] * szx_bp;

          real EspIm = (qIMinus[o][XX][i]+epsInitxx) + (qIMinus[o][YY][i]+epsInityy) + (qIMinus[o][ZZ][i]+epsInitzz);
          real EspIIm = (qIMinus[o][XX][i]+epsInitxx)*(qIMinus[o][XX][i]+epsInitxx)
            + (qIMinus[o][YY][i]+epsInityy)*(qIMinus[o][YY][i]+epsInityy)
            + (qIMinus[o][ZZ][i]+epsInitzz)*(qIMinus[o][ZZ][i]+epsInitzz)
            + 2*(qIMinus[o][XY][i]+epsInitxy)*(qIMinus[o][XY][i]+epsInitxy)
            + 2*(qIMinus[o][YZ][i]+epsInityz)*(qIMinus[o][YZ][i]+epsInityz)
            + 2*(qIMinus[o][XZ][i]+epsInitzx)*(qIMinus[o][XZ][i]+epsInitzx);
          real alpham = qIMinus[o][DAM][i];
          real xim;
          if (EspIIm > 1e-30){
            xim = EspIm / std::sqrt(EspIIm);
          } else{
            xim = 0.0;
          }

          // damage stress minus
          mu_eff = mu0M - alpham*impAndEta[ltsFace].gammaRM*impAndEta[ltsFace].xi0M
              - 0.5*alpham*impAndEta[ltsFace].gammaRM*xim;
          real sxx_sm = lambda0M*EspIm
                        - alpham*impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm)
                        + 2*mu_eff*(qIMinus[o][XX][i]+epsInitxx);
          real syy_sm = lambda0M*EspIm
                        - alpham*impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm)
                        + 2*mu_eff*(qIMinus[o][YY][i]+epsInityy);
          real szz_sm = lambda0M*EspIm
                        - alpham*impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm)
                        + 2*mu_eff*(qIMinus[o][ZZ][i]+epsInitzz);

          real sxy_sm = 2*mu_eff*(qIMinus[o][XY][i]+epsInitxy);
          real syz_sm = 2*mu_eff*(qIMinus[o][YZ][i]+epsInityz);
          real szx_sm = 2*mu_eff*(qIMinus[o][XZ][i]+epsInitzx);

          // breakage stress
          real sxx_bm = (2.0*aB2 + 3.0*xim*aB3)*EspIm
                        + aB1 * std::sqrt(EspIIm)
                        + (2.0*aB0 + aB1*xim - aB3*xim*xim*xim)*(qIMinus[o][XX][i]+epsInitxx);
          real syy_bm = (2.0*aB2 + 3.0*xim*aB3)*EspIm
                        + aB1 * std::sqrt(EspIIm)
                        + (2.0*aB0 + aB1*xim - aB3*xim*xim*xim)*(qIMinus[o][YY][i]+epsInityy);
          real szz_bm = (2.0*aB2 + 3.0*xim*aB3)*EspIm
                        + aB1 * std::sqrt(EspIIm)
                        + (2.0*aB0 + aB1*xim - aB3*xim*xim*xim)*(qIMinus[o][ZZ][i]+epsInitzz);

          real sxy_bm = (2.0*aB0 + aB1*xim - aB3*xim*xim*xim)*(qIMinus[o][XY][i]+epsInitxy);
          real syz_bm = (2.0*aB0 + aB1*xim - aB3*xim*xim*xim)*(qIMinus[o][YZ][i]+epsInityz);
          real szx_bm = (2.0*aB0 + aB1*xim - aB3*xim*xim*xim)*(qIMinus[o][XZ][i]+epsInitzx);


          qStressIMinus[o][XX][i] = 
            (1-qIMinus[o][BRE][i]) * sxx_sm
            + qIMinus[o][BRE][i] * sxx_bm;

          qStressIMinus[o][YY][i] =
            (1-qIMinus[o][BRE][i]) * syy_sm
            + qIMinus[o][BRE][i] * syy_bm;

          qStressIMinus[o][ZZ][i] =
            (1-qIMinus[o][BRE][i]) * szz_sm
            + qIMinus[o][BRE][i] * szz_bm;

          qStressIMinus[o][XY][i] = 
            (1-qIMinus[o][BRE][i]) * sxy_sm
            + qIMinus[o][BRE][i] * sxy_bm;

          qStressIMinus[o][YZ][i] =
            (1-qIMinus[o][BRE][i]) * syz_sm
            + qIMinus[o][BRE][i] * syz_bm;

          qStressIMinus[o][XZ][i] = 
            (1-qIMinus[o][BRE][i]) * szx_sm
            + qIMinus[o][BRE][i] * szx_bm;

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

      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                qStressInterpolatedPlus,
                                                qStressInterpolatedMinus,
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace]);
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
