#ifndef SEISSOL_GPU_LINEARSLIPWEAKENING_H
#define SEISSOL_GPU_LINEARSLIPWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "utils/logger.h"


namespace seissol::dr::friction_law::gpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
class LinearSlipWeakeningBase : public GpuBaseFrictionLaw {
  public:
  LinearSlipWeakeningBase(dr::DRParameters& drParameters) : GpuBaseFrictionLaw(drParameters){};

  /**
   * critical velocity at which slip rate is considered as being zero for instaneous healing
   */
  static constexpr real u0 = 10e-14;

  real (*dC)[misc::numPaddedPoints];
  real (*muS)[misc::numPaddedPoints];
  real (*muD)[misc::numPaddedPoints];
  real (*cohesion)[misc::numPaddedPoints];

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             TractionResults& tractionResults,
                             std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::numPaddedPoints>& strengthBuffer,
                             unsigned int ltsFace,
                             unsigned int timeIndex) override;

  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture* dynRup,
                                      real fullUpdateTime) override;

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  void calcSlipRateAndTraction(FaultStresses& faultStresses,
                               TractionResults& tractionResults,
                               std::array<real, misc::numPaddedPoints>& strength,
                               unsigned int timeIndex,
                               unsigned int ltsFace);

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                            unsigned int ltsFace);

  /**
   * Instantaneous healing option:
   * Reset Mu and Slip, if slipRateMagnitude drops below threshold
   * This function is currently not used, as we miss an appropriate benchmark.
   */
  void instantaneousHealing(unsigned int ltsFace);

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput(unsigned int ltsFace) override;
};


class LinearSlipWeakeningLaw : public LinearSlipWeakeningBase {
  public:
  LinearSlipWeakeningLaw(dr::DRParameters& drParameters) : LinearSlipWeakeningBase(drParameters){};

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
               unsigned int ltsFace) override{};
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                unsigned int ltsFace) override{};

  void calcStrengthHook(FaultStresses& faultStresses,
                        std::array<real, misc::numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace) override;

  void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace) override;
};

class LinearSlipWeakeningLawForcedRuptureTime : public LinearSlipWeakeningLaw {
  public:
  LinearSlipWeakeningLawForcedRuptureTime(dr::DRParameters& drParameters)
      : LinearSlipWeakeningLaw(drParameters){};

  real (*forcedRuptureTime)[misc::numPaddedPoints];
  real* tn;

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
               unsigned int ltsFace) override;
  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture* dynRup,
                                      real fullUpdateTime) override;

  void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace) override;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_LINEARSLIPWEAKENING_H
