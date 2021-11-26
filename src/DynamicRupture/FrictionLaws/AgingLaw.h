#ifndef SEISSOL_AGINGLAW_H
#define SEISSOL_AGINGLAW_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class AgingLaw : public RateAndStateBase<AgingLaw> {
  public:
  using RateAndStateBase<AgingLaw>::RateAndStateBase;
  using RateAndStateBase<AgingLaw>::copyLtsTreeToLocal;

  real calcStateVariableHook(real SV0, real tmp, real time_inc, real rs_sl0);
  void preCalcTime(){};
  void setInitialValues(std::array<real, numPaddedPoints>& localStateVariable,
                        unsigned int ltsFace) {}
  void calcInitialSlipRate(std::array<real, numPaddedPoints>& totalShearStressYZ,
                           FaultStresses& faultStresses,
                           std::array<real, numPaddedPoints>& stateVarZero,
                           std::array<real, numPaddedPoints>& localStateVariable,
                           std::array<real, numPaddedPoints>& temporarySlipRate,
                           unsigned int timeIndex,
                           unsigned int ltsFace) {}
  void hookSetInitialP_f(std::array<real, numPaddedPoints>& P_f, unsigned int ltsFace) {}

  void hookCalcP_f(std::array<real, numPaddedPoints>& P_f,
                   FaultStresses& faultStresses,
                   bool saveTmpInTP,
                   unsigned int timeIndex,
                   unsigned int ltsFace) {}
  void updateStateVariableIterative(bool& has_converged,
                                    std::array<real, numPaddedPoints>& stateVarZero,
                                    std::array<real, numPaddedPoints>& SR_tmp,
                                    std::array<real, numPaddedPoints>& LocSV,
                                    std::array<real, numPaddedPoints>& P_f,
                                    std::array<real, numPaddedPoints>& normalStress,
                                    std::array<real, numPaddedPoints>& TotalShearStressYZ,
                                    std::array<real, numPaddedPoints>& SRtest,
                                    FaultStresses& faultStresses,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace) {}

  void executeIfNotConverged(std::array<real, numPaddedPoints>& LocSV, unsigned ltsFace) {}
  void calcSlipRateAndTraction(std::array<real, numPaddedPoints>& stateVarZero,
                               std::array<real, numPaddedPoints>& SR_tmp,
                               std::array<real, numPaddedPoints>& localStateVariable,
                               std::array<real, numPaddedPoints>& normalStress,
                               std::array<real, numPaddedPoints>& TotalShearStressYZ,
                               FaultStresses& faultStresses,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {}
  void resampleStateVar(std::array<real, numPaddedPoints>& deltaStateVar, unsigned int ltsFace) {}
  void saveDynamicStressOutput(unsigned int face) {}
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_AGINGLAW_H
