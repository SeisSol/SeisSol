#ifndef SEISSOL_NOFAULT_H
#define SEISSOL_NOFAULT_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * No friction computation input stress equals output
 */
template <typename Config>
class NoFault : public BaseFrictionLaw<Config, NoFault<Config>> {
  public:
  using RealT = typename Config::RealT;
  using BaseFrictionLaw<Config, NoFault<Config>>::BaseFrictionLaw;

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime) {}

  void updateFrictionAndSlip(
      FaultStresses<Config> const& faultStresses,
      TractionResults<Config>& tractionResults,
      std::array<RealT, misc::numPaddedPoints<Config>> const& stateVariableBuffer,
      std::array<RealT, misc::numPaddedPoints<Config>> const& strengthBuffer,
      unsigned ltsFace,
      unsigned timeIndex);

  void preHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
               unsigned ltsFace){};
  void postHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
                unsigned ltsFace){};
  void saveDynamicStressOutput(unsigned int ltsFace){};
};
} // namespace seissol::dr::friction_law
#endif // SEISSOL_NOFAULT_H
