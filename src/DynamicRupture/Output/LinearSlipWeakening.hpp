#ifndef SEISSOL_DR_OUTPUT_LSW_HPP
#define SEISSOL_DR_OUTPUT_LSW_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class LinearSlipWeakening : public ReceiverBasedOutput {
  protected:
  real computeLocalStrength() override {
    using DrLtsDescrT = seissol::initializers::LTS_LinearSlipWeakening;
    auto* cohesions = local.layer->var(static_cast<DrLtsDescrT*>(drDescr)->cohesion);
    auto cohesion = cohesions[local.ltsId][local.nearestGpIndex];

    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
               std::min(effectiveNormalStress, static_cast<real>(0.0)) -
           cohesion;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_LSW_HPP
