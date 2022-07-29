#ifndef SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP
#define SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  real computeLocalStrength() override {
    using DrLtsDescrT = seissol::initializers::LTS_LinearSlipWeakeningBimaterial;
    auto const* const regularisedStrengths =
        local.layer->var(static_cast<DrLtsDescrT*>(drDescr)->regularisedStrength);
    return regularisedStrengths[local.ltsId][local.nearestGpIndex];
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP
