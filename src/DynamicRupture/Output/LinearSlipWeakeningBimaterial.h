#ifndef SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP
#define SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  real computeLocalStrength(LocalInfo& local) override {
    using DrLtsDescrType = seissol::initializer::LTSLinearSlipWeakeningBimaterial;
    const auto* const regularisedStrengths =
        getCellData(local, static_cast<DrLtsDescrType*>(drDescr)->regularisedStrength);
    return regularisedStrengths[local.nearestGpIndex];
  }

  std::vector<std::size_t> getOutputVariables() const override {
    using DrLtsDescrType = seissol::initializer::LTSLinearSlipWeakeningBimaterial;
    auto baseVector = LinearSlipWeakening::getOutputVariables();
    baseVector.push_back(static_cast<DrLtsDescrType*>(drDescr)->regularisedStrength.index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP
