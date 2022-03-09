#ifndef SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
#define SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class ImposedSlipRates : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* drDescr,
                   seissol::Interoperability& eInteroperability) override {
    Base::tiePointers(layerData, drDescr, eInteroperability);
  }

  void postCompute(seissol::initializers::DynamicRupture& drDescr) override {
    // do nothing
  }
  real computeLocalStrength() override {
    return -1.0 * local.mu * std::min(local.p + local.p0 - local.pf, static_cast<real>(0.0));
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
