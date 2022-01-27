#ifndef SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
#define SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class ImposedSlipRates : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& eInteroperability) override {
    Base::tiePointers(layerData, dynRup, eInteroperability);
  }

  void postCompute(seissol::initializers::DynamicRupture& dynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
