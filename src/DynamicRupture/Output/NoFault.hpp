#ifndef SEISSOL_DR_OUTPUT_NO_FAULT_HPP
#define SEISSOL_DR_OUTPUT_NO_FAULT_HPP

#include "DynamicRupture/Output/Base.hpp"
#include <Solver/Interoperability.h>

namespace seissol::dr::output {
class NoFault : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* drDescr,
                   seissol::Interoperability& eInteroperability) override {
    Base::tiePointers(layerData, drDescr, eInteroperability);
  }

  void postCompute(seissol::initializers::DynamicRupture& drDescr) override {
    // do nothing
  }

  real computeLocalStrength() override { return 0.0; }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_NO_FAULT_HPP
