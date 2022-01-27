#ifndef SEISSOL_DR_OUTPUT_NO_FAULT_HPP
#define SEISSOL_DR_OUTPUT_NO_FAULT_HPP

#include "DynamicRupture/Output/Base.hpp"
#include <Solver/Interoperability.h>

namespace seissol::dr::output {
class NoFault : public Base {
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

#endif // SEISSOL_DR_OUTPUT_NO_FAULT_HPP
