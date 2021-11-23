#ifndef SEISSOL_DR_OUTOUT_NO_OUTPUT_HPP
#define SEISSOL_DR_OUTOUT_NO_OUTPUT_HPP

#include <Solver/Interoperability.h>
#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class NoFault : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& e_interoperability) override {
    Base::tiePointers(layerData, dynRup, e_interoperability);
  }

  void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTOUT_NO_OUTPUT_HPP
