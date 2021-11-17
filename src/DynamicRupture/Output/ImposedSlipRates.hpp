#ifndef SEISSOL_DR_OUTOUT_IMPOSED_RS_HPP
#define SEISSOL_DR_OUTOUT_IMPOSED_RS_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class ImposedSlipRates : public Base {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    Base::tiePointers(layerData, dynRup, e_interoperability);
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTOUT_IMPOSED_RS_HPP
