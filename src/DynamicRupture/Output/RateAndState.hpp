#ifndef SEISSOL_DR_OUTPUT_RS_HPP
#define SEISSOL_DR_OUTPUT_RS_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class RateAndState : public ReceiverOutput {
  protected:
  real computeLocalStrength(LocalInfo& local) override {
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
           std::min(effectiveNormalStress, static_cast<real>(0.0));
  }

  real computeStateVariable(LocalInfo& local) override {
    const auto* descr = reinterpret_cast<seissol::initializer::LTSRateAndState*>(drDescr);
    assert((descr != nullptr) && "dr descr. must be a subtype of LTS_RateAndState");
    return (local.layer->var(descr->stateVariable))[local.ltsId][local.nearestGpIndex];
  }
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_RS_HPP
