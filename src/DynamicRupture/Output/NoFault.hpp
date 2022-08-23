#ifndef SEISSOL_DR_OUTPUT_NO_FAULT_HPP
#define SEISSOL_DR_OUTPUT_NO_FAULT_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class NoFault : public ReceiverOutput {
  protected:
  real computeLocalStrength(LocalInfo& local) override { return 0.0; }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_NO_FAULT_HPP
