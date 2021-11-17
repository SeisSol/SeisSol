#ifndef SEISSOL_DR_OUTOUT_RS_TP_HPP
#define SEISSOL_DR_OUTOUT_RS_TP_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class RateAndStateThermalPressurisation : public RateAndState {
  public:
  using RateAndState::postCompute;
  using RateAndState::RateAndState;
  using RateAndState::tiePointers;
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTOUT_RS_TP_HPP
