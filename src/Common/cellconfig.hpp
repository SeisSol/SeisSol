#ifndef SEISSOL_COMMON_CELLCONFIG_HPP_
#define SEISSOL_COMMON_CELLCONFIG_HPP_

namespace seissol {
// TODO(David): once C++20 hits, make CellConfig implement a concept
// dummy template type
template <typename MaterialTP, typename RealTP, int ConvergenceOrderP, bool PlasticityP>
struct CellConfig {
  using MaterialT = MaterialTP;
  using RealT = RealTP;
  static constexpr bool Plasticity = PlasticityP;
  static constexpr int ConvergenceOrder = ConvergenceOrderP;
};
} // namespace seissol

#endif
