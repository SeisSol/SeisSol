// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_STABLESUM_H_
#define SEISSOL_SRC_NUMERICAL_STABLESUM_H_

namespace seissol::numerical {

/**
A Kahan Summation implementation for a more numerically-stable summation.
(beware of aggressive compilers)

Cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm .

How to use:

auto acc = StableAccumulator();
acc += 1e100;
acc += 1e-100;
acc -= 1e100;
std::cout << acc.result() << std::endl;
*/
template <typename RealT>
struct StableAccumulator {
  StableAccumulator() = default;
  StableAccumulator(RealT start) : acc(start), corr(0) {}

  auto operator+(RealT number) -> StableAccumulator<RealT> {
    StableAccumulator<RealT> newacc;
    const auto numberC = number - corr;
    newacc.acc = acc + numberC;
    newacc.corr = (newacc.acc - acc) - numberC;
    return newacc;
  }

  auto operator+=(RealT number) -> StableAccumulator<RealT>& {
    const auto tempnew = *this + number;
    this->acc = tempnew.acc;
    this->corr = tempnew.corr;
    return *this;
  }

  auto operator-(RealT number) -> StableAccumulator<RealT> { return (*this + (-number)); }

  auto operator-=(RealT number) -> StableAccumulator<RealT> { return (*this += (-number)); }

  [[nodiscard]] RealT result() const { return acc; }

  private:
  RealT acc{0};
  RealT corr{0};
};

} // namespace seissol::numerical

#endif // SEISSOL_SRC_NUMERICAL_STABLESUM_H_
