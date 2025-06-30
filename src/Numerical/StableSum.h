// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_STABLESUM_H_
#define SEISSOL_SRC_NUMERICAL_STABLESUM_H_

namespace seissol::numerical {

template <typename RealT>
struct StableAccumulator {
  StableAccumulator() = default;
  StableAccumulator(RealT start) : acc(start), corr(0) {}

  auto operator+(RealT number) -> StableAccumulator<RealT> {
    StableAccumulator<RealT> newacc;
    const auto numberC = number - corr;
    newacc.acc = acc + number;
    newacc.corr = (newacc.acc - number) - numberC;
    return newacc;
  }

  auto operator+=(RealT number) -> StableAccumulator<RealT>& {
    const auto tempnew = *this + number;
    this->acc = tempnew.acc;
    this->corr = tempnew.corr;
    return *this;
  }

  RealT result() const { return acc; }

  private:
  RealT acc{0};
  RealT corr{0};
};

} // namespace seissol::numerical

#endif // SEISSOL_SRC_NUMERICAL_STABLESUM_H_
