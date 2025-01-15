// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
//
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_NUMERICAL_STATISTICS_H_
#define SEISSOL_SRC_NUMERICAL_STATISTICS_H_

#include <vector>

namespace seissol::statistics {
struct Summary {
  Summary(double value = 0.0);
  Summary(const std::vector<double>& values);

  double mean;
  double std;
  double min;
  double median;
  double max;
};

auto parallelSummary(double value) -> Summary;
} // namespace seissol::statistics

#endif // SEISSOL_SRC_NUMERICAL_STATISTICS_H_
