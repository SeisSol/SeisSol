// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 **/

#ifndef SEISSOL_SRC_NUMERICAL_STATISTICS_H_
#define SEISSOL_SRC_NUMERICAL_STATISTICS_H_

#include <vector>

namespace seissol {
namespace statistics {
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
} // namespace statistics
} // namespace seissol

#endif // SEISSOL_SRC_NUMERICAL_STATISTICS_H_
