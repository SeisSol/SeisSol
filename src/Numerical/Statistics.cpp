// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Lukas Krenz
// SPDX-FileContributor: Carsten Uphoff

#include "Statistics.h"

#include "Parallel/MPI.h"
#include <algorithm>
#include <cmath>
#include <vector>

seissol::statistics::Summary::Summary(double value)
    : mean(value), std(0.0), min(value), median(value), max(value) {}

seissol::statistics::Summary::Summary(const std::vector<double>& values) : median(-1) {
  std::vector<double> sortedValues(values);
  std::sort(sortedValues.begin(), sortedValues.end());

  auto n = sortedValues.size();

  if (n % 2 == 1) {
    median = sortedValues[n / 2];
  } else {
    // Median not uniq. defined, take mean of two candidates.
    median = 0.5 * (sortedValues[n / 2 - 1] + sortedValues[n / 2]);
  }
  min = sortedValues[0];
  max = sortedValues[n - 1];

  sum = 0.0;
  for (const auto num : sortedValues) {
    sum += num;
  }
  mean = sum / n;

  auto variance = 0.0;
  for (const auto num : sortedValues) {
    const auto diff = mean - num;
    variance += diff * diff;
  }

  variance /= n;

  std = std::sqrt(variance);
}

auto seissol::statistics::parallelSummary(double value) -> Summary {
#ifdef USE_MPI
  auto collect = seissol::MPI::mpi.collect(value);
  const int rank = seissol::MPI::mpi.rank();
  if (rank == 0) {
    return Summary(collect);
  }
  return {};
#else
  return Summary(value);
#endif
}
