/**
 * @author Lukas Krenz (lukas.krenz AT tum.de)
 * @author Carsten Uphoff (c.uphoff AT tum.de)
 **/

#include "Statistics.h"

#include "Parallel/MPI.h"
#include <algorithm>
#include <cmath>
#include <vector>

seissol::statistics::Summary::Summary(double value)
    : mean(value), std(0.0), min(value), median(value), max(value) {}

seissol::statistics::Summary::Summary(const std::vector<double>& values) {
  std::vector<double> sortedValues(values);
  std::sort(sortedValues.begin(), sortedValues.end());

  auto n = sortedValues.size();

  median = -1;
  if (n % 2 == 1) {
    median = sortedValues[n / 2];
  } else {
    // Median not uniq. defined, take mean of two candidates.
    median = 0.5 * (sortedValues[n / 2 - 1] + sortedValues[n / 2]);
  }
  min = sortedValues[0];
  max = sortedValues[n - 1];

  mean = 0.0;
  auto meanOfSquares = 0.0;
  for (const auto num : sortedValues) {
    mean += num;
    meanOfSquares += num * num;
  }

  mean /= n;
  meanOfSquares /= n;

  // Note that this computation is numerically unstable!
  const auto variance = meanOfSquares - mean * mean;
  std = std::sqrt(variance);
}

auto seissol::statistics::parallelSummary(double value) -> Summary {
#ifdef USE_MPI
  auto collect = seissol::MPI::mpi.collect(value);
  const int rank = seissol::MPI::mpi.rank();
  if (rank == 0) {
    return Summary(collect);
  }
  return Summary();
#else
  return Summary(value);
#endif
}
