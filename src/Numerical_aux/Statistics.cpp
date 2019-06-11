/**
 * @author Lukas Krenz (lukas.krenz AT tum.de)
 * @author Carsten Uphoff (c.uphoff AT tum.de)
 **/

#include "Statistics.h"

#include <algorithm>
#include <cmath>
#include "Parallel/MPI.h"

seissol::statistics::Summary::Summary(double value)
  : mean(value), std(0.0), min(value), median(value), max(value)
{
}

seissol::statistics::Summary::Summary(std::vector<double> const& values)
{
  std::vector<double> sortedValues(values);
	std::sort(sortedValues.begin(), sortedValues.end());

  auto N = sortedValues.size();

  median = -1;
  if (N % 2 == 1) {
    median = sortedValues[N / 2];
  } else {
    // Median not uniq. defined, take mean of two candidates.
    median = 0.5 *
      (sortedValues[N / 2 - 1] +
       sortedValues[N / 2]);
  }
  min = sortedValues[0];
  max = sortedValues[N - 1];
  
  mean = 0.0;
  auto meanOfSquares = 0.0;
  for (const auto num : sortedValues) {
    mean += num;
    meanOfSquares += num * num;
  }

  mean /= N;
  meanOfSquares /= N;
  
  // Note that this computation is numerically unstable!
  const auto variance = meanOfSquares - mean * mean;
  std = std::sqrt(variance);
}

auto seissol::statistics::parallelSummary(double value) -> Summary {
#ifdef USE_MPI
	const int rank = seissol::MPI::mpi.rank();
	auto collect = std::vector<double>(seissol::MPI::mpi.size());

	MPI_Gather(&value, 1, MPI_DOUBLE,
		   collect.data(), 1, MPI_DOUBLE,
		   0, seissol::MPI::mpi.comm());

	if (rank == 0) {
    return Summary(collect);
  }
  return Summary();
#else
  return Summary(value);
#endif
}
