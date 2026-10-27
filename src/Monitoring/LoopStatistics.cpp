// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "LoopStatistics.h"

#include "Common/Json.h"
#include "Monitoring/RunSummary.h"
#include "Monitoring/Stopwatch.h"
#include "Numerical/Statistics.h"
#include "Unit.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iterator>
#include <mpi.h>
#include <string>
#include <time.h>
#include <utils/logger.h>
#include <vector>

namespace seissol {

LoopStatistics::Region::Region(const std::string& name, bool includeInSummary)
    : name(name), includeInSummary(includeInSummary) {}

void LoopStatistics::addRegion(const std::string& name, bool includeInSummary) {
  regions_.emplace_back(name, includeInSummary);
}

unsigned LoopStatistics::getRegion(const std::string& name) const {
  const auto first = regions_.cbegin();
  const auto it =
      std::find_if(first, regions_.cend(), [&name](const auto& elem) { return elem.name == name; });
  assert(it != regions_.end());
  return std::distance(first, it);
}

void LoopStatistics::begin(unsigned region) {
  (void)clock_gettime(CLOCK_MONOTONIC, &regions_[region].begin);
}

void LoopStatistics::end(unsigned region, unsigned numIterations, unsigned subRegion) {
  timespec endTime{};
  (void)clock_gettime(CLOCK_MONOTONIC, &endTime);
  addSample(region, numIterations, subRegion, regions_[region].begin, endTime);
}

void LoopStatistics::addSample(
    unsigned region, unsigned numIterations, unsigned /*subRegion*/, timespec begin, timespec end) {
  if (numIterations > 0) {
    auto& vars = regions_[region].variables;
    const auto time = seconds(difftime(begin, end));
    vars.x += numIterations;
    vars.x2 += static_cast<double>(numIterations) * static_cast<double>(numIterations);
    vars.xy += static_cast<double>(numIterations) * time;
    vars.y += time;
    vars.y2 += time * time;
    ++vars.n;
  }
}

void LoopStatistics::reset() {
  for (auto& region : regions_) {
    region.variables = StatisticVariables();
    // (region.begin is not reset)
  }
}

double LoopStatistics::kernelTime() const {
  double time = 0.0;
  for (const auto& region : regions_) {
    // the regions that are not, such as the actor states, overlap with the kernels
    if (region.includeInSummary) {
      time += region.variables.y;
    }
  }
  return time;
}

json::Value LoopStatistics::printSummary(MPI_Comm comm) {
  const auto nRegions = regions_.size();
  constexpr int NumSumComponents = 6;
  auto sums = std::vector<double>(NumSumComponents * nRegions);
  const double totalTimePerRank = kernelTime();

  // Helper functions to access sums
  const auto getNumIters = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 0];
  };
  const auto getNumItersSquared = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 1];
  };
  const auto getTimeTimesNumIters = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 2];
  };
  const auto getTime = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 3];
  };
  const auto getTimeSquared = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 4];
  };
  const auto getNumberOfSamples = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 5];
  };

  for (unsigned region = 0; region < nRegions; ++region) {
    getNumIters(region) = regions_[region].variables.x;
    getNumItersSquared(region) = regions_[region].variables.x2;
    getTimeTimesNumIters(region) = regions_[region].variables.xy;
    getTime(region) = regions_[region].variables.y;
    getTimeSquared(region) = regions_[region].variables.y2;
    getNumberOfSamples(region) = regions_[region].variables.n;
  }

  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  const auto summary = seissol::statistics::parallelSummary(totalTimePerRank);
  logInfo() << "Time spent in compute kernels: mean =" << summary.mean << " std =" << summary.std
            << " min =" << summary.min << " median =" << summary.median << " max =" << summary.max;

  const auto loadImbalance = 1.0 - summary.mean / summary.max;
  logInfo() << "Load imbalance:" << 100.0 * loadImbalance << "%";

  // the time of every region per rank, which the sums below fold together
  std::vector<seissol::statistics::Summary> regionTimes;
  regionTimes.reserve(nRegions);
  for (unsigned region = 0; region < nRegions; ++region) {
    regionTimes.push_back(seissol::statistics::parallelSummary(regions_[region].variables.y));
  }

  MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM, comm);

  auto regressionCoeffs = std::vector<double>(2 * nRegions);
  auto stderror = std::vector<double>(nRegions, 0.0);
  for (unsigned region = 0; region < nRegions; ++region) {
    const double x = getNumIters(region);
    const double x2 = getNumItersSquared(region);
    const double xy = getTimeTimesNumIters(region);
    const double y = getTime(region);
    const double y2 = getTimeSquared(region);
    const double n = getNumberOfSamples(region);

    const double det = n * x2 - x * x;
    const double constant = (x2 * y - x * xy) / det;
    const double slope = (-x * y + n * xy) / det;
    regressionCoeffs[2 * region + 0] = constant;
    regressionCoeffs[2 * region + 1] = slope;

    // for each sample, do:
    // error = time - (constant + slope * numIters)
    // stderror[region] += error * error
    // sum up over all samples over all ranks (i.e. MPI_Reduce)

    // that can be replaced by (including the MPI operation):
    // error * error = constant * constant + slope * slope * numIters * numIters + time * time - 2 *
    // constant * time - 2 * slope * numIters * time + 2 * constant * slope * numIters thus, by
    // summing up error * error over all samples (and using x, x2, xy, y, y2):

    stderror[region] = constant * constant + slope * slope * x2 + y2 - 2 * constant * y -
                       2 * slope * xy + 2 * constant * slope * x;
  }

  if (rank == 0) {
    double totalTime = 0.0;
    logInfo() << "Regression analysis of compute kernels:";
    for (unsigned region = 0; region < nRegions; ++region) {
      if (!regions_[region].includeInSummary) {
        continue;
      }
      const double x = getNumIters(region);
      const double x2 = getNumItersSquared(region);
      const double y = getTime(region);
      const double n = getNumberOfSamples(region);

      const double xm = x / n;
      const double xv = x2 - 2 * x * xm + xm * xm;

      // https://en.wikipedia.org/wiki/Simple_linear_regression#Normality_assumption
      const double se = std::sqrt((stderror[region] / (n - 2)) / xv);

      const auto names = std::vector<std::string>{"constant", "per element"};
      logInfo() << regions_[region].name << "(total time):" << y
                << "s ( =" << UnitTime.formatTime(y).c_str() << ")";
      for (unsigned c = 0; c < 2; ++c) {
        logInfo() << regions_[region].name << "(" << names[c]
                  << "):" << regressionCoeffs[2 * region + c] << "(sample size:" << n
                  << ", standard error:" << se << ")";
      }
      totalTime += y;
    }

    logInfo() << "Total time spent in compute kernels:" << totalTime
              << "s ( =" << UnitTime.formatTime(totalTime).c_str() << ")";
  }

  json::Value result;
  result["total"] = monitoring::toJson(summary);
  result["loadImbalance"] = loadImbalance;
  auto& regions = result["regions"] = json::Value::object();
  for (unsigned region = 0; region < nRegions; ++region) {
    const double x = getNumIters(region);
    const double x2 = getNumItersSquared(region);
    const double n = getNumberOfSamples(region);
    const double xm = x / n;
    const double xv = x2 - 2 * x * xm + xm * xm;
    auto& entry = regions[regions_[region].name];
    entry["includedInTotal"] = regions_[region].includeInSummary;
    entry["time"] = monitoring::toJson(regionTimes[region]);
    entry["iterations"] = x;
    entry["samples"] = n;
    // the regions outside of the total time count no elements, but states
    if (regions_[region].includeInSummary) {
      entry["regression"]["constant"] = regressionCoeffs[2 * region + 0];
      entry["regression"]["perElement"] = regressionCoeffs[2 * region + 1];
      entry["regression"]["standardError"] = std::sqrt((stderror[region] / (n - 2)) / xv);
    }
  }
  return result;
}

} // namespace seissol
