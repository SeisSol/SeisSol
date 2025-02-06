// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "LoopStatistics.h"
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
#ifdef USE_NETCDF
#include <fstream>
#include <netcdf.h>
#include <ostream>
#include <sstream>
#ifdef USE_MPI
#include <netcdf_par.h>
#endif // USE_MPI
#endif // USE_NETCDF

#include "Monitoring/Stopwatch.h"
#include "Numerical/Statistics.h"

#ifdef USE_NETCDF
namespace {

void check_err(const int stat, const int line, const char* file) {
  if (stat != NC_NOERR) {
    logError() << "line" << line << "of" << file << ":" << nc_strerror(stat) << std::endl;
  }
}

template <typename T>
nc_type type2nc() {
  if constexpr (std::is_signed_v<T>) {
    static_assert(std::is_integral_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
                  "type2nc supports 32 or 64 bit integral types only");
    if constexpr (sizeof(T) == 4) {
      return NC_INT;
    } else {
      return NC_INT64;
    }
  } else {
    if constexpr (sizeof(T) == 4) {
      return NC_UINT;
    } else {
      return NC_UINT64;
    }
  }
}

} // namespace
#endif

namespace seissol {

void LoopStatistics::enableSampleOutput(bool enabled) { outputSamples = enabled; }

LoopStatistics::Region::Region(const std::string& name, bool includeInSummary)
    : name(name), includeInSummary(includeInSummary) {}

void LoopStatistics::addRegion(const std::string& name, bool includeInSummary) {
  regions.emplace_back(name, includeInSummary);
}

unsigned LoopStatistics::getRegion(const std::string& name) const {
  auto first = regions.cbegin();
  auto it =
      std::find_if(first, regions.cend(), [&name](const auto& elem) { return elem.name == name; });
  assert(it != regions.end());
  return std::distance(first, it);
}

void LoopStatistics::begin(unsigned region) {
  clock_gettime(CLOCK_MONOTONIC, &regions[region].begin);
}

void LoopStatistics::end(unsigned region, unsigned numIterations, unsigned subRegion) {
  timespec endTime{};
  clock_gettime(CLOCK_MONOTONIC, &endTime);
  addSample(region, numIterations, subRegion, regions[region].begin, endTime);
}

void LoopStatistics::addSample(
    unsigned region, unsigned numIterations, unsigned subRegion, timespec begin, timespec end) {
  if (outputSamples) {
    Sample sample{};
    sample.begin = begin;
    sample.end = end;
    sample.numIters = numIterations;
    sample.subRegion = subRegion;
    regions[region].times.emplace_back(sample);
  }
  if (numIterations > 0) {
    auto& vars = regions[region].variables;
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
  for (auto& region : regions) {
    region.times.resize(0);
    region.variables = StatisticVariables();
    // (region.begin is not reset)
  }
}

void LoopStatistics::printSummary(MPI_Comm comm) {
  const auto nRegions = regions.size();
  constexpr int NumSumComponents = 6;
  auto sums = std::vector<double>(NumSumComponents * nRegions);
  double totalTimePerRank = 0.0;

  // Helper functions to access sums
  auto getNumIters = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 0];
  };
  auto getNumItersSquared = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 1];
  };
  auto getTimeTimesNumIters = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 2];
  };
  auto getTime = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 3];
  };
  auto getTimeSquared = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 4];
  };
  auto getNumberOfSamples = [&sums](std::size_t region) -> double& {
    return sums[NumSumComponents * region + 5];
  };

  for (unsigned region = 0; region < nRegions; ++region) {
    getNumIters(region) = regions[region].variables.x;
    getNumItersSquared(region) = regions[region].variables.x2;
    getTimeTimesNumIters(region) = regions[region].variables.xy;
    getTime(region) = regions[region].variables.y;
    getTimeSquared(region) = regions[region].variables.y2;
    getNumberOfSamples(region) = regions[region].variables.n;

    // Make sure that events that lead to duplicate accounting are ignored
    if (regions[region].includeInSummary) {
      totalTimePerRank += regions[region].variables.y;
    }
  }

  int rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(comm, &rank);
#else
  rank = 0;
#endif

  const auto summary = seissol::statistics::parallelSummary(totalTimePerRank);
  logInfo() << "Time spent in compute kernels: mean =" << summary.mean << " std =" << summary.std
            << " min =" << summary.min << " median =" << summary.median << " max =" << summary.max;

  const auto loadImbalance = 1.0 - summary.mean / summary.max;
  logInfo() << "Load imbalance:" << 100.0 * loadImbalance << "%";

#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM, comm);
#endif

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
      if (!regions[region].includeInSummary) {
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

      const char* names[] = {"constant", "per element"};
      logInfo() << regions[region].name << "(total time):" << y
                << "s ( =" << UnitTime.formatTime(y).c_str() << ")";
      for (unsigned c = 0; c < 2; ++c) {
        logInfo() << regions[region].name << "(" << names[c]
                  << "):" << regressionCoeffs[2 * region + c] << "(sample size:" << n
                  << ", standard error:" << se << ")";
      }
      totalTime += y;
    }

    logInfo() << "Total time spent in compute kernels:" << totalTime
              << "s ( =" << UnitTime.formatTime(totalTime).c_str() << ")";
  }
}

void LoopStatistics::writeSamples(const std::string& outputPrefix,
                                  bool isLoopStatisticsNetcdfOutputOn) {
  if (isLoopStatisticsNetcdfOutputOn) {
    const auto loopStatFile = outputPrefix + "-loopStat-";
    const auto rank = MPI::mpi.rank();
#if defined(USE_NETCDF) && defined(USE_MPI)
    logInfo() << "Starting to write loop statistics samples to disk.";
    const unsigned nRegions = regions.size();
    for (unsigned region = 0; region < nRegions; ++region) {
      const std::ofstream file;
      std::stringstream ss;
      ss << loopStatFile << regions[region].name << ".nc";
      const std::string fileName = ss.str();

      long nSamples = regions[region].times.size();
      long sampleOffset = 0;
      MPI_Scan(&nSamples, &sampleOffset, 1, MPI_LONG, MPI_SUM, MPI::mpi.comm());

      int ncid = 0;
      int stat = 0;
      stat = nc_create_par(fileName.c_str(),
                           NC_MPIIO | NC_CLOBBER | NC_NETCDF4,
                           MPI::mpi.comm(),
                           MPI_INFO_NULL,
                           &ncid);
      check_err(stat, __LINE__, __FILE__);

      int sampledim = 0;
      int rankdim = 0;
      int timespectyp = 0;
      int sampletyp = 0;
      int offsetid = 0;
      int sampleid = 0;

      stat = nc_def_dim(ncid, "rank", 1 + MPI::mpi.size(), &rankdim);
      check_err(stat, __LINE__, __FILE__);
      stat = nc_def_dim(ncid, "sample", NC_UNLIMITED, &sampledim);
      check_err(stat, __LINE__, __FILE__);

      stat = nc_def_compound(ncid, sizeof(timespec), "timespec", &timespectyp);
      check_err(stat, __LINE__, __FILE__);
      {
        stat = nc_insert_compound(ncid,
                                  timespectyp,
                                  "sec",
                                  NC_COMPOUND_OFFSET(timespec, tv_sec),
                                  type2nc<decltype(timespec::tv_sec)>());
        check_err(stat, __LINE__, __FILE__);
        stat = nc_insert_compound(ncid,
                                  timespectyp,
                                  "nsec",
                                  NC_COMPOUND_OFFSET(timespec, tv_nsec),
                                  type2nc<decltype(timespec::tv_nsec)>());
        check_err(stat, __LINE__, __FILE__);
      }

      stat = nc_def_compound(ncid, sizeof(Sample), "Sample", &sampletyp);
      check_err(stat, __LINE__, __FILE__);
      {
        stat = nc_insert_compound(
            ncid, sampletyp, "begin", NC_COMPOUND_OFFSET(Sample, begin), timespectyp);
        check_err(stat, __LINE__, __FILE__);
        stat = nc_insert_compound(
            ncid, sampletyp, "end", NC_COMPOUND_OFFSET(Sample, end), timespectyp);
        check_err(stat, __LINE__, __FILE__);
        stat = nc_insert_compound(
            ncid, sampletyp, "loopLength", NC_COMPOUND_OFFSET(Sample, numIters), NC_UINT);
        check_err(stat, __LINE__, __FILE__);
        stat = nc_insert_compound(
            ncid, sampletyp, "subRegion", NC_COMPOUND_OFFSET(Sample, subRegion), NC_UINT);
        check_err(stat, __LINE__, __FILE__);
      }

      stat = nc_def_var(ncid, "offset", NC_INT64, 1, &rankdim, &offsetid);
      check_err(stat, __LINE__, __FILE__);
      stat = nc_def_var(ncid, "sample", sampletyp, 1, &sampledim, &sampleid);
      check_err(stat, __LINE__, __FILE__);

      stat = nc_enddef(ncid);
      check_err(stat, __LINE__, __FILE__);

      stat = nc_var_par_access(ncid, offsetid, NC_COLLECTIVE);
      check_err(stat, __LINE__, __FILE__);
      stat = nc_var_par_access(ncid, sampleid, NC_COLLECTIVE);
      check_err(stat, __LINE__, __FILE__);

      std::size_t start = 0;
      std::size_t count = 0;
      long offsetData[2];
      if (rank == 0) {
        start = 0;
        count = 2;
        offsetData[0] = 0;
      } else {
        start = 1 + rank;
        count = 1;
      }
      offsetData[count - 1] = sampleOffset;

      stat = nc_put_vara_long(ncid, offsetid, &start, &count, offsetData);
      check_err(stat, __LINE__, __FILE__);

      start = sampleOffset - nSamples;
      count = nSamples;
      stat = nc_put_vara(ncid, sampleid, &start, &count, regions[region].times.data());
      check_err(stat, __LINE__, __FILE__);

      stat = nc_close(ncid);
      check_err(stat, __LINE__, __FILE__);
    }
    logInfo() << "Finished writing loop statistics samples.";
#else
    logWarning() << "Writing loop statistics requires NetCDF and MPI.";
#endif
  }
}

} // namespace seissol
