// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MONITORING_LOOPSTATISTICS_H_
#define SEISSOL_SRC_MONITORING_LOOPSTATISTICS_H_

#include "Parallel/MPI.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <unordered_map>
#include <vector>

namespace seissol {
class LoopStatistics {
  public:
  void enableSampleOutput(bool enabled);

  void addRegion(const std::string& name, bool includeInSummary = true);

  [[nodiscard]] unsigned getRegion(const std::string& name) const;

  void begin(unsigned region);

  void end(unsigned region, unsigned numIterations, unsigned subRegion);

  void addSample(
      unsigned region, unsigned numIterations, unsigned subRegion, timespec begin, timespec end);

  void reset();

  void printSummary(MPI_Comm comm);

  void writeSamples(const std::string& outputPrefix, bool isLoopStatisticsNetcdfOutputOn);

  private:
  struct Sample {
    timespec begin;
    timespec end;
    unsigned numIters;
    unsigned subRegion;
  };

  struct StatisticVariables {
    double x = 0;
    double x2 = 0;
    double xy = 0;
    double y = 0;
    double y2 = 0;
    unsigned long long n = 0;
  };

  struct Region {
    std::string name;
    std::vector<Sample> times;
    bool includeInSummary;
    timespec begin{};
    StatisticVariables variables;

    Region(const std::string& name, bool includeInSummary);
  };

  std::vector<Region> regions;
  bool outputSamples = false;
};
} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_LOOPSTATISTICS_H_
