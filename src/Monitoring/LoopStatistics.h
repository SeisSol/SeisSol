// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MONITORING_LOOPSTATISTICS_H_
#define SEISSOL_SRC_MONITORING_LOOPSTATISTICS_H_

#include "Common/Json.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <time.h>
#include <unordered_map>
#include <vector>

namespace seissol {
class LoopStatistics {
  public:
  void addRegion(const std::string& name, bool includeInSummary = true);

  [[nodiscard]] unsigned getRegion(const std::string& name) const;

  /**
   * The time a loop of a region starts at, for end(). Loops of several clusters can run at the same
   * time (e.g. in a task graph), so each one keeps its own.
   */
  [[nodiscard]] static timespec begin();

  void end(unsigned region, unsigned numIterations, unsigned subRegion, const timespec& begin);

  void addSample(
      unsigned region, unsigned numIterations, unsigned subRegion, timespec begin, timespec end);

  void reset();

  /**
   * @brief Prints the time spent in the regions, and a regression of it over the iterations.
   *
   * Collective. Returns what it printed, for the run summary (complete on rank 0 only).
   */
  json::Value printSummary(MPI_Comm comm);

  //! @brief The time this rank spent in the regions that count towards the total, the kernels.
  [[nodiscard]] double kernelTime() const;

  private:
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
    bool includeInSummary;
    StatisticVariables variables;

    Region(const std::string& name, bool includeInSummary);
  };

  std::vector<Region> regions_;
  // the samples of loops that run at the same time
  std::mutex samplesMutex_;
};
} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_LOOPSTATISTICS_H_
