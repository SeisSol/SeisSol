// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_RUNSUMMARY_H_
#define SEISSOL_SRC_MONITORING_RUNSUMMARY_H_

#include "Common/Json.h"
#include "IO/Instance/Point/Csv.h"
#include "Numerical/Statistics.h"

#include <chrono>
#include <optional>
#include <string>
#include <string_view>

namespace seissol::monitoring {

//! @brief What one simulation epoch, from one synchronization point to the next, took.
struct EpochRecord {
  double timeFrom{0};
  double timeTo{0};
  //! the wall time of the epoch, its compute part (advancing the clusters, waiting included) and
  //! its blocking IO, over the ranks (complete on rank 0 only)
  statistics::Summary wall;
  statistics::Summary compute;
  statistics::Summary io;
  //! the rates of the epoch, summed over the ranks: hardware and non-zero FLOP/s, kernel bytes/s
  double hardwareFlopPerSecond{0};
  double nonzeroFlopPerSecond{0};
  double kernelBytesPerSecond{0};
};

/**
 * @brief What a run is, how it went and how long it took, for programs to read.
 *
 * A JSON file next to the output, <prefix>-summary.json, which rank 0 writes. It is rewritten as
 * the run moves on: at every stage of the initialization, during the simulation and at the end,
 * each time replacing the file in one step, so that a reader always finds a complete file. A run
 * that ends early, by an error or by the time limit of its job, leaves the summary of the last
 * stage it reached: its status stays "running", and its stage says how far it got.
 *
 * Next to it, <prefix>-epochs.csv holds a row per simulation epoch, which the summary only counts.
 *
 * Every rank may fill in content; only rank 0's is written. What the statistics over the ranks
 * give is complete on rank 0 only, which is all that is needed.
 */
class RunSummary {
  public:
  //! @brief The name and version of the format, which a reader checks first.
  static constexpr std::string_view Schema = "seissol-run-summary/1";

  /**
   * @brief Starts the summary, in the files next to the output with @p prefix , and writes it.
   *
   * What was filled in before is kept, behind the header. Collective: gathers the layout of the
   * run over all ranks.
   */
  void start(const std::string& prefix, const std::string& runId);

  //! @brief Whether start() was called; until then, nothing is written.
  [[nodiscard]] bool started() const { return started_; }

  //! @brief The member @p key of the top level of the summary.
  json::Value& operator[](std::string_view key) { return root_[key]; }

  //! @brief The summary as it stands.
  [[nodiscard]] const json::Value& content() const { return root_; }

  //! @brief Records that the run has reached @p stage , and writes the summary.
  void setStage(std::string_view stage);

  /**
   * @brief Records a simulation epoch: a row of the epoch table, and the progress in the summary.
   *
   * Both are written out every so often rather than after every epoch, which can be frequent; a
   * stage, finish() and flush() write everything out.
   */
  void addEpoch(const EpochRecord& epoch, double endTime);

  //! @brief Records that the run has ended as planned, and writes everything out.
  void finish();

  //! @brief Writes out the summary and the epoch rows not written yet, from rank 0.
  void flush();

  //! @brief Where the summary is written.
  [[nodiscard]] const std::string& path() const { return path_; }

  private:
  void writeSummary();

  std::string path_;
  std::string epochPath_;
  std::string runId_;
  json::Value root_;
  bool started_{false};
  std::size_t epochCount_{0};
  std::optional<io::instance::point::Csv> epochs_;
  std::chrono::steady_clock::time_point lastFlush_;
};

/**
 * @brief The minimum, maximum, mean, median, standard deviation and sum of @p summary.
 *
 * As statistics::parallelSummary produces it: complete on rank 0 only.
 */
json::Value toJson(const statistics::Summary& summary);

/**
 * @brief The statistics over all ranks of @p localValue , as toJson() writes them.
 *
 * Collective; complete on rank 0 only.
 */
json::Value statisticsOverRanks(double localValue);

/**
 * @brief The current time, in UTC, as ISO 8601 (2026-10-16T12:34:56Z).
 */
std::string utcTimestamp();

/**
 * @brief The environment variables that steer a run of SeisSol, by name.
 *
 * Those of SeisSol itself and of the libraries and runtimes it runs on (OpenMP, MPI, ASYNC, the
 * device runtimes, the job scheduler); not the whole environment, which can hold anything.
 */
json::Value relevantEnvironment();

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_RUNSUMMARY_H_
