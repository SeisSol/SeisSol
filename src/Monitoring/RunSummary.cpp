// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RunSummary.h"

#include "BuildInfo.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Filesystem.h"
#include "Common/Json.h"
#include "Common/RunId.h"
#include "IO/Instance/Point/Csv.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"
#include "Parallel/OpenMP.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <limits>
#include <mpi.h>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef __APPLE__
#include <crt_externs.h>
#else
// as POSIX declares it
// NOLINTNEXTLINE(readability-redundant-declaration,cppcoreguidelines-avoid-non-const-global-variables)
extern char** environ;
#endif

namespace seissol::monitoring {

namespace {
//! @brief The environment of the process, as NAME=value entries ending in a null pointer.
const char* const* environmentVariables() {
#ifdef __APPLE__
  // a shared library cannot refer to environ on macOS
  return *_NSGetEnviron();
#else
  return environ;
#endif
}

//! @brief The first line of what MPI_Get_library_version says, which names the library.
std::string mpiLibrary() {
  std::array<char, MPI_MAX_LIBRARY_VERSION_STRING> buffer{};
  int length = 0;
  MPI_Get_library_version(buffer.data(), &length);
  std::string version(buffer.data(), static_cast<std::size_t>(length));
  const auto end = version.find_first_of(std::string("\n\0", 2));
  if (end != std::string::npos) {
    version.resize(end);
  }
  while (!version.empty() && (version.back() == ' ' || version.back() == ',')) {
    version.pop_back();
  }
  return version;
}

} // namespace

json::Value buildInformation() {
  json::Value build;
  build["version"] = BuildInfo::VersionString;
  build["commit"] = BuildInfo::CommitHash;
  build["commitTimestamp"] = BuildInfo::CommitTimestamp;
  build["hostArch"] = BuildInfo::SeisSolHostArch;
  build["deviceBackend"] = BuildInfo::SeisSolDeviceBackend;
  build["deviceArch"] = BuildInfo::SeisSolDeviceArch;
  build["compiler"]["id"] = BuildInfo::CompilerId;
  build["compiler"]["version"] = BuildInfo::CompilerVersion;
  build["buildType"] = BuildInfo::BuildType;
  auto& configurations = build["configurations"] = json::Value::array();
  for (std::size_t id = 0; id < builtConfigCount(); ++id) {
    configurations.push(configName(configValue(static_cast<ConfigId>(id))));
  }
  return build;
}

namespace {
//! @brief The layout of the run over nodes and threads. Collective.
json::Value parallelSection() {
  json::Value parallel;
  const auto& mpi = seissol::Mpi::mpi;
  parallel["ranks"] = mpi.size();
  if (mpi.rank() == 0) {
    const auto& hostNames = seissol::Mpi::mpi.getHostNames();
    parallel["nodes"] = std::set<std::string>(hostNames.begin(), hostNames.end()).size();
  }
  parallel["ranksPerNode"] = statisticsOverRanks(static_cast<double>(mpi.sharedMemMpiSize()));
  parallel["threadsPerRank"] = statisticsOverRanks(static_cast<double>(OpenMP::threadCount()));
  parallel["mpiLibrary"] = mpiLibrary();
  return parallel;
}
} // namespace

json::Value toJson(const statistics::Summary& summary) {
  json::Value value;
  value["min"] = summary.min;
  value["max"] = summary.max;
  value["mean"] = summary.mean;
  value["median"] = summary.median;
  value["std"] = summary.std;
  value["sum"] = summary.sum;
  return value;
}

json::Value statisticsOverRanks(double localValue) {
  return toJson(statistics::parallelSummary(localValue));
}

std::string utcTimestamp() {
  const auto now = std::time(nullptr);
  std::tm utc{};
#ifdef _WIN32
  gmtime_s(&utc, &now);
#else
  gmtime_r(&now, &utc);
#endif
  std::array<char, 32> buffer{};
  const auto length = std::strftime(buffer.data(), buffer.size(), "%Y-%m-%dT%H:%M:%SZ", &utc);
  return {buffer.data(), length};
}

json::Value relevantEnvironment() {
  // the knobs of SeisSol and of what it runs on; the job variables that say where it ran
  constexpr std::array<std::string_view, 20> Prefixes{
      "SEISSOL_", "OMP_",      "KMP_",  "GOMP_",   "ASYNC_", "XDMFWRITER_", "MPICH_",
      "I_MPI_",   "OMPI_MCA_", "UCX_",  "FI_",     "HSA_",   "HIP_",        "ROCR_",
      "CUDA_",    "ZE_",       "SYCL_", "ONEAPI_", "NCCL_",  "LIKWID_"};
  constexpr std::array<std::string_view, 12> Names{"SLURM_JOB_ID",
                                                   "SLURM_JOB_NAME",
                                                   "SLURM_JOB_NUM_NODES",
                                                   "SLURM_JOB_PARTITION",
                                                   "SLURM_CLUSTER_NAME",
                                                   "SLURM_NTASKS",
                                                   "SLURM_NTASKS_PER_NODE",
                                                   "SLURM_CPUS_PER_TASK",
                                                   "SLURM_CPU_BIND",
                                                   "PBS_JOBID",
                                                   "LSB_JOBID",
                                                   "FLUX_JOB_ID"};
  // what an MPI launcher sets for its own runtime rather than for the run
  constexpr std::array<std::string_view, 5> RuntimePrefixes{"OMPI_MCA_orte_",
                                                            "OMPI_MCA_ess",
                                                            "OMPI_MCA_pmix",
                                                            "OMPI_MCA_initial_wdir",
                                                            "OMPI_MCA_shmem_RUNTIME_QUERY_hint"};
  json::Value environment = json::Value::object();
  for (const char* const* entry = environmentVariables(); entry != nullptr && *entry != nullptr;
       ++entry) {
    const std::string_view variable(*entry);
    const auto equals = variable.find('=');
    if (equals == std::string_view::npos) {
      continue;
    }
    const auto name = variable.substr(0, equals);
    bool relevant = false;
    for (const auto prefix : Prefixes) {
      relevant = relevant || name.starts_with(prefix);
    }
    for (const auto known : Names) {
      relevant = relevant || name == known;
    }
    for (const auto prefix : RuntimePrefixes) {
      relevant = relevant && !name.starts_with(prefix);
    }
    if (relevant) {
      environment[name] = std::string(variable.substr(equals + 1));
    }
  }
  return environment;
}

void RunSummary::start(const std::string& prefix, const std::string& runId) {
  path_ = prefix + "-summary.json";
  epochPath_ = prefix + "-epochs.csv";
  rankPath_ = prefix + "-ranks.csv";
  runId_ = runId;

  // the header first, so that a reader sees what the file is before anything else
  json::Value root;
  root["schema"] = std::string(Schema);
  root["runId"] = runId;
  root["status"] = "running";
  root["stage"] = "init";
  root["started"] = utcTimestamp();
  root["updated"] = utcTimestamp();
  root["build"] = buildInformation();
  root["parallel"] = parallelSection();
  root["environment"] = relevantEnvironment();
  for (const auto& [key, value] : root_.isObject() ? root_.asObject() : json::Value::Object{}) {
    root[key] = value;
  }
  root_ = std::move(root);

  auto& epochs = epochs_.emplace("epochs");
  epochs.addColumn<std::uint64_t>("epoch");
  epochs.addColumn<double>("timeFrom");
  epochs.addColumn<double>("timeTo");
  epochs.addColumn<double>("wallMin");
  epochs.addColumn<double>("wallMean");
  epochs.addColumn<double>("wallMax");
  epochs.addColumn<double>("computeMean");
  epochs.addColumn<double>("computeMax");
  epochs.addColumn<double>("ioMean");
  epochs.addColumn<double>("ioMax");
  epochs.addColumn<double>("hardwareFlopPerSecond");
  epochs.addColumn<double>("nonzeroFlopPerSecond");
  epochs.addColumn<double>("kernelBytesPerSecond");
  epochs.addTextColumn("runId", RunIdLength);

  started_ = true;
  flush();
}

void RunSummary::setStage(std::string_view stage) {
  root_["stage"] = std::string(stage);
  flush();
}

void RunSummary::setRankValue(std::string_view column, double localValue) {
  for (auto& [name, value] : rankValues_) {
    if (name == column) {
      value = localValue;
      return;
    }
  }
  rankValues_.emplace_back(std::string(column), localValue);
}

json::Value RunSummary::overRanks(std::string_view column, double localValue) {
  setRankValue(column, localValue);
  return statisticsOverRanks(localValue);
}

void RunSummary::writeRanks() {
  auto& mpi = seissol::Mpi::mpi;

  // the columns of rank 0, which every rank fills in, in its order
  std::vector<char> names;
  for (const auto& [name, _] : rankValues_) {
    names.insert(names.end(), name.begin(), name.end());
    names.push_back('\n');
  }
  mpi.broadcastContainer(names, 0);
  std::vector<std::string> columns;
  std::string current;
  for (const auto character : names) {
    if (character == '\n') {
      columns.push_back(current);
      current.clear();
    } else {
      current.push_back(character);
    }
  }

  std::vector<double> values(columns.size(), std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < columns.size(); ++i) {
    for (const auto& [name, value] : rankValues_) {
      if (name == columns[i]) {
        values[i] = value;
      }
    }
  }
  const auto rankValues = mpi.collectContainer(values);
  const auto localRanks = mpi.collect(mpi.sharedMemMpiRank());

  if (!started_ || mpi.rank() != 0) {
    return;
  }
  const auto& hostNames = mpi.getHostNames();
  std::size_t hostWidth = 1;
  for (const auto& hostName : hostNames) {
    hostWidth = std::max(hostWidth, hostName.size());
  }

  io::instance::point::Csv table("ranks");
  table.addColumn<int>("rank");
  table.addTextColumn("host", hostWidth);
  table.addColumn<int>("localRank");
  for (const auto& column : columns) {
    table.addColumn<double>(column);
  }
  table.addTextColumn("runId", RunIdLength);
  for (int rank = 0; rank < mpi.size(); ++rank) {
    table.addCell<int>(rank);
    table.addText(static_cast<std::size_t>(rank) < hostNames.size() ? hostNames[rank] : "");
    table.addCell<int>(localRanks[rank]);
    for (const auto value : rankValues[rank]) {
      table.addCell<double>(value);
    }
    table.addText(runId_);
  }
  table.writeFile(rankPath_);
}

void RunSummary::addEpoch(const EpochRecord& epoch, double endTime) {
  ++epochCount_;
  auto& simulation = root_["simulation"];
  simulation["currentTime"] = epoch.timeTo;
  simulation["progress"] = endTime > 0 ? epoch.timeTo / endTime : 1.0;
  simulation["epochs"] = epochCount_;

  if (started_ && seissol::Mpi::mpi.rank() == 0) {
    auto& epochs = epochs_.value();
    epochs.addCell<std::uint64_t>(epochCount_);
    epochs.addCell<double>(epoch.timeFrom);
    epochs.addCell<double>(epoch.timeTo);
    epochs.addCell<double>(epoch.wall.min);
    epochs.addCell<double>(epoch.wall.mean);
    epochs.addCell<double>(epoch.wall.max);
    epochs.addCell<double>(epoch.compute.mean);
    epochs.addCell<double>(epoch.compute.max);
    epochs.addCell<double>(epoch.io.mean);
    epochs.addCell<double>(epoch.io.max);
    epochs.addCell<double>(epoch.hardwareFlopPerSecond);
    epochs.addCell<double>(epoch.nonzeroFlopPerSecond);
    epochs.addCell<double>(epoch.kernelBytesPerSecond);
    epochs.addText(runId_);
  }

  // an epoch can take a fraction of a second; a file system should not see a write for each
  constexpr auto FlushInterval = std::chrono::seconds(10);
  if (std::chrono::steady_clock::now() - lastFlush_ >= FlushInterval) {
    flush();
  }
}

void RunSummary::finish() {
  root_["status"] = "completed";
  root_["stage"] = "finished";
  flush();
}

void RunSummary::flush() {
  lastFlush_ = std::chrono::steady_clock::now();
  if (!started_ || seissol::Mpi::mpi.rank() != 0) {
    return;
  }
  if (epochCount_ > 0) {
    epochs_.value().appendFile(epochPath_);
  }
  writeSummary();
}

void RunSummary::writeSummary() {
  root_["updated"] = utcTimestamp();
  if (!writeFileAtomically(path_, root_.dump() + "\n")) {
    logWarning() << "Could not write the run summary" << path_;
  }
}

} // namespace seissol::monitoring
