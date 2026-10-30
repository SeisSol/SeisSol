// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "BuildInfo.h"
#include "Common/ConfigRegistry.h"
#include "Common/Filesystem.h"
#include "Common/Json.h"
#include "Common/RunId.h"
#include "Initializer/InitProcedure/Init.h"
#include "Initializer/Parameters/ParameterOverrides.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Preflight.h"
#include "Monitoring/Instrumentation.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/Runtime/Stream.h"
#include "SeisSol.h"

#include <cstdlib>
#include <ctime>
#include <exception>
#include <fty/fty.hpp>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utils/args.h>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/timeutils.h>
#include <vector>
#include <yaml-cpp/yaml.h>

#ifdef USE_ASAGI
#include "Reader/AsagiModule.h"
#endif

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

#if defined(__GNUC__) || defined(__linux__)
#include <cfenv>
#include <fenv.h>
#endif

namespace {
using namespace seissol;

/**
 * The ID of this run: the one SEISSOL_RUN_ID gives, so that whoever started the run can find its
 * files again, or a new one. Collective; rank 0 decides it for everyone.
 */
std::string determineRunId(utils::Env& env) {
  std::string runId;
  if (seissol::Mpi::mpi.rank() == 0) {
    const auto given = env.getOptional<std::string>("RUN_ID");
    if (given.has_value()) {
      const auto parsed = seissol::parseRunId(given.value());
      if (!parsed.has_value()) {
        logError() << "SEISSOL_RUN_ID is set to" << given.value()
                   << ", which is not a UUID (e.g. 123e4567-e89b-42d3-a456-426614174000).";
      }
      runId = parsed.value_or(std::string{});
    } else {
      runId = seissol::newRunId();
    }
  }
  seissol::Mpi::mpi.broadcastContainer(runId, 0);
  return runId;
}

/**
 * Takes the --set options out of the command line, in the order they were given.
 *
 * utils::Args keeps only the last value of an option given several times, and an override is
 * meant to be repeated; what remains is left for utils::Args.
 */
std::vector<std::string> extractOverrides(int& argc, char** argv) {
  constexpr std::string_view Option = "--set";
  std::vector<std::string> overrides;
  int kept = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string_view argument(argv[i]);
    if (argument == Option) {
      if (i + 1 == argc) {
        logError() << "--set needs a value: --set section.field=value";
      }
      overrides.emplace_back(argv[i + 1]);
      ++i;
    } else if (argument.substr(0, Option.size() + 1) == std::string(Option) + "=") {
      overrides.emplace_back(argument.substr(Option.size() + 1));
    } else {
      argv[kept] = argv[i];
      ++kept;
    }
  }
  argc = kept;
  argv[argc] = nullptr;
  return overrides;
}

/**
 * How the run was started, for its summary.
 */
json::Value invocation(const char* executable,
                       const std::string& parameterFile,
                       const std::vector<initializer::parameters::ParameterOverride>& overrides,
                       const std::optional<std::string>& checkpoint) {
  json::Value value;
  std::error_code error;
  // the executable as it is on disk, which argv[0] need not name
  const auto self = filesystem::read_symlink("/proc/self/exe", error);
  value["executable"] = error ? std::string(executable) : self.string();
  const auto absolute = filesystem::absolute(parameterFile, error);
  value["parameterFile"] = error ? parameterFile : absolute.string();
  const auto workingDirectory = filesystem::current_path(error);
  value["workingDirectory"] = error ? std::string() : workingDirectory.string();
  auto& given = value["overrides"] = json::Value::object();
  for (const auto& entry : overrides) {
    given[entry.name()] = entry.value;
  }
  value["checkpoint"] = checkpoint.has_value() ? json::Value(checkpoint.value()) : json::Value();
  return value;
}

//! @brief The command line of SeisSol, for parsing it and for its help.
utils::Args commandLine() {
  utils::Args args("SeisSol is a scientific software for the numerical simulation of seismic wave "
                   "phenomena and earthquake dynamics. This version of SeisSol was built for the "
                   "following configurations:\n" +
                   describeBuiltConfigs());
  args.addAdditionalOption("parameterfile", "The parameter file", false);
  args.addOption(
      "checkpoint", 'c', "The checkpoint file to restart from", utils::Args::Required, false);
  // taken out of the command line before parsing it (see extractOverrides); here for the help
  args.addOption("set",
                 0,
                 "Overrides a parameter of the parameter file, as section.field=value; may be "
                 "given several times",
                 utils::Args::Required,
                 false);
  // likewise taken out before (see extractCheckOption)
  args.addOption("check",
                 0,
                 "Checks the run instead of running it: the parameters, the files they name, the "
                 "mesh and the output directory, on one process; --check=<file> writes a report",
                 utils::Args::Optional,
                 false);
  return args;
}

/**
 * Checks the run that the command line describes, without starting it: on one process, without
 * MPI and devices.
 */
int checkMain(int argc, char* argv[], const std::string& reportPath) {
  logInfo() << "Checking a run of SeisSol" << BuildInfo::VersionString.c_str()
            << "without running it.";
  const auto overrides = extractOverrides(argc, argv);
  auto args = commandLine();
  if (args.parse(argc, argv) != utils::Args::Success) {
    return 1;
  }
  const auto* parameterFile = args.getAdditionalArgument("parameterfile", "parameters.par");
  const auto checkpoint =
      args.isSet("checkpoint")
          ? std::make_optional<std::string>(args.getArgument<const char*>("checkpoint"))
          : std::nullopt;
  return seissol::initializer::preflight(parameterFile, overrides, checkpoint, reportPath);
}

std::shared_ptr<YAML::Node> readYamlParams(const std::string& parameterFile) {
  // Read parameter file input from file
  fty::Loader<fty::AsLowercase> loader{};
  std::shared_ptr<YAML::Node> inputParams = nullptr;
  try {
    inputParams = std::make_shared<YAML::Node>(loader.load(parameterFile));
  } catch (const std::exception& error) {
    logError() << "Error while reading the parameter file:" << std::string(error.what());
  }
  return inputParams;
}
} // namespace

int main(int argc, char* argv[]) {
  try {
    // a check runs on its own, before anything of a run starts
    const auto checkReport = seissol::initializer::extractCheckOption(argc, argv);
    if (checkReport.has_value()) {
      return checkMain(argc, argv, checkReport.value());
    }

#ifdef ACL_DEVICE
    if (seissol::Mpi::requestedDataTransferMode() == seissol::Mpi::DataTransferMode::ActorMpi) {
      // the device then waits for the host, and the host for other processes: CUDA must not load
      // a kernel on its first launch, as it may wait for the device meanwhile (unless set already)
      setenv("CUDA_MODULE_LOADING", "EAGER", 0);
    }
    // queues of their own for the stream runtimes; not with stream-mpi, which is untried with them
    // (CCL and SHMEM get streams of HIP, which wait for the own streams on the device, see
    // StreamExchangeScheduler)
    utils::Env deviceEnv("SEISSOL_");
    const auto requestedMode = seissol::Mpi::requestedDataTransferMode();
    const bool streamLibraries = requestedMode.has_value() &&
                                 requestedMode == seissol::Mpi::DataTransferMode::DirectStreamMpi;
    auto ownQueues = streamLibraries ? 0 : seissol::ownQueueCount(deviceEnv);
    if (ownQueues > 0 && requestedMode == seissol::Mpi::DataTransferMode::ActorMpi &&
        ownQueues % 2 != 0) {
      // the copy layers wait for the halo exchange there; with the interior layers on the same
      // queues, runs hung (TPV13 on LUMI, one queue for all clusters): an even number keeps them
      // apart (see ownQueueOf())
      ++ownQueues;
    }
    if (ownQueues > 0) {
      // the streams of HIP then only record graphs and hold allocations: one hardware queue each
      // leaves the others of the device to the own queues
      setenv("GPU_MAX_HW_QUEUES", "1", 0);
      // the queue of the other stream runtimes and that of the device timeline are idle most of
      // the time: they may share a pipe of the command processor with a queue of the clusters
      // (the others get pipes of their own as long as there are, see the Device library)
      const auto light = "0," + std::to_string(seissol::parallel::runtime::TimelineOwnQueue);
      setenv("DEVICE_OWN_QUEUE_LIGHT", light.c_str(), 0);
    }
    seissol::Mpi::mpi.bindAcceleratorDevice();
    const device::DeviceInstance& device = device::DeviceInstance::instance();
    device.api().initialize();
    if (ownQueues > 0 && device.api().enableOwnQueues()) {
      seissol::parallel::runtime::ownQueueCount() = ownQueues;
      seissol::parallel::runtime::ownFaceQueueCount() = seissol::ownFaceQueueCount(deviceEnv);
      seissol::parallel::runtime::ownBranchQueueCount() = seissol::ownBranchQueueCount(deviceEnv);
      const auto priorities = deviceEnv.get<std::string>("OWN_QUEUE_PRIORITIES", "");
      seissol::parallel::runtime::ownQueuePriorities() = {
          priorities.find("copy") != std::string::npos,
          priorities.find("faces") != std::string::npos};
      if (deviceEnv.get<std::string>("OWN_QUEUE_LAYOUT", "layers") == "halo") {
        seissol::parallel::runtime::ownQueueLayout() =
            seissol::parallel::runtime::OwnQueueLayout::Halo;
      }
    }
#endif // ACL_DEVICE

    utils::Env env("SEISSOL_");

#ifdef USE_ASAGI
    // Construct an instance of AsagiModule, to initialize it.
    // Its MPI mode is derived from the environment at construction time, i.e. before MPI
    // is set up.
    seissol::asagi::AsagiModule::initInstance(env);
#endif

    seissol::Mpi::mpi.init(argc, argv);
    const int rank = seissol::Mpi::mpi.rank();

    utils::Logger::setRank(rank);
#ifdef ACL_DEVICE
    if (seissol::parallel::runtime::ownQueueCount() > 0) {
      logInfo() << "The clusters run their work on" << seissol::parallel::runtime::ownQueueCount()
                << "queues of their own, the other stream runtimes on one more"
                << "(SEISSOL_OWN_QUEUES).";
      if (seissol::parallel::runtime::ownBranchQueueCount() > 0) {
        logInfo() << "The branches of the work run on"
                  << seissol::parallel::runtime::ownBranchQueueCount()
                  << "more (SEISSOL_OWN_QUEUE_BRANCHES).";
      }
      if (seissol::parallel::runtime::ownFaceQueueCount() > 0) {
        logInfo() << "The face work runs on" << seissol::parallel::runtime::ownFaceQueueCount()
                  << "more (SEISSOL_OWN_QUEUES_FACES).";
      }
    } else if (seissol::ownQueueCount(env) > 0) {
      logWarning() << "Queues of their own (SEISSOL_OWN_QUEUES) are not available: they need HIP,"
                   << "and a transfer mode other than stream-mpi.";
    }
#endif

    const auto runId = determineRunId(env);

    LIKWID_MARKER_INIT;
#pragma omp parallel
    {
      LIKWID_MARKER_THREADINIT;
      LIKWID_MARKER_REGISTER("SeisSol");
      LIKWID_MARKER_REGISTER("computeDynamicRuptureFrictionLaw");
      LIKWID_MARKER_REGISTER("computeDynamicRupturePostHook");
      LIKWID_MARKER_REGISTER("computeDynamicRupturePreHook");
      LIKWID_MARKER_REGISTER("computeDynamicRuptureFinalizeImposedState");
      LIKWID_MARKER_REGISTER("computeDynamicRuptureSpaceTimeInterpolation");
      LIKWID_MARKER_REGISTER("computeDynamicRuptureTimeStepLoop");
    }

#pragma omp parallel
    {
      LIKWID_MARKER_START("SeisSol");
    }

    EPIK_TRACER("SeisSol");
    SCOREP_USER_REGION("SeisSol", SCOREP_USER_REGION_TYPE_FUNCTION);

    // Print welcome message
    logInfo() << "Welcome to SeisSol";
    logInfo() << "Copyright (c) 2012 -" << BuildInfo::CommitYear.c_str() << " SeisSol Group";
    logInfo() << "Version:" << BuildInfo::VersionString.c_str();
    logInfo() << "Built on:" << __DATE__ << __TIME__;
    logInfo() << "Last commit:" << BuildInfo::CommitHash.c_str() << "at"
              << BuildInfo::CommitTimestamp.c_str();
    logInfo() << "Compiled with HOST_ARCH =" << BuildInfo::SeisSolHostArch.c_str();
#ifdef ACL_DEVICE
    logInfo() << "Compiled with DEVICE_BACKEND =" << BuildInfo::SeisSolDeviceBackend.c_str();
    logInfo() << "Compiled with DEVICE_ARCH =" << BuildInfo::SeisSolDeviceArch.c_str();
#endif
    logInfo() << "Run ID:" << runId.c_str();

    if (env.get<bool>("FLOATING_POINT_EXCEPTION", false)) {
      // Check if on a GNU system (Linux) or other platform
#if defined(__GNUC__) && defined(__linux__)
      // a CPU that cannot trap on floating-point exceptions refuses to (RISC-V, many ARM CPUs)
      if (feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT) == -1) {
        logWarning() << "Floating-point exception handling was requested, but the CPU does not "
                        "trap on floating-point exceptions.";
      } else {
        logInfo() << "Enabling floating point exception handlers.";
      }
#else
      logInfo() << "Floating-point exception handling was requested, but is not supported on this "
                   "platform.";
#endif
    }

    const auto overrideTexts = extractOverrides(argc, argv);
    std::vector<seissol::initializer::parameters::ParameterOverride> overrides;
    for (const auto& text : overrideTexts) {
      const auto parsed = seissol::initializer::parameters::parseOverride(text);
      if (!parsed.has_value()) {
        logError() << "The override" << text << "is not of the form section.field=value.";
      }
      overrides.push_back(parsed.value_or(seissol::initializer::parameters::ParameterOverride{}));
    }

    auto args = commandLine();
    switch (args.parse(argc, argv)) {
    case utils::Args::Help: {
      [[fallthrough]];
    }
    case utils::Args::Error: {
      seissol::Mpi::finalize();
      return 1;
    }
    case utils::Args::Success: {
      break;
    }
    }
    const auto* parameterFile = args.getAdditionalArgument("parameterfile", "parameters.par");
    logInfo() << "Using the parameter file" << parameterFile;
    // read parameter file input
    const auto yamlParams = readYamlParams(parameterFile);
    for (const auto& entry : overrides) {
      logInfo() << "Overriding" << entry.name().c_str() << "with" << entry.value;
    }
    seissol::initializer::parameters::applyOverrides(*yamlParams, overrides);
    seissol::initializer::parameters::ParameterReader parameterReader(
        *yamlParams, parameterFile, false);
    auto parameters = [&]() {
      try {
        return seissol::initializer::parameters::readSeisSolParameters(&parameterReader);
      } catch (const std::exception& error) {
        logError() << "Uncaught error while parsing the parameter file:"
                   << std::string(error.what());
      }
      return initializer::parameters::SeisSolParameters();
    }();
    parameterReader.warnUnknown();
    for (const auto& entry :
         seissol::initializer::parameters::unknownOverrides(parameterReader, overrides)) {
      logError() << "The override" << entry.name()
                 << "names a parameter that SeisSol does not read; check its spelling.";
    }

    {
      // Initialize SeisSol
      seissol::SeisSol seissolInstance(parameters, env);
#ifdef USE_ASAGI
      auto& asagiModule = seissol::asagi::AsagiModule::getInstance();
      asagiModule.setPinning(seissolInstance.pinning());
      asagiModule.registerHooks(seissolInstance.modules());
#endif
      seissolInstance.setRunId(runId);
      seissolInstance.runSummary()["invocation"] = invocation(
          argv[0],
          parameterFile,
          overrides,
          args.isSet("checkpoint")
              ? std::make_optional<std::string>(args.getArgument<const char*>("checkpoint"))
              : std::nullopt);

      if (args.isSet("checkpoint")) {
        const auto* checkpointFile = args.getArgument<const char*>("checkpoint");
        seissolInstance.loadCheckpoint(checkpointFile);
      }

      // run SeisSol
      const bool runSeisSol = seissolInstance.init();

      const auto stamp = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr));
      seissolInstance.setBackupTimeStamp(stamp);

      // Run SeisSol
      if (runSeisSol) {
        seissol::initializer::initprocedure::seissolMain(seissolInstance);
      }

#pragma omp parallel
      {
        LIKWID_MARKER_STOP("SeisSol");
      }

      LIKWID_MARKER_CLOSE;
      // Finalize SeisSol
      seissolInstance.finalize();
    }

#ifdef ACL_DEVICE
    device.api().finalize();
#endif
    return 0;
  } catch (const std::exception& error) {
    logError() << "An uncaught error occurred while running SeisSol:" << std::string(error.what());
  }
  return 1;
}
