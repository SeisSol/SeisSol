/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/PreProcessorMacros.h"
#include "Modules/Modules.h"
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fty/fty.hpp>
#include <memory>
#include <ostream>
#include <string>
#include <utils/logger.h>
#include <utils/timeutils.h>
#include <xdmfwriter/scorep_wrapper.h>
#include <yaml-cpp/yaml.h>

#include "Initializer/InitProcedure/Init.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "SeisSol.h"
#include "utils/args.h"

#ifdef USE_ASAGI
#include "Reader/AsagiModule.h"
#endif

#ifdef ACL_DEVICE
#include "device.h"
#endif

std::shared_ptr<YAML::Node> readYamlParams(const std::string& parameterFile) {
  // Read parameter file input from file
  fty::Loader<fty::AsLowercase> loader{};
  std::shared_ptr<YAML::Node> inputParams = nullptr;
  try {
    inputParams = std::make_shared<YAML::Node>(loader.load(parameterFile));
  } catch (const std::exception& error) {
    logError() << "Error while reading the parameter file:" << std::string(error.what())
               << std::endl;
  }
  return inputParams;
}

int main(int argc, char* argv[]) {
#ifdef ACL_DEVICE
#ifdef USE_MPI
  seissol::MPI::mpi.bindAcceleratorDevice();
#endif // USE_MPI
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->initialize();
  device.api->allocateStackMem();
#endif // ACL_DEVICE

#ifdef USE_ASAGI
  // Construct an instance of AsagiModule, to initialize it.
  // It needs to be done here, as it registers PRE_MPI hooks
  seissol::asagi::AsagiModule::getInstance();
#endif
  // Call pre MPI hooks
  seissol::Modules::callHook<ModuleHook::PreMPI>();

  seissol::MPI::mpi.init(argc, argv);
  const int rank = seissol::MPI::mpi.rank();

  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("SeisSol");
    LIKWID_MARKER_REGISTER("computeDynamicRuptureFrictionLaw");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePostHook");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePostcomputeImposedState");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePreHook");
    LIKWID_MARKER_REGISTER("computeDynamicRupturePrecomputeStress");
    LIKWID_MARKER_REGISTER("computeDynamicRuptureSpaceTimeInterpolation");
    LIKWID_MARKER_REGISTER("computeDynamicRuptureUpdateFrictionAndSlip");
  }

#pragma omp parallel
  { LIKWID_MARKER_START("SeisSol"); }

  EPIK_TRACER("SeisSol");
  SCOREP_USER_REGION("SeisSol", SCOREP_USER_REGION_TYPE_FUNCTION);

  // TODO Read parameters here
  // Parse command line arguments
  utils::Args args("SeisSol is a scientific software for the numerical simulation of seismic wave "
                   "phenomena and earthquake dynamics.");
  args.addAdditionalOption("file", "The parameter file", false);
  args.addOption(
      "checkpoint", 'c', "The checkpoint file to restart from", utils::Args::Optional, false);
  switch (args.parse(argc, argv)) {
  case utils::Args::Help: {
    [[fallthrough]];
  }
  case utils::Args::Error: {
    seissol::MPI::mpi.finalize();
    exit(1);
    break;
  }
  case utils::Args::Success: {
    break;
  }
  }
  const auto parameterFile = args.getAdditionalArgument("file", "parameters.par");
  logInfo(rank) << "Using the parameter file" << parameterFile;
  // read parameter file input
  const auto yamlParams = readYamlParams(parameterFile);
  seissol::initializer::parameters::ParameterReader parameterReader(
      *yamlParams.get(), parameterFile, false);
  auto parameters = seissol::initializer::parameters::readSeisSolParameters(&parameterReader);
  parameterReader.warnUnknown();

  // Initialize SeisSol
  seissol::SeisSol seissolInstance(parameters);

  if (args.isSet("checkpoint")) {
    const auto checkpointFile = args.getArgument<const char*>("checkpoint");
    seissolInstance.loadCheckpoint(checkpointFile);
    logInfo(rank) << "Using the checkpoint file" << checkpointFile;
  }

  // run SeisSol
  const bool runSeisSol = seissolInstance.init(argc, argv);

  const auto stamp = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(0L));
  seissolInstance.setBackupTimeStamp(stamp);

  // Run SeisSol
  if (runSeisSol) {
    seissol::initializer::initprocedure::seissolMain(seissolInstance);
  }

#pragma omp parallel
  { LIKWID_MARKER_STOP("SeisSol"); }

  LIKWID_MARKER_CLOSE;
  // Finalize SeisSol
  seissolInstance.finalize();

#ifdef ACL_DEVICE
  device.api->finalize();
#endif
  return 0;
}
