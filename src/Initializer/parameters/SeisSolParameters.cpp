#include "SeisSolParameters.h"

namespace seissol::initializers::parameters {

SeisSolParameters readSeisSolParameters(ParameterReader* parameterReader) {
  // ParameterReader baseReader(baseNode, false);
  const CubeGeneratorParameters cubeGeneratorParameters =
      readCubeGeneratorParameters(parameterReader);
  const DRParameters drParameters = readDRParameters(parameterReader);
  const InitializationParameters initializationParameters =
      readInitializationParameters(parameterReader);
  const MeshParameters meshParameters = readMeshParameters(parameterReader);
  const ModelParameters modelParameters = readModelParameters(parameterReader);
  const OutputParameters outputParameters = readOutputParameters(parameterReader);
  const SourceParameters sourceParameters = readSourceParameters(parameterReader);
  const TimeSteppingParameters timeSteppingParameters = readTimeSteppingParameters(parameterReader);

  // TODO: Add sanity check again

  // if (isModelViscoelastic()) {
  //   double maxTimestepWidthDefault =
  //       0.25 / (seissolParams.model.freqCentral * std::sqrt(seissolParams.model.freqRatio));
  //   if (reader->hasField("fixtimestep")) {
  //     seissolParams.timeStepping.maxTimestepWidth =
  //         reader->readWithDefault("fixtimestep", maxTimestepWidthDefault);
  //     if (seissolParams.timeStepping.maxTimestepWidth > maxTimestepWidthDefault) {
  //       logWarning(seissol::MPI::mpi.rank())
  //           << "The given maximum timestep width (fixtimestep) is set to"
  //           << seissolParams.timeStepping.maxTimestepWidth
  //           << "which is larger than the recommended value of" << maxTimestepWidthDefault
  //           << " for visco-elastic material (as specified in the documentation). Please be aware
  //           "
  //              "that a too large maximum timestep width may cause the solution to become
  //              unstable.";
  //     } else {
  //       logInfo(seissol::MPI::mpi.rank())
  //           << "Maximum timestep width (fixtimestep) given as"
  //           << seissolParams.timeStepping.maxTimestepWidth << "(less or equal to reference
  //           timestep"
  //           << maxTimestepWidthDefault << ")";
  //     }
  //   } else {
  //     seissolParams.timeStepping.maxTimestepWidth = maxTimestepWidthDefault;
  //     logInfo(seissol::MPI::mpi.rank())
  //         << "Setting maximum timestep width to" << maxTimestepWidthDefault
  //         << " for visco-elastic material (as specified in the documentation).";
  //   }
  // } else {
  //   seissolParams.timeStepping.maxTimestepWidth = reader->readWithDefault("fixtimestep", 5000.0);
  // }
  logInfo(seissol::MPI::mpi.rank()) << "Reading SeisSol parameter file...";

  parameterReader->warnDeprecated({"boundaries",
                                   "rffile",
                                   "inflowbound",
                                   "inflowboundpwfile",
                                   "inflowbounduin",
                                   "source110",
                                   "source15",
                                   "source1618",
                                   "source17",
                                   "source19",
                                   "spongelayer",
                                   "sponges",
                                   "analysis",
                                   "analysisfields",
                                   "debugging"});

  logInfo(seissol::MPI::mpi.rank()) << "SeisSol parameter file read successfully.";

  auto printYesNo = [](bool yesno) { return yesno ? "yes" : "no"; };

  logInfo(seissol::MPI::mpi.rank()) << "Model information:";
  logInfo(seissol::MPI::mpi.rank()) << "Elastic model:" << printYesNo(isModelElastic());
  logInfo(seissol::MPI::mpi.rank()) << "Viscoelastic model:" << printYesNo(isModelViscoelastic());
  logInfo(seissol::MPI::mpi.rank()) << "Anelastic model:" << printYesNo(isModelAnelastic());
  logInfo(seissol::MPI::mpi.rank()) << "Poroelastic model:" << printYesNo(isModelPoroelastic());
  logInfo(seissol::MPI::mpi.rank()) << "Anisotropic model:" << printYesNo(isModelAnisotropic());
  logInfo(seissol::MPI::mpi.rank()) << "Plasticity:" << printYesNo(modelParameters.plasticity);

  return SeisSolParameters{cubeGeneratorParameters,
                           drParameters,
                           initializationParameters,
                           meshParameters,
                           modelParameters,
                           outputParameters,
                           sourceParameters,
                           timeSteppingParameters};
}
} // namespace seissol::initializers::parameters
