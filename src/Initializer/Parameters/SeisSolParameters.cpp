#include "SeisSolParameters.h"

namespace seissol::initializer::parameters {

SeisSolParameters readSeisSolParameters(ParameterReader* parameterReader) {
  logInfo(seissol::MPI::mpi.rank()) << "Reading SeisSol parameter file...";

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
  logInfo(seissol::MPI::mpi.rank()) << "Non Linear Damage Model:" << printYesNo(isModelDamagedElastic());

  return SeisSolParameters{cubeGeneratorParameters,
                           drParameters,
                           initializationParameters,
                           meshParameters,
                           modelParameters,
                           outputParameters,
                           sourceParameters,
                           timeSteppingParameters};
}
} // namespace seissol::initializer::parameters
