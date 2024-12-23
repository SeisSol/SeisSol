#include "OutputParameters.h"
#include <Equations/Datastructures.h>
#include <Initializer/InputAux.h>
#include <Initializer/Parameters/ParameterReader.h>
#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::parameters {

void warnIntervalAndDisable(bool& enabled,
                            double interval,
                            const std::string& valName,
                            const std::string& intName) {
  if (enabled && interval <= 0) {
    auto intPhrase = valName + " = 0";
    logInfo(seissol::MPI::mpi.rank())
        << "In your parameter file, you have specified a non-positive interval for" << intName
        << ". This mechanism is deprecated and may be removed in a future version of "
           "SeisSol. Consider disabling the whole module by setting"
        << valName << "to 0 instead by adding" << intPhrase
        << "to the \"output\" section of your parameter file instead.";
    // still, replicate the old behavior.
    enabled = false;
  }
}

CheckpointParameters readCheckpointParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");

  auto enabled = reader->readWithDefault("checkpoint", true);
  double interval = 0.0;
  if (enabled) {
    interval = reader->readWithDefault("checkpointinterval", 0.0);
    warnIntervalAndDisable(enabled, interval, "checkpoint", "checkpointinterval");
  } else {
    reader->markUnused({"checkpointinterval"});
  }

  reader->warnDeprecated({"checkpointbackend", "checkpointfile"});

  return CheckpointParameters{enabled, interval};
}

ElementwiseFaultParameters readElementwiseParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("elementwise");

  const auto printTimeIntervalSec = reader->readWithDefault("printtimeinterval_sec", 1.0);
  const auto outputMaskString =
      reader->readWithDefault<std::string>("outputmask", "1 1 1 1 1 1 0 0 0 0 0 0");
  const std::array<bool, 12> outputMask = convertStringToArray<bool, 12>(outputMaskString, false);
  const auto refinementStrategy = reader->readWithDefaultEnum<FaultRefinement>(
      "refinement_strategy",
      FaultRefinement::None,
      {FaultRefinement::Triple, FaultRefinement::Quad, FaultRefinement::None});
  const int refinement = reader->readWithDefault("refinement", 2);
  reader->warnDeprecated({"printintervalcriterion"});

  const auto vtkorder = reader->readWithDefault("vtkorder", -1);

  return ElementwiseFaultParameters{
      printTimeIntervalSec, outputMask, refinementStrategy, refinement, vtkorder};
}

EnergyOutputParameters readEnergyParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");

  bool enabled = reader->readWithDefault("energyoutput", false);
  const auto interval = reader->readWithDefault("energyoutputinterval", VeryLongTime);
  warnIntervalAndDisable(enabled, interval, "energyoutput", "energyoutputinterval");

  const auto computeVolumeEnergiesEveryOutput =
      reader->readWithDefault("computevolumeenergieseveryoutput", 1);
  const auto terminalOutput = reader->readWithDefault("energyterminaloutput", false);
  const auto terminalPrecision = reader->readWithDefault("energyterminalprecision", 6);

  auto* abortCriteriaReader = baseReader->readSubNode("abortcriteria");
  const auto terminatorMaxTimePostRupture = abortCriteriaReader->readWithDefault(
      "terminatormaxtimepostrupture", std::numeric_limits<double>::infinity());
  const auto terminatorMomentRateThreshold = abortCriteriaReader->readWithDefault(
      "terminatormomentratethreshold", -std::numeric_limits<double>::infinity());

  return EnergyOutputParameters{enabled,
                                computeVolumeEnergiesEveryOutput,
                                interval,
                                terminalOutput,
                                terminalPrecision,
                                terminatorMaxTimePostRupture,
                                terminatorMomentRateThreshold};
}

FreeSurfaceOutputParameters readFreeSurfaceParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");

  auto enabled = reader->readWithDefault("surfaceoutput", false);
  const auto interval = reader->readWithDefault("surfaceoutputinterval", VeryLongTime);
  warnIntervalAndDisable(enabled, interval, "surfaceoutput", "surfaceoutputinterval");

  const auto refinement = reader->readWithDefault("surfaceoutputrefinement", 0U);

  const auto vtkorder = reader->readWithDefault("surfacevtkorder", -1);

  return FreeSurfaceOutputParameters{enabled, refinement, interval, vtkorder};
}

PickpointParameters readPickpointParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("pickpoint");

  const auto printTimeInterval = reader->readWithDefault("printtimeinterval", 1);
  const auto maxPickStore = reader->readWithDefault("maxpickstore", 50);

  const auto outputMaskString =
      reader->readWithDefault<std::string>("outputmask", "1 1 1 1 1 1 0 0 0 0 0 0");
  const std::array<bool, 12> outputMask = convertStringToArray<bool, 12>(outputMaskString, false);

  const auto pickpointFileName = reader->readWithDefault("ppfilename", std::string(""));

  const auto collectiveio = reader->readWithDefault("receivercollectiveio", false);

  reader->warnDeprecated({"noutpoints"});

  return PickpointParameters{
      printTimeInterval, maxPickStore, outputMask, pickpointFileName, collectiveio};
}

ReceiverOutputParameters readReceiverParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");

  const auto interval = reader->readWithDefault("receiveroutputinterval", VeryLongTime);
  auto enabled = reader->readWithDefault("receiveroutput", true);
  warnIntervalAndDisable(enabled, interval, "receiveroutput", "receiveroutputinterval");

  const auto computeRotation = reader->readWithDefault("receivercomputerotation", false);
  const auto computeStrain = reader->readWithDefault("receivercomputestrain", false);
  const auto samplingInterval = reader->readWithDefault("pickdt", 0.005);
  const auto fileName = reader->readWithDefault("rfilename", std::string(""));

  warnIntervalAndDisable(enabled, samplingInterval, "receiveroutput", "pickdt");

  const auto collectiveio = reader->readWithDefault("receivercollectiveio", false);

  return ReceiverOutputParameters{
      enabled, computeRotation, computeStrain, interval, samplingInterval, fileName, collectiveio};
}

WaveFieldOutputParameters readWaveFieldParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");

  auto enabled = reader->readWithDefault("wavefieldoutput", true);
  const auto interval = reader->readWithDefault("timeinterval", VeryLongTime);
  warnIntervalAndDisable(enabled, interval, "wavefieldoutput", "timeinterval");
  const auto refinement =
      reader->readWithDefaultEnum<VolumeRefinement>("refinement",
                                                    VolumeRefinement::NoRefine,
                                                    {VolumeRefinement::NoRefine,
                                                     VolumeRefinement::Refine4,
                                                     VolumeRefinement::Refine8,
                                                     VolumeRefinement::Refine32});

  const auto boundsString =
      reader->readWithDefault("outputregionbounds", std::string("0.0 0.0 0.0 0.0 0.0 0.0"));
  const auto boundsRaw = convertStringToArray<double, 6>(boundsString);
  const auto boundsEnabled =
      std::none_of(boundsRaw.begin(), boundsRaw.end(), [](double d) { return d == 0; });
  const OutputInterval intervalX = {boundsRaw[0], boundsRaw[1]};
  const OutputInterval intervalY = {boundsRaw[2], boundsRaw[3]};
  const OutputInterval intervalZ = {boundsRaw[4], boundsRaw[5]};
  const OutputBounds bounds(boundsEnabled, intervalX, intervalY, intervalZ);

  const auto format = reader->readWithDefaultEnum<OutputFormat>(
      "format", OutputFormat::None, {OutputFormat::None, OutputFormat::Xdmf});
  if (enabled && format == OutputFormat::None) {
    logInfo(seissol::MPI::mpi.rank())
        << "Disabling the wavefield output by setting \"outputformat = 10\" is deprecated "
           "and may be removed in a future version of SeisSol. Consider using the parameter "
           "\"wavefieldoutput\" instead. To disable wavefield output, add \"wavefieldoutput "
           "= 0\" to the \"output\" section of your parameters file.";

    enabled = false;
  }

  const auto outputMaskString =
      reader->readOrFail<std::string>("ioutputmask", "No output mask given.");
  const std::array<bool, seissol::model::MaterialT::NumQuantities> outputMask =
      convertStringToArray<bool, seissol::model::MaterialT::NumQuantities>(outputMaskString, false);

  const auto plasticityMaskString =
      reader->readWithDefault("iplasticitymask", std::string("0 0 0 0 0 0 1"));
  const std::array<bool, 7> plasticityMask =
      convertStringToArray<bool, 7>(plasticityMaskString, false);

  const auto integrationMaskString =
      reader->readWithDefault("integrationmask", std::string("0 0 0 0 0 0 0 0 0"));
  const std::array<bool, 9> integrationMask =
      convertStringToArray<bool, 9>(integrationMaskString, false);

  const auto groupsRaw = reader->readWithDefault("outputgroups", std::vector<int>());
  const auto groups = std::unordered_set<int>(groupsRaw.begin(), groupsRaw.end());

  const auto vtkorder = reader->readWithDefault("wavefieldvtkorder", -1);

  return WaveFieldOutputParameters{enabled,
                                   vtkorder,
                                   interval,
                                   refinement,
                                   bounds,
                                   outputMask,
                                   plasticityMask,
                                   integrationMask,
                                   groups};
}

OutputParameters readOutputParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");

  const auto loopStatisticsNetcdfOutput =
      reader->readWithDefault("loopstatisticsnetcdfoutput", false);
  const auto format = reader->readWithDefaultEnum<OutputFormat>(
      "format", OutputFormat::None, {OutputFormat::None, OutputFormat::Xdmf});
  const auto xdmfWriterBackend = reader->readWithDefaultStringEnum<xdmfwriter::BackendType>(
      "xdmfwriterbackend",
      "posix",
      {
          {"posix", xdmfwriter::BackendType::POSIX},
#ifdef USE_HDF
          {"hdf5", xdmfwriter::BackendType::H5},
#endif
      });
  const auto prefix =
      reader->readOrFail<std::string>("outputfile", "Output file prefix not defined.");

  const auto checkpointParameters = readCheckpointParameters(baseReader);
  const auto elementwiseParameters = readElementwiseParameters(baseReader);
  const auto energyParameters = readEnergyParameters(baseReader);
  const auto freeSurfaceParameters = readFreeSurfaceParameters(baseReader);
  const auto pickpointParameters = readPickpointParameters(baseReader);
  const auto receiverParameters = readReceiverParameters(baseReader);
  const auto waveFieldParameters = readWaveFieldParameters(baseReader);

  reader->warnDeprecated({"rotation",
                          "interval",
                          "nrecordpoints",
                          "printintervalcriterion",
                          "pickdttype",
                          "ioutputmaskmaterial",
                          "faultoutputflag"});

  return OutputParameters(loopStatisticsNetcdfOutput,
                          format,
                          xdmfWriterBackend,
                          prefix,
                          checkpointParameters,
                          elementwiseParameters,
                          energyParameters,
                          freeSurfaceParameters,
                          pickpointParameters,
                          receiverParameters,
                          waveFieldParameters);
}
} // namespace seissol::initializer::parameters
