// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "OutputParameters.h"

#include "Core/Common/ConfigLayout.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/InputAux.h"
#include "Core/IO/WriterParameters.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/Plasticity/Physics/PlasticityQuantities.h"
#include "Simulation/Parameters/ParameterReader.h"
#include "WavePropagation/Output/ReceiverParameters.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::parameters {

namespace {
/**
 * A mask over the quantities of the material of every configuration in `configs`, by its id: the
 * first entries of `mask`, as many as the material has quantities; those `mask` does not give are
 * off.
 */
std::vector<std::vector<bool>> quantityMasks(const std::string& mask,
                                             const std::vector<ConfigId>& configs) {
  std::vector<std::vector<bool>> masks(builtConfigCount());
  for (const auto config : configs) {
    masks[config] = convertStringToVector<bool>(mask, configLayout(config).numQuantities, false);
  }
  return masks;
}

ProjectionMethod readProjectionMethod(ParameterReader* reader,
                                      const std::string& field,
                                      ProjectionMethod defaultValue) {
  return reader->readWithDefaultStringEnum<ProjectionMethod>(
      field,
      defaultValue,
      {
          {"pointwise", ProjectionMethod::Pointwise},
          {"l2", ProjectionMethod::L2},
      });
}

// the mode that `field` chooses, or else `fallback`
std::optional<TimeSeriesMode> readTimeSeriesMode(ParameterReader* reader,
                                                 const std::string& field,
                                                 std::optional<TimeSeriesMode> fallback) {
  if (!reader->hasField(field)) {
    return fallback;
  }
  return reader->readWithDefaultStringEnum<TimeSeriesMode>(
      field,
      TimeSeriesMode::Snapshot,
      {
          {"snapshot", TimeSeriesMode::Snapshot},
          {"incremental", TimeSeriesMode::Incremental},
          {"monolith", TimeSeriesMode::Monolith},
      });
}
// the mask of the twelve quantities of a fault output, from its text
std::array<bool, 12> faultOutputMask(const std::string& text) {
  return convertStringToArray<bool, 12>(text, false);
}
} // namespace

void warnIntervalAndDisable(bool& enabled,
                            double interval,
                            const std::string& valName,
                            const std::string& intName) {
  if (enabled && interval <= 0) {
    const auto intPhrase = valName + " = 0";
    logInfo() << "In your parameter file, you have specified a non-positive interval for" << intName
              << ". This mechanism is deprecated and may be removed in a future version of "
                 "SeisSol. Consider disabling the whole module by setting"
              << valName << "to 0 instead by adding" << intPhrase
              << "to the \"output\" section of your parameter file instead.";
    // still, replicate the old behavior.
    enabled = false;
  }
}

std::optional<CheckpointParameters> readCheckpointParameters(ParameterReader* baseReader) {
  const CheckpointParameters defaults{};
  auto* reader = baseReader->readSubNode("output");

  // a run writes checkpoints if it says how often, unless it switches them off
  auto enabled = reader->readWithDefault("checkpoint", reader->hasField("checkpointinterval"));
  double interval = defaults.interval;
  if (enabled) {
    interval = reader->readWithDefault("checkpointinterval", defaults.interval);
    warnIntervalAndDisable(enabled, interval, "checkpoint", "checkpointinterval");
  } else {
    reader->markUnused({"checkpointinterval"});
  }

  reader->warnDeprecated({"checkpointbackend", "checkpointfile"});

  if (!enabled) {
    return std::nullopt;
  }
  return CheckpointParameters{.interval = interval};
}

ElementwiseFaultParameters
    readElementwiseParameters(ParameterReader* baseReader,
                              std::optional<TimeSeriesMode> defaultTimeSeries) {
  const ElementwiseFaultParameters defaults{};
  auto* reader = baseReader->readSubNode("elementwise");

  const auto printTimeIntervalSec =
      reader->readWithDefault("printtimeinterval_sec", defaults.printTimeIntervalSec);
  const auto outputMask =
      reader->readWithDefaultConverted("outputmask", defaults.outputMask, faultOutputMask);
  const auto refinementStrategy = reader->readWithDefaultEnum<FaultRefinement>(
      "refinement_strategy",
      defaults.refinementStrategy,
      {FaultRefinement::Triple, FaultRefinement::Quad, FaultRefinement::None});
  const int refinement = reader->readWithDefault("refinement", defaults.refinement);
  reader->warnDeprecated({"printintervalcriterion"});

  const auto vtkorder = reader->readWithDefault("vtkorder", defaults.vtkorder);

  const auto timeSeries = readTimeSeriesMode(reader, "timeseries", defaultTimeSeries);

  return ElementwiseFaultParameters{.printTimeIntervalSec = printTimeIntervalSec,
                                    .outputMask = outputMask,
                                    .refinementStrategy = refinementStrategy,
                                    .refinement = refinement,
                                    .vtkorder = vtkorder,
                                    .timeSeries = timeSeries};
}

std::optional<EnergyOutputParameters> readEnergyParameters(ParameterReader* baseReader) {
  const EnergyOutputParameters defaults{};
  auto* reader = baseReader->readSubNode("output");

  bool enabled = reader->readWithDefault("energyoutput", false);
  const auto interval = reader->readWithDefault("energyoutputinterval", defaults.interval);
  warnIntervalAndDisable(enabled, interval, "energyoutput", "energyoutputinterval");

  const auto computeVolumeEnergiesEveryOutput = reader->readWithDefault(
      "computevolumeenergieseveryoutput", defaults.computeVolumeEnergiesEveryOutput);
  const auto terminalOutput =
      reader->readWithDefault("energyterminaloutput", defaults.terminalOutput);
  const auto terminalPrecision =
      reader->readWithDefault("energyterminalprecision", defaults.terminalPrecision);

  auto* abortCriteriaReader = baseReader->readSubNode("abortcriteria");
  const auto terminatorMaxTimePostRupture = abortCriteriaReader->readWithDefault(
      "terminatormaxtimepostrupture", defaults.terminatorMaxTimePostRupture);
  const auto terminatorMomentRateThreshold = abortCriteriaReader->readWithDefault(
      "terminatormomentratethreshold", defaults.terminatorMomentRateThreshold);

  if (!enabled) {
    return std::nullopt;
  }
  return EnergyOutputParameters{.computeVolumeEnergiesEveryOutput =
                                    computeVolumeEnergiesEveryOutput,
                                .interval = interval,
                                .terminalOutput = terminalOutput,
                                .terminalPrecision = terminalPrecision,
                                .terminatorMaxTimePostRupture = terminatorMaxTimePostRupture,
                                .terminatorMomentRateThreshold = terminatorMomentRateThreshold};
}

std::optional<FreeSurfaceOutputParameters>
    readFreeSurfaceParameters(ParameterReader* baseReader,
                              std::optional<TimeSeriesMode> defaultTimeSeries,
                              const std::vector<ConfigId>& configs) {
  const FreeSurfaceOutputParameters defaults{};
  auto* reader = baseReader->readSubNode("output");

  auto enabled = reader->readWithDefault("surfaceoutput", false);
  const auto interval = reader->readWithDefault("surfaceoutputinterval", defaults.interval);
  warnIntervalAndDisable(enabled, interval, "surfaceoutput", "surfaceoutputinterval");

  const auto refinement = reader->readWithDefault("surfaceoutputrefinement", defaults.refinement);

  // Every level splits each output triangle into four, and the whole refined surface is held in
  // memory while the output is assembled. There is no level at which that stops working, so this
  // states the price rather than refusing to pay it.
  constexpr unsigned RefinementWarningDepth = 3;
  if (refinement > RefinementWarningDepth) {
    logWarning() << "surfaceoutputrefinement =" << refinement << "puts" << (1U << (2 * refinement))
                 << "output cells on every surface face. Memory and output size grow by a factor "
                    "of four with every level.";
  }

  const auto vtkorder = reader->readWithDefault("surfacevtkorder", defaults.vtkorder);

  // the velocities by default, wherever the material has them (e.g. acoustic or elastic)
  const auto surfaceOutputMaskString = reader->read<std::string>("surfaceoutputmask");
  std::vector<std::vector<bool>> surfaceOutputMask(builtConfigCount());
  for (const auto config : configs) {
    const auto& layout = configLayout(config);
    std::string velocities;
    for (std::size_t i = 0; i < layout.velocityOffset; ++i) {
      velocities += "0 ";
    }
    velocities += "1 1 1";
    surfaceOutputMask[config] = convertStringToVector<bool>(
        surfaceOutputMaskString.value_or(velocities), layout.numQuantities, false);
  }

  // The free-surface output has always been an average over each output subcell (cf. the former
  // FreeSurfaceIntegrator::computeSubTriangleAverages), i.e. an L2 projection.
  const auto projection = readProjectionMethod(reader, "surfaceprojection", defaults.projection);

  const auto timeSeries = readTimeSeriesMode(reader, "surfacetimeseries", defaultTimeSeries);

  if (!enabled) {
    return std::nullopt;
  }
  return FreeSurfaceOutputParameters{.refinement = refinement,
                                     .interval = interval,
                                     .vtkorder = vtkorder,
                                     .outputMask = surfaceOutputMask,
                                     .projection = projection,
                                     .timeSeries = timeSeries};
}

PickpointParameters readPickpointParameters(ParameterReader* baseReader) {
  const PickpointParameters defaults{};
  auto* reader = baseReader->readSubNode("pickpoint");

  const auto printTimeInterval =
      reader->readWithDefault("printtimeinterval", defaults.printTimeInterval);
  if (printTimeInterval < 1) {
    logError() << "printtimeinterval of the on-fault receivers has to be at least 1, but is"
               << printTimeInterval;
  }

  // zero: every printtimeinterval time steps of the finest time cluster with fault faces
  const auto samplingInterval = reader->readWithDefault("pickdt", defaults.samplingInterval);
  if (samplingInterval < 0) {
    logError() << "pickdt of the on-fault receivers must not be negative, but is"
               << samplingInterval;
  }

  const auto interval = reader->readWithDefault("outputinterval", defaults.writeInterval);

  const auto outputMask =
      reader->readWithDefaultConverted("outputmask", defaults.outputMask, faultOutputMask);

  const auto pickpointFileName = reader->readPath("ppfilename");

  const auto collectiveio = reader->readWithDefault("receivercollectiveio", defaults.collectiveio);
  const auto format = reader->readWithDefaultStringEnum<ReceiverOutputFormat>(
      "format",
      defaults.format,
      {
          {"csv", ReceiverOutputFormat::Csv},
          {"hdf5", ReceiverOutputFormat::Hdf5},
      });
  const auto samplechunk = reader->readWithDefault("samplechunk", defaults.samplechunk);
  if (reader->readWithDefault("aggregateperrank", false)) {
    logWarning() << "aggregateperrank is no longer supported: the text output writes a file per "
                    "on-fault receiver, and format = 'hdf5' writes all of them into a single file.";
  }

  if (collectiveio && format != ReceiverOutputFormat::Hdf5) {
    logError() << "Collective IO for the on-fault receiver output is still under construction.";
  }

  reader->warnDeprecated({"noutpoints", "maxpickstore"});

  return PickpointParameters{.printTimeInterval = printTimeInterval,
                             .writeInterval = interval,
                             .outputMask = outputMask,
                             .pickpointFileName = pickpointFileName,
                             .collectiveio = collectiveio,
                             .format = format,
                             .samplechunk = samplechunk,
                             .samplingInterval = samplingInterval};
}

std::optional<ReceiverOutputParameters> readReceiverParameters(ParameterReader* baseReader) {
  const ReceiverOutputParameters defaults{};
  auto* reader = baseReader->readSubNode("output");

  const auto interval = reader->readWithDefault("receiveroutputinterval", defaults.writeInterval);
  // a run writes the off-fault receivers if it gives their points, unless it switches them off
  auto enabled = reader->readWithDefault("receiveroutput", reader->hasField("rfilename"));
  warnIntervalAndDisable(enabled, interval, "receiveroutput", "receiveroutputinterval");

  const auto format = reader->readWithDefaultStringEnum<ReceiverOutputFormat>(
      "receiverformat",
      defaults.format,
      {
          {"csv", ReceiverOutputFormat::Csv},
          {"hdf5", ReceiverOutputFormat::Hdf5},
      });

  const auto computeRotation =
      reader->readWithDefault("receivercomputerotation", defaults.computeRotation);
  const auto computeStrain =
      reader->readWithDefault("receivercomputestrainrate", defaults.computeStrain);
  const auto samplingInterval = reader->readWithDefault("pickdt", defaults.samplingInterval);
  const auto fileName = reader->readPath("rfilename");

  warnIntervalAndDisable(enabled, samplingInterval, "receiveroutput", "pickdt");

  const auto collectiveio = reader->readWithDefault("receivercollectiveio", defaults.collectiveio);
  const auto samplechunk = reader->readWithDefault("receiversamplechunk", defaults.samplechunk);

  if (enabled && !fileName.has_value()) {
    logError() << "The off-fault receiver output is enabled, but no receiver point file was given.";
  }

  if (!enabled) {
    return std::nullopt;
  }
  return ReceiverOutputParameters{.format = format,
                                  .computeRotation = computeRotation,
                                  .computeStrain = computeStrain,
                                  .writeInterval = interval,
                                  .samplingInterval = samplingInterval,
                                  .fileName = fileName.value_or(defaults.fileName),
                                  .collectiveio = collectiveio,
                                  .samplechunk = samplechunk};
}

std::optional<WaveFieldOutputParameters>
    readWaveFieldParameters(ParameterReader* baseReader,
                            std::optional<TimeSeriesMode> defaultTimeSeries,
                            const std::vector<ConfigId>& configs) {
  const WaveFieldOutputParameters defaults{};
  auto* reader = baseReader->readSubNode("output");

  bool enabled = false;
  const auto enabledPre = reader->read<bool>("wavefieldoutput");
  if (enabledPre.has_value()) {
    enabled = enabledPre.value();
  } else if (reader->hasField("format")) {
    const auto format = reader->readWithDefaultEnum<OutputFormat>(
        "format", OutputFormat::None, {OutputFormat::None, OutputFormat::Xdmf});

    // TODO: deprecate the "format" value at some point (not yet though)
    logInfo()
        << "Disabling/enabling the wavefield output via the \"format\" option is deprecated "
           "and may be removed in a future version of SeisSol. Consider using the parameter "
           "\"wavefieldoutput\" instead. To enable/disable wavefield output, add \"wavefieldoutput "
           "= 1\" or \"wavefieldoutput = 0\" to the \"output\" section of your parameters file, "
           "respectively.";

    enabled = format == OutputFormat::Xdmf;
  }

  const auto interval = reader->readWithDefault("timeinterval", defaults.interval);
  warnIntervalAndDisable(enabled, interval, "wavefieldoutput", "timeinterval");
  const auto refinement =
      reader->readWithDefaultEnum<VolumeRefinement>("refinement",
                                                    defaults.refinement,
                                                    {VolumeRefinement::NoRefine,
                                                     VolumeRefinement::Refine4,
                                                     VolumeRefinement::Refine8,
                                                     VolumeRefinement::Refine32});

  // the region is on only if none of its bounds is zero
  const auto bounds = reader->readWithDefaultConverted(
      "outputregionbounds", defaults.bounds, [](const std::string& text) {
        const auto boundsRaw = convertStringToArray<double, 6>(text);
        const auto boundsEnabled = std::ranges::none_of(boundsRaw, [](double d) { return d == 0; });
        const OutputInterval intervalX = {.lower = boundsRaw[0], .upper = boundsRaw[1]};
        const OutputInterval intervalY = {.lower = boundsRaw[2], .upper = boundsRaw[3]};
        const OutputInterval intervalZ = {.lower = boundsRaw[4], .upper = boundsRaw[5]};
        return OutputBounds(boundsEnabled, intervalX, intervalY, intervalZ);
      });

  const auto outputMaskString =
      reader->readOrFail<std::string>("ioutputmask", "No output mask given.");
  const auto outputMask = quantityMasks(outputMaskString, configs);

  const auto plasticityMask = reader->readWithDefaultConverted(
      "iplasticitymask", defaults.plasticityMask, [](const std::string& text) {
        return convertStringToArray<bool, seissol::model::PlasticityQuantityCount>(text, false);
      });

  const auto integrationMaskString =
      reader->readWithDefault("integrationmask", std::string("0 0 0 0 0 0 0 0 0"));
  const auto integrationMask = quantityMasks(integrationMaskString, configs);

  const auto groupsRaw = reader->readWithDefault(
      "outputgroups", std::vector<int>(defaults.groups.begin(), defaults.groups.end()));
  const auto groups = std::unordered_set<int>(groupsRaw.begin(), groupsRaw.end());

  const auto vtkorder = reader->readWithDefault("wavefieldvtkorder", defaults.vtkorder);

  const auto computeRotation =
      reader->readWithDefault("wavefieldcomputerotation", defaults.computeRotation);
  const auto computeStrain =
      reader->readWithDefault("wavefieldcomputestrain", defaults.computeStrain);

  // The wavefield output has always been a point evaluation at the output points (cf. the former
  // refinement::VariableSubsampler), unlike the free-surface output.
  const auto projection = readProjectionMethod(reader, "wavefieldprojection", defaults.projection);

  const auto timeSeries = readTimeSeriesMode(reader, "wavefieldtimeseries", defaultTimeSeries);

  if (enabledPre.has_value()) {
    reader->warnDeprecated({"format"});
  }

  if (!enabled) {
    return std::nullopt;
  }
  return WaveFieldOutputParameters{.vtkorder = vtkorder,
                                   .interval = interval,
                                   .refinement = refinement,
                                   .bounds = bounds,
                                   .outputMask = outputMask,
                                   .plasticityMask = plasticityMask,
                                   .integrationMask = integrationMask,
                                   .groups = groups,
                                   .computeRotation = computeRotation,
                                   .computeStrain = computeStrain,
                                   .projection = projection,
                                   .timeSeries = timeSeries};
}

OutputParameters readOutputParameters(ParameterReader* baseReader,
                                      const std::vector<ConfigId>& configs,
                                      bool withFault) {
  const OutputParameters defaults{};
  auto* reader = baseReader->readSubNode("output");

  const auto hdfcompress = reader->readWithDefault("hdfcompress", defaults.hdfcompress);

  // Each output has a field of its own that falls back to this one; where neither is given, the
  // writer chooses by the format of the output (see TimeSeriesMode).
  const auto defaultTimeSeries = readTimeSeriesMode(reader, "outputtimeseries", std::nullopt);
  const auto xdmfWriterBackend =
      reader->readWithDefaultStringEnum<XdmfBackend>("xdmfwriterbackend",
                                                     defaults.xdmfWriterBackend,
                                                     {
                                                         {"posix", XdmfBackend::Posix},
#ifdef USE_HDF
                                                         {"hdf5", XdmfBackend::Hdf5},
#endif
                                                     });

  // the directory and the name of the outputs: outdir and outname, or the older outputfile, which
  // gives both at once; each of the former takes the place of its part of the latter (e.g. outdir
  // given on the command line, for a parameter file with outputfile)
  const auto outputDirR = reader->read<std::string>("outdir");
  const auto outputNameR = reader->read<std::string>("outname");
  std::string outputDir;
  std::string outputName;
  if (outputDirR.has_value() && outputNameR.has_value()) {
    reader->warnDeprecatedSingle("outputfile");
  } else if (const auto prefix = reader->read<std::string>("outputfile"); prefix.has_value()) {
    const auto prefixPath = std::filesystem::path(prefix.value());
    outputDir = prefixPath.parent_path().string();
    outputName = prefixPath.filename().string();
  }
  outputDir = outputDirR.value_or(outputDir);
  outputName = outputNameR.value_or(outputName);
  if (outputName.empty()) {
    logError() << "The outputs have no name: give outname (and outdir), or outputfile.";
  }

  const auto dirPerSim = reader->readWithDefault("outsimdir", defaults.outputDirPerSim);
  const auto epochInterval = reader->readWithDefault("epochinterval", defaults.epochInterval);
  const auto vertexWeightOutput =
      reader->readWithDefault("vertexweightoutput", defaults.vertexWeightOutput);

  const auto checkpoint = readCheckpointParameters(baseReader);
  if (checkpoint.has_value() && configs.size() > 1) {
    logError() << "Checkpoints are not supported yet for runs whose cells compute in several "
                  "configurations.";
  }
  const auto energy = readEnergyParameters(baseReader);

  // the outputs, each named as the files it has always written
  std::vector<OutputEntry> outputs;
  if (auto waveField = readWaveFieldParameters(baseReader, defaultTimeSeries, configs)) {
    outputs.push_back({.name = "wavefield", .settings = std::move(waveField.value())});
  }
  if (auto surface = readFreeSurfaceParameters(baseReader, defaultTimeSeries, configs)) {
    outputs.push_back({.name = "surface", .settings = std::move(surface.value())});
  }
  if (auto receivers = readReceiverParameters(baseReader)) {
    outputs.push_back({.name = "receivers", .settings = std::move(receivers.value())});
  }

  // the outputs of the fault, which the section of the rupture switches on
  const auto faultOutputs =
      baseReader->readSubNode("dynamicrupture")
          ->readWithDefaultEnum<OutputType>("outputpointtype",
                                            OutputType::None,
                                            {OutputType::None,
                                             OutputType::AtPickpoint,
                                             OutputType::Elementwise,
                                             OutputType::AtPickpointAndElementwise});
  const auto elementwise = readElementwiseParameters(baseReader, defaultTimeSeries);
  auto pickpoints = readPickpointParameters(baseReader);
  if (withFault && (faultOutputs == OutputType::Elementwise ||
                    faultOutputs == OutputType::AtPickpointAndElementwise)) {
    outputs.push_back({.name = "fault", .settings = elementwise});
  }
  if (withFault && (faultOutputs == OutputType::AtPickpoint ||
                    faultOutputs == OutputType::AtPickpointAndElementwise)) {
    outputs.push_back({.name = "faultreceivers", .settings = std::move(pickpoints)});
  }

  reader->warnDeprecated({"projection",
                          "rotation",
                          "interval",
                          "nrecordpoints",
                          "printintervalcriterion",
                          "pickdttype",
                          "ioutputmaskmaterial",
                          "faultoutputflag",
                          "loopstatisticsnetcdfoutput"});

  return OutputParameters{.xdmfWriterBackend = xdmfWriterBackend,
                          .hdfcompress = hdfcompress,
                          .outputPath = outputDir,
                          .outputName = outputName,
                          .outputDirPerSim = dirPerSim,
                          .checkpoint = checkpoint,
                          .energy = energy,
                          .outputs = std::move(outputs),
                          .epochInterval = epochInterval,
                          .vertexWeightOutput = vertexWeightOutput};
}
} // namespace seissol::initializer::parameters
