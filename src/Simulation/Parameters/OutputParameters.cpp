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
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <utils/logger.h>
#include <variant>
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

// The names of the fields of the outputs of each kind: in a block of their own
// (&Output<Kind><n>), or as the parameter file has always given them, in &Output, &Elementwise
// and &Pickpoint.

struct WaveFieldFields {
  const char* interval;
  const char* refinement;
  const char* bounds;
  const char* mask;
  const char* plasticityMask;
  const char* integrationMask;
  const char* groups;
  const char* vtkorder;
  const char* computeRotation;
  const char* computeStrain;
  const char* projection;
  const char* timeSeries;
};

constexpr WaveFieldFields LegacyWaveFieldFields{.interval = "timeinterval",
                                                .refinement = "refinement",
                                                .bounds = "outputregionbounds",
                                                .mask = "ioutputmask",
                                                .plasticityMask = "iplasticitymask",
                                                .integrationMask = "integrationmask",
                                                .groups = "outputgroups",
                                                .vtkorder = "wavefieldvtkorder",
                                                .computeRotation = "wavefieldcomputerotation",
                                                .computeStrain = "wavefieldcomputestrain",
                                                .projection = "wavefieldprojection",
                                                .timeSeries = "wavefieldtimeseries"};

constexpr WaveFieldFields BlockWaveFieldFields{.interval = "interval",
                                               .refinement = "refinement",
                                               .bounds = "bounds",
                                               .mask = "mask",
                                               .plasticityMask = "plasticitymask",
                                               .integrationMask = "integrationmask",
                                               .groups = "groups",
                                               .vtkorder = "vtkorder",
                                               .computeRotation = "computerotation",
                                               .computeStrain = "computestrain",
                                               .projection = "projection",
                                               .timeSeries = "timeseries"};

struct SurfaceFields {
  const char* interval;
  const char* refinement;
  const char* vtkorder;
  const char* mask;
  const char* projection;
  const char* timeSeries;
};

constexpr SurfaceFields LegacySurfaceFields{.interval = "surfaceoutputinterval",
                                            .refinement = "surfaceoutputrefinement",
                                            .vtkorder = "surfacevtkorder",
                                            .mask = "surfaceoutputmask",
                                            .projection = "surfaceprojection",
                                            .timeSeries = "surfacetimeseries"};

constexpr SurfaceFields BlockSurfaceFields{.interval = "interval",
                                           .refinement = "refinement",
                                           .vtkorder = "vtkorder",
                                           .mask = "mask",
                                           .projection = "projection",
                                           .timeSeries = "timeseries"};

struct ReceiverFields {
  const char* file;
  const char* interval;
  const char* samplingInterval;
  const char* format;
  const char* computeRotation;
  const char* computeStrain;
  const char* collectiveio;
  const char* samplechunk;
};

constexpr ReceiverFields LegacyReceiverFields{.file = "rfilename",
                                              .interval = "receiveroutputinterval",
                                              .samplingInterval = "pickdt",
                                              .format = "receiverformat",
                                              .computeRotation = "receivercomputerotation",
                                              .computeStrain = "receivercomputestrainrate",
                                              .collectiveio = "receivercollectiveio",
                                              .samplechunk = "receiversamplechunk"};

constexpr ReceiverFields BlockReceiverFields{.file = "file",
                                             .interval = "interval",
                                             .samplingInterval = "samplinginterval",
                                             .format = "format",
                                             .computeRotation = "computerotation",
                                             .computeStrain = "computestrainrate",
                                             .collectiveio = "collectiveio",
                                             .samplechunk = "samplechunk"};

struct FaultFields {
  const char* interval;
  const char* mask;
  const char* refinementStrategy;
  const char* refinement;
  const char* vtkorder;
  const char* timeSeries;
};

constexpr FaultFields LegacyFaultFields{.interval = "printtimeinterval_sec",
                                        .mask = "outputmask",
                                        .refinementStrategy = "refinement_strategy",
                                        .refinement = "refinement",
                                        .vtkorder = "vtkorder",
                                        .timeSeries = "timeseries"};

constexpr FaultFields BlockFaultFields{.interval = "interval",
                                       .mask = "mask",
                                       .refinementStrategy = "refinementstrategy",
                                       .refinement = "refinement",
                                       .vtkorder = "vtkorder",
                                       .timeSeries = "timeseries"};

struct FaultReceiverFields {
  const char* file;
  const char* interval;
  const char* samplingInterval;
  const char* samplingSteps;
  const char* mask;
  const char* format;
  const char* collectiveio;
  const char* samplechunk;
};

constexpr FaultReceiverFields LegacyFaultReceiverFields{.file = "ppfilename",
                                                        .interval = "outputinterval",
                                                        .samplingInterval = "pickdt",
                                                        .samplingSteps = "printtimeinterval",
                                                        .mask = "outputmask",
                                                        .format = "format",
                                                        .collectiveio = "receivercollectiveio",
                                                        .samplechunk = "samplechunk"};

constexpr FaultReceiverFields BlockFaultReceiverFields{.file = "file",
                                                       .interval = "interval",
                                                       .samplingInterval = "samplinginterval",
                                                       .samplingSteps = "samplingsteps",
                                                       .mask = "mask",
                                                       .format = "format",
                                                       .collectiveio = "collectiveio",
                                                       .samplechunk = "samplechunk"};

ReceiverOutputFormat readReceiverFormat(ParameterReader* reader,
                                        const char* field,
                                        ReceiverOutputFormat defaultValue) {
  return reader->readWithDefaultStringEnum<ReceiverOutputFormat>(
      field,
      defaultValue,
      {
          {"csv", ReceiverOutputFormat::Csv},
          {"hdf5", ReceiverOutputFormat::Hdf5},
      });
}

/// The settings of a wave field output; all but whether it is on.
WaveFieldOutputParameters readWaveField(ParameterReader* reader,
                                        const WaveFieldFields& fields,
                                        std::optional<TimeSeriesMode> defaultTimeSeries,
                                        const std::vector<ConfigId>& configs) {
  const WaveFieldOutputParameters defaults{};
  const auto interval = reader->readWithDefault(fields.interval, defaults.interval);
  const auto refinement =
      reader->readWithDefaultEnum<VolumeRefinement>(fields.refinement,
                                                    defaults.refinement,
                                                    {VolumeRefinement::NoRefine,
                                                     VolumeRefinement::Refine4,
                                                     VolumeRefinement::Refine8,
                                                     VolumeRefinement::Refine32});

  // the region is on only if none of its bounds is zero
  const auto bounds =
      reader->readWithDefaultConverted(fields.bounds, defaults.bounds, [](const std::string& text) {
        const auto boundsRaw = convertStringToArray<double, 6>(text);
        const auto boundsEnabled = std::ranges::none_of(boundsRaw, [](double d) { return d == 0; });
        const OutputInterval intervalX = {.lower = boundsRaw[0], .upper = boundsRaw[1]};
        const OutputInterval intervalY = {.lower = boundsRaw[2], .upper = boundsRaw[3]};
        const OutputInterval intervalZ = {.lower = boundsRaw[4], .upper = boundsRaw[5]};
        return OutputBounds(boundsEnabled, intervalX, intervalY, intervalZ);
      });

  const auto outputMaskString =
      reader->readOrFail<std::string>(fields.mask, "No output mask given.");
  const auto outputMask = quantityMasks(outputMaskString, configs);

  const auto plasticityMask = reader->readWithDefaultConverted(
      fields.plasticityMask, defaults.plasticityMask, [](const std::string& text) {
        return convertStringToArray<bool, seissol::model::PlasticityQuantityCount>(text, false);
      });

  const auto integrationMaskString =
      reader->readWithDefault(fields.integrationMask, std::string("0 0 0 0 0 0 0 0 0"));
  const auto integrationMask = quantityMasks(integrationMaskString, configs);

  const auto groupsRaw = reader->readWithDefault(
      fields.groups, std::vector<int>(defaults.groups.begin(), defaults.groups.end()));
  const auto groups = std::unordered_set<int>(groupsRaw.begin(), groupsRaw.end());

  const auto vtkorder = reader->readWithDefault(fields.vtkorder, defaults.vtkorder);

  const auto computeRotation =
      reader->readWithDefault(fields.computeRotation, defaults.computeRotation);
  const auto computeStrain = reader->readWithDefault(fields.computeStrain, defaults.computeStrain);

  // The wavefield output has always been a point evaluation at the output points (cf. the former
  // refinement::VariableSubsampler), unlike the free-surface output.
  const auto projection = readProjectionMethod(reader, fields.projection, defaults.projection);

  const auto timeSeries = readTimeSeriesMode(reader, fields.timeSeries, defaultTimeSeries);

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

/// The settings of a free surface output; all but whether it is on.
FreeSurfaceOutputParameters readSurface(ParameterReader* reader,
                                        const SurfaceFields& fields,
                                        std::optional<TimeSeriesMode> defaultTimeSeries,
                                        const std::vector<ConfigId>& configs) {
  const FreeSurfaceOutputParameters defaults{};
  const auto interval = reader->readWithDefault(fields.interval, defaults.interval);
  const auto refinement = reader->readWithDefault(fields.refinement, defaults.refinement);

  // Every level splits each output triangle into four, and the whole refined surface is held in
  // memory while the output is assembled. There is no level at which that stops working, so this
  // states the price rather than refusing to pay it.
  constexpr unsigned RefinementWarningDepth = 3;
  if (refinement > RefinementWarningDepth) {
    logWarning() << std::string(fields.refinement) + " =" << refinement << "puts"
                 << (1U << (2 * refinement))
                 << "output cells on every surface face. Memory and output size grow by a factor "
                    "of four with every level.";
  }

  const auto vtkorder = reader->readWithDefault(fields.vtkorder, defaults.vtkorder);

  // the velocities by default, wherever the material has them (e.g. acoustic or elastic)
  const auto surfaceOutputMaskString = reader->read<std::string>(fields.mask);
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
  const auto projection = readProjectionMethod(reader, fields.projection, defaults.projection);

  const auto timeSeries = readTimeSeriesMode(reader, fields.timeSeries, defaultTimeSeries);

  return FreeSurfaceOutputParameters{.refinement = refinement,
                                     .interval = interval,
                                     .vtkorder = vtkorder,
                                     .outputMask = surfaceOutputMask,
                                     .projection = projection,
                                     .timeSeries = timeSeries};
}

/// The settings of a receiver output; all but its file and whether it is on.
ReceiverOutputParameters readReceivers(ParameterReader* reader, const ReceiverFields& fields) {
  const ReceiverOutputParameters defaults{};
  const auto interval = reader->readWithDefault(fields.interval, defaults.writeInterval);
  const auto format = readReceiverFormat(reader, fields.format, defaults.format);
  const auto computeRotation =
      reader->readWithDefault(fields.computeRotation, defaults.computeRotation);
  const auto computeStrain = reader->readWithDefault(fields.computeStrain, defaults.computeStrain);
  const auto samplingInterval =
      reader->readWithDefault(fields.samplingInterval, defaults.samplingInterval);
  const auto collectiveio = reader->readWithDefault(fields.collectiveio, defaults.collectiveio);
  const auto samplechunk = reader->readWithDefault(fields.samplechunk, defaults.samplechunk);
  return ReceiverOutputParameters{.format = format,
                                  .computeRotation = computeRotation,
                                  .computeStrain = computeStrain,
                                  .writeInterval = interval,
                                  .samplingInterval = samplingInterval,
                                  .fileName = defaults.fileName,
                                  .collectiveio = collectiveio,
                                  .samplechunk = samplechunk};
}

/// The settings of an elementwise fault output.
ElementwiseFaultParameters readFault(ParameterReader* reader,
                                     const FaultFields& fields,
                                     std::optional<TimeSeriesMode> defaultTimeSeries) {
  const ElementwiseFaultParameters defaults{};
  const auto printTimeIntervalSec =
      reader->readWithDefault(fields.interval, defaults.printTimeIntervalSec);
  const auto outputMask =
      reader->readWithDefaultConverted(fields.mask, defaults.outputMask, faultOutputMask);
  const auto refinementStrategy = reader->readWithDefaultEnum<FaultRefinement>(
      fields.refinementStrategy,
      defaults.refinementStrategy,
      {FaultRefinement::Triple, FaultRefinement::Quad, FaultRefinement::None});
  const int refinement = reader->readWithDefault(fields.refinement, defaults.refinement);
  const auto vtkorder = reader->readWithDefault(fields.vtkorder, defaults.vtkorder);
  const auto timeSeries = readTimeSeriesMode(reader, fields.timeSeries, defaultTimeSeries);
  return ElementwiseFaultParameters{.printTimeIntervalSec = printTimeIntervalSec,
                                    .outputMask = outputMask,
                                    .refinementStrategy = refinementStrategy,
                                    .refinement = refinement,
                                    .vtkorder = vtkorder,
                                    .timeSeries = timeSeries};
}

/// The settings of an output of on-fault receivers.
PickpointParameters readFaultReceivers(ParameterReader* reader, const FaultReceiverFields& fields) {
  const PickpointParameters defaults{};
  const auto printTimeInterval =
      reader->readWithDefault(fields.samplingSteps, defaults.printTimeInterval);
  if (printTimeInterval < 1) {
    logError() << fields.samplingSteps << "of the on-fault receivers has to be at least 1, but is"
               << printTimeInterval;
  }

  // zero: every printtimeinterval time steps of the finest time cluster with fault faces
  const auto samplingInterval =
      reader->readWithDefault(fields.samplingInterval, defaults.samplingInterval);
  if (samplingInterval < 0) {
    logError() << fields.samplingInterval
               << "of the on-fault receivers must not be negative, but is" << samplingInterval;
  }

  const auto interval = reader->readWithDefault(fields.interval, defaults.writeInterval);
  const auto outputMask =
      reader->readWithDefaultConverted(fields.mask, defaults.outputMask, faultOutputMask);
  const auto pickpointFileName = reader->readPath(fields.file);
  const auto collectiveio = reader->readWithDefault(fields.collectiveio, defaults.collectiveio);
  const auto format = readReceiverFormat(reader, fields.format, defaults.format);
  const auto samplechunk = reader->readWithDefault(fields.samplechunk, defaults.samplechunk);

  if (collectiveio && format != ReceiverOutputFormat::Hdf5) {
    logError() << "Collective IO for the on-fault receiver output is still under construction.";
  }

  return PickpointParameters{.printTimeInterval = printTimeInterval,
                             .writeInterval = interval,
                             .outputMask = outputMask,
                             .pickpointFileName = pickpointFileName,
                             .collectiveio = collectiveio,
                             .format = format,
                             .samplechunk = samplechunk,
                             .samplingInterval = samplingInterval};
}

// The kinds of outputs that have blocks of their own (&Output<Kind><n>), in the order the list of
// outputs holds their blocks, with the names their blocks end in. A kind whose name begins that of
// another one comes after it.
enum class BlockKind { WaveField, Surface, Receivers, FaultReceivers, Fault };

struct BlockKindName {
  std::string_view name;
  BlockKind kind;
};

constexpr std::array<BlockKindName, 5> BlockKinds = {{
    {.name = "wavefield", .kind = BlockKind::WaveField},
    {.name = "surface", .kind = BlockKind::Surface},
    {.name = "receivers", .kind = BlockKind::Receivers},
    {.name = "faultreceivers", .kind = BlockKind::FaultReceivers},
    {.name = "fault", .kind = BlockKind::Fault},
}};

// a block of an output: its kind and number, and the name of its section
struct OutputBlock {
  BlockKind kind{};
  std::string kindName;
  std::string number;
  std::string section;
};

// the block of an output that the section `section` is, if it is one: output<kind><number>, the
// number made of digits, or none
std::optional<OutputBlock> outputBlockOf(const std::string& section) {
  constexpr std::string_view Prefix = "output";
  if (!section.starts_with(Prefix)) {
    return std::nullopt;
  }
  const auto rest = std::string_view(section).substr(Prefix.size());
  for (const auto& [name, kind] : BlockKinds) {
    if (!rest.starts_with(name)) {
      continue;
    }
    const auto number = rest.substr(name.size());
    if (std::ranges::all_of(number, [](char c) { return c >= '0' && c <= '9'; })) {
      return OutputBlock{.kind = kind,
                         .kindName = std::string(name),
                         .number = std::string(number),
                         .section = section};
    }
  }
  return std::nullopt;
}

// the number of the outputs of the kind T
template <typename T>
std::size_t countOf(const std::vector<OutputEntry>& outputs) {
  return static_cast<std::size_t>(std::ranges::count_if(
      outputs, [](const OutputEntry& entry) { return std::holds_alternative<T>(entry.settings); }));
}

/// The outputs in blocks of their own (&Output<Kind><n>), by kind and then by number; those of the
/// fault only with a fault (`withFault`).
std::vector<OutputEntry> readOutputBlocks(ParameterReader* baseReader,
                                          std::optional<TimeSeriesMode> defaultTimeSeries,
                                          const std::vector<ConfigId>& configs,
                                          bool withFault) {
  std::vector<OutputBlock> blocks;
  for (const auto& section : baseReader->fieldNames()) {
    if (auto block = outputBlockOf(section)) {
      blocks.push_back(std::move(block.value()));
    }
  }
  std::ranges::sort(blocks, [](const OutputBlock& a, const OutputBlock& b) {
    return std::tuple(a.kind, a.number.size(), a.number) <
           std::tuple(b.kind, b.number.size(), b.number);
  });

  std::vector<OutputEntry> outputs;
  for (const auto& block : blocks) {
    auto* reader = baseReader->readSubNode(block.section);
    const auto name = reader->readWithDefault("name", block.kindName + block.number);
    // a block switches its output on, so its times have to be positive
    const auto requirePositive = [&](const char* field, double value) {
      if (!(value > 0)) {
        logError() << "The field" << field << "of the output" << name
                   << "has to be positive, but is" << value;
      }
    };
    const auto requireFile = [&](const std::optional<std::string>& file) {
      if (!file.has_value()) {
        logError() << "The output" << name << "names no file of points (file).";
      }
      return file.value_or("");
    };
    // the outputs of the fault are read either way, so that their fields count as known
    const auto withoutFault = [&] {
      if (!withFault) {
        logWarning() << "The run has no fault, so it leaves out the output" << name << ".";
      }
      return !withFault;
    };
    switch (block.kind) {
    case BlockKind::WaveField: {
      auto settings = readWaveField(reader, BlockWaveFieldFields, defaultTimeSeries, configs);
      requirePositive(BlockWaveFieldFields.interval, settings.interval);
      outputs.push_back({.name = name, .settings = std::move(settings)});
      break;
    }
    case BlockKind::Surface: {
      auto settings = readSurface(reader, BlockSurfaceFields, defaultTimeSeries, configs);
      requirePositive(BlockSurfaceFields.interval, settings.interval);
      outputs.push_back({.name = name, .settings = std::move(settings)});
      break;
    }
    case BlockKind::Receivers: {
      auto settings = readReceivers(reader, BlockReceiverFields);
      requirePositive(BlockReceiverFields.interval, settings.writeInterval);
      requirePositive(BlockReceiverFields.samplingInterval, settings.samplingInterval);
      settings.fileName = requireFile(reader->readPath(BlockReceiverFields.file));
      outputs.push_back({.name = name, .settings = std::move(settings)});
      break;
    }
    case BlockKind::FaultReceivers: {
      auto settings = readFaultReceivers(reader, BlockFaultReceiverFields);
      requirePositive(BlockFaultReceiverFields.interval, settings.writeInterval);
      requireFile(settings.pickpointFileName);
      if (!withoutFault()) {
        outputs.push_back({.name = name, .settings = std::move(settings)});
      }
      break;
    }
    case BlockKind::Fault: {
      const auto settings = readFault(reader, BlockFaultFields, defaultTimeSeries);
      requirePositive(BlockFaultFields.interval, settings.printTimeIntervalSec);
      if (!withoutFault()) {
        outputs.push_back({.name = name, .settings = settings});
      }
      break;
    }
    }
  }
  return outputs;
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
  auto* reader = baseReader->readSubNode("elementwise");
  auto settings = readFault(reader, LegacyFaultFields, defaultTimeSeries);
  reader->warnDeprecated({"printintervalcriterion"});
  return settings;
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
  auto* reader = baseReader->readSubNode("output");
  auto enabled = reader->readWithDefault("surfaceoutput", false);
  auto settings = readSurface(reader, LegacySurfaceFields, defaultTimeSeries, configs);
  warnIntervalAndDisable(enabled, settings.interval, "surfaceoutput", LegacySurfaceFields.interval);
  if (!enabled) {
    return std::nullopt;
  }
  return settings;
}

PickpointParameters readPickpointParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("pickpoint");
  auto settings = readFaultReceivers(reader, LegacyFaultReceiverFields);
  if (reader->readWithDefault("aggregateperrank", false)) {
    logWarning() << "aggregateperrank is no longer supported: the text output writes a file per "
                    "on-fault receiver, and format = 'hdf5' writes all of them into a single file.";
  }
  reader->warnDeprecated({"noutpoints", "maxpickstore"});
  return settings;
}

std::optional<ReceiverOutputParameters> readReceiverParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("output");
  // a run writes the off-fault receivers if it gives their points, unless it switches them off
  auto enabled =
      reader->readWithDefault("receiveroutput", reader->hasField(LegacyReceiverFields.file));
  auto settings = readReceivers(reader, LegacyReceiverFields);
  warnIntervalAndDisable(
      enabled, settings.writeInterval, "receiveroutput", LegacyReceiverFields.interval);
  warnIntervalAndDisable(
      enabled, settings.samplingInterval, "receiveroutput", LegacyReceiverFields.samplingInterval);
  const auto fileName = reader->readPath(LegacyReceiverFields.file);
  if (enabled && !fileName.has_value()) {
    logError() << "The off-fault receiver output is enabled, but no receiver point file was given.";
  }
  if (!enabled) {
    return std::nullopt;
  }
  settings.fileName = fileName.value_or(settings.fileName);
  return settings;
}

std::optional<WaveFieldOutputParameters>
    readWaveFieldParameters(ParameterReader* baseReader,
                            std::optional<TimeSeriesMode> defaultTimeSeries,
                            const std::vector<ConfigId>& configs) {
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

  auto settings = readWaveField(reader, LegacyWaveFieldFields, defaultTimeSeries, configs);
  warnIntervalAndDisable(
      enabled, settings.interval, "wavefieldoutput", LegacyWaveFieldFields.interval);

  if (enabledPre.has_value()) {
    reader->warnDeprecated({"format"});
  }

  if (!enabled) {
    return std::nullopt;
  }
  return settings;
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

  // the outputs as the parameter file has always given them, each named as the files it has always
  // written
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

  // then those in blocks of their own
  for (auto& entry : readOutputBlocks(baseReader, defaultTimeSeries, configs, withFault)) {
    outputs.push_back(std::move(entry));
  }

  // every output writes files of its name
  std::unordered_set<std::string> names;
  for (const auto& entry : outputs) {
    if (!names.insert(entry.name).second) {
      logError() << "Two outputs are named" << entry.name
                 << "; give one of them another one (name = ...).";
    }
  }

  // the clusters of the fault sample its outputs into one place each, so far
  const auto atMostOne = [&](std::size_t count, const char* kind) {
    if (count > 1) {
      logError() << "A run can have only one" << kind << "output for now, but this one has" << count
                 << ".";
    }
  };
  atMostOne(countOf<ElementwiseFaultParameters>(outputs), "fault");
  atMostOne(countOf<PickpointParameters>(outputs), "on-fault receiver");

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
