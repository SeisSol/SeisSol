// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ModelParameters.h"

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/ConfigValue.h"
#include "Core/Common/RunConfigs.h"
#include "Features/Plasticity/PlasticityParameters.h"
#include "Features/TimeMirror/ITMParameters.h"
#include "Simulation/Parameters/ParameterReader.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::initializer::parameters {

ITMParameters readITMParameters(ParameterReader* baseReader) {
  const ITMParameters defaults{};
  auto* reader = baseReader->readSubNode("equations");
  const auto itmEnabled = reader->readWithDefault("itmenable", defaults.itmEnabled);
  const auto itmStartingTime = reader->readWithDefault("itmstartingtime", defaults.itmStartingTime);
  const auto itmDuration = reader->readWithDefault("itmtime", defaults.itmDuration);
  const auto itmVelocityScalingFactor =
      reader->readWithDefault("itmvelocityscalingfactor", defaults.itmVelocityScalingFactor);
  const auto reflectionType =
      reader->readWithDefaultEnum<ReflectionType>("itmreflectiontype",
                                                  defaults.itmReflectionType,
                                                  {ReflectionType::BothWaves,
                                                   ReflectionType::BothWavesVelocity,
                                                   ReflectionType::Pwave,
                                                   ReflectionType::Swave});
  if (itmEnabled) {
    if (itmDuration <= 0.0) {
      logError() << "ITM Time is not positive. It should be positive!";
    }
    if (itmVelocityScalingFactor < 0.0) {
      logError() << "ITM Velocity Scaling Factor is less than zero. It should be positive!";
    }
    if (itmStartingTime < 0.0) {
      logError() << "ITM Starting Time can not be less than zero";
    }
  } else {
    reader->markUnused(
        {"itmstartingtime", "itmtime", "itmvelocityscalingfactor", "itmreflectiontype"});
  }
  return ITMParameters{.itmEnabled = itmEnabled,
                       .itmStartingTime = itmStartingTime,
                       .itmDuration = itmDuration,
                       .itmVelocityScalingFactor = itmVelocityScalingFactor,
                       .itmReflectionType = reflectionType};
}

namespace {

ConfigId configByName(std::string name) {
  sanitize(name);
  const auto config = findConfig(name);
  if (!config.has_value()) {
    logError() << "The configuration" << name
               << "is not built into this executable. It has:" << builtConfigNames();
  }
  return config.value();
}

} // namespace

ConfigId readConfig(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("equations");
  const auto name = reader->read<std::string>("configuration");
  if (!name.has_value()) {
    return defaultConfig();
  }
  const auto config = configByName(name.value());
  // the parameter file is the more specific of the two
  const auto fromEnvironment = environmentConfig();
  if (fromEnvironment.has_value() && fromEnvironment.value() != config) {
    logWarning() << "SEISSOL_CONFIGURATION names"
                 << configName(configValue(fromEnvironment.value()))
                 << "but the parameter file names" << configName(configValue(config))
                 << "which takes precedence.";
  }
  return config;
}

std::unordered_map<int, ConfigId> readGroupConfigs(ParameterReader* baseReader, ConfigId config) {
  auto* reader = baseReader->readSubNode("equations");
  const auto configMap = reader->readWithDefault<std::string>("configmap", "");

  std::unordered_map<int, ConfigId> groupConfigs;
  for (const auto& entry : utils::StringUtils::split(configMap, ';')) {
    auto trimmed = entry;
    if (utils::StringUtils::trim(trimmed).empty()) {
      continue;
    }
    const auto groupsAndName = utils::StringUtils::split(entry, ':');
    if (groupsAndName.size() != 2) {
      logError() << "The configmap entry" << entry
                 << "does not have the form \"group,group,...:configuration\".";
    }
    const auto groupConfig = configByName(groupsAndName[1]);
    if (configValue(groupConfig).numSimulations != configValue(config).numSimulations) {
      logError() << "The configuration" << configName(configValue(groupConfig))
                 << "fuses another number of simulations than" << configName(configValue(config))
                 << "; all configurations of a run fuse the same number.";
    }
    for (auto group : utils::StringUtils::split(groupsAndName[0], ',')) {
      utils::StringUtils::trim(group);
      int groupId = 0;
      const auto* groupEnd = group.data() + group.size();
      const auto parsed = std::from_chars(group.data(), groupEnd, groupId);
      if (group.empty() || parsed.ec != std::errc{} || parsed.ptr != groupEnd) {
        logError() << "The configmap entry" << entry << "names" << group
                   << "as a mesh group, which is not an integer.";
      }
      const auto [found, inserted] = groupConfigs.emplace(groupId, groupConfig);
      if (!inserted && found->second != groupConfig) {
        logError() << "The configmap gives the mesh group" << groupId
                   << "more than one configuration.";
      }
    }
  }
  return groupConfigs;
}

std::size_t readSimulationCount(ParameterReader* baseReader, ConfigId config) {
  auto* reader = baseReader->readSubNode("equations");
  const auto simulationCount =
      reader->readWithDefault<std::size_t>("simcount", configValue(config).numSimulations);
  if (simulationCount == 0) {
    logError() << "The number of simulations (simcount) has to be positive.";
  }
  return simulationCount;
}

RunConfigs readRunConfigs(ParameterReader* baseReader) {
  const ConfigId config = readConfig(baseReader);
  auto groupConfigs = readGroupConfigs(baseReader, config);
  const auto simulationCount = readSimulationCount(baseReader, config);
  return RunConfigs{.config = config,
                    .groupConfigs = std::move(groupConfigs),
                    .simulationCount = simulationCount};
}

ModelParameters readModelParameters(ParameterReader* baseReader, const RunConfigs& runConfigs) {
  auto* reader = baseReader->readSubNode("equations");

  const auto simulationCount = runConfigs.simulationCount;

  // the easi boundary condition of every simulation, the one of the first if not given
  std::vector<std::string> boundaryFileNames(simulationCount);
  boundaryFileNames[0] = reader->readPath("boundaryfilename").value_or("");
  for (std::size_t i = 1; i < simulationCount; ++i) {
    boundaryFileNames[i] =
        reader->readPath("boundaryfilename" + std::to_string(i)).value_or(boundaryFileNames[0]);
  }

  // The material of every simulation, the one of the first if not given. The plasticity of every
  // simulation likewise, the one of the first if not given, and else its own material file. An
  // empty file name counts as none.
  std::vector<std::string> materialFileNames(simulationCount);
  std::vector<std::string> plasticityFileNames(simulationCount);
  const auto plasticityFileName = reader->readPath("plasticityfilename").value_or("");
  for (std::size_t i = 0; i < simulationCount; ++i) {
    if (i == 0) {
      materialFileNames[i] = reader->readPathOrFail("materialfilename", "No material file given.");
    } else {
      materialFileNames[i] =
          reader->readPath("materialfilename" + std::to_string(i)).value_or(materialFileNames[0]);
    }
    plasticityFileNames[i] =
        i == 0 ? plasticityFileName
               : reader->readPath("plasticityfilename" + std::to_string(i)).value_or("");
    if (plasticityFileNames[i].empty()) {
      plasticityFileNames[i] = plasticityFileName;
    }
    if (plasticityFileNames[i].empty()) {
      plasticityFileNames[i] = materialFileNames[i];
    }
  }

  // The simulations fused into a cell share its material, and so have to name the same file.
  const auto fusedCount = configValue(runConfigs.config).numSimulations;
  for (std::size_t i = 0; i < simulationCount; ++i) {
    const auto first = i - i % fusedCount;
    if (materialFileNames[i] != materialFileNames[first]) {
      logError() << "The simulations" << first << "and" << i
                 << "are fused into the same cell, which has one material, but they name "
                    "different material files:"
                 << materialFileNames[first] << "for the first and" << materialFileNames[i]
                 << "for the second.";
    }
  }

  const bool hasBoundaryFile =
      std::ranges::any_of(boundaryFileNames, [](const auto& name) { return !name.empty(); });

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  return ModelParameters{.hasBoundaryFile = hasBoundaryFile,
                         .boundaryFileNames = boundaryFileNames,
                         .materialFileNames = materialFileNames,
                         .plasticityFileNames = plasticityFileNames};
}

WavePropagationParameters readWavePropagationParameters(ParameterReader* baseReader,
                                                        const RunConfigs& runConfigs) {
  const WavePropagationParameters defaults{};
  auto* reader = baseReader->readSubNode("equations");

  const bool useCellHomogenizedMaterial =
      reader->readWithDefault("usecellhomogenizedmaterial", defaults.useCellHomogenizedMaterial);

  const double gravitationalAcceleration =
      reader->readWithDefault("gravitationalacceleration", defaults.gravitationalAcceleration);

  const bool isAnelastic = std::ranges::any_of(runConfigs.configs(), [](ConfigId config) {
    return configValue(config).relaxationMechanisms > 0;
  });

  const auto freqCentral = reader->readIfRequired("freqcentral", isAnelastic, defaults.freqCentral);
  const auto freqRatio = reader->readIfRequired("freqratio", isAnelastic, defaults.freqRatio);
  if (isAnelastic && (freqRatio <= 0)) {
    logError() << "The freqratio parameter must be positive; but that is currently not the case.";
  }

  const auto flux =
      reader->readWithDefaultStringEnum<NumericalFlux>("numflux",
                                                       defaults.flux,
                                                       {
                                                           {"godunov", NumericalFlux::Godunov},
                                                           {"rusanov", NumericalFlux::Rusanov},
                                                       });

  const auto fluxNearFault =
      reader->readWithDefaultStringEnum<NumericalFlux>("numfluxnearfault",
                                                       defaults.fluxNearFault,
                                                       {
                                                           {"godunov", NumericalFlux::Godunov},
                                                           {"rusanov", NumericalFlux::Rusanov},
                                                       });

  return WavePropagationParameters{.flux = flux,
                                   .fluxNearFault = fluxNearFault,
                                   .freqCentral = freqCentral,
                                   .freqRatio = freqRatio,
                                   .useCellHomogenizedMaterial = useCellHomogenizedMaterial,
                                   .gravitationalAcceleration = gravitationalAcceleration};
}

PlasticityParameters readPlasticityParameters(ParameterReader* baseReader) {
  const PlasticityParameters defaults{};
  auto* reader = baseReader->readSubNode("equations");

  const bool plasticity = reader->readWithDefault("plasticity", defaults.enabled);

  const bool plasticityPointwise =
      reader->readWithDefault("plasticitypointwise", defaults.pointwise);

  const auto plasticityDisabledGroups = reader->readWithDefaultConverted(
      "plasticitydisabledgroups", defaults.disabledGroups, [](const std::string& text) {
        std::unordered_set<int> groups;
        for (const auto& group : utils::StringUtils::split(text, ',')) {
          groups.emplace(std::stoi(group));
        }
        return groups;
      });

  const double tv = reader->readWithDefault("tv", defaults.tv);

  return PlasticityParameters{.enabled = plasticity,
                              .pointwise = plasticityPointwise,
                              .disabledGroups = plasticityDisabledGroups,
                              .tv = tv};
}

} // namespace seissol::initializer::parameters
