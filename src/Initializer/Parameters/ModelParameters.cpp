// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ModelParameters.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Initializer/Parameters/ParameterReader.h"

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

ConfigId ModelParameters::configOfGroup(int group) const {
  const auto found = groupConfigs.find(group);
  return found == groupConfigs.end() ? config : found->second;
}

std::vector<ConfigId> ModelParameters::configs() const {
  std::vector<ConfigId> ofGroups;
  ofGroups.reserve(groupConfigs.size());
  for (const auto& [group, groupConfig] : groupConfigs) {
    ofGroups.push_back(groupConfig);
  }
  std::sort(ofGroups.begin(), ofGroups.end());
  std::vector<ConfigId> result{config};
  for (const auto groupConfig : ofGroups) {
    if (std::find(result.begin(), result.end(), groupConfig) == result.end()) {
      result.push_back(groupConfig);
    }
  }
  return result;
}

ITMParameters readITMParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("equations");
  const auto itmEnabled = reader->readWithDefault<bool>("itmenable", false);
  const auto itmStartingTime = reader->readWithDefault<double>("itmstartingtime", 0.0);
  const auto itmDuration = reader->readWithDefault<double>("itmtime", 0.0);
  const auto itmVelocityScalingFactor =
      reader->readWithDefault<double>("itmvelocityscalingfactor", 1.0);
  const auto reflectionType =
      reader->readWithDefaultEnum<ReflectionType>("itmreflectiontype",
                                                  ReflectionType::BothWaves,
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
  return ITMParameters{
      itmEnabled, itmStartingTime, itmDuration, itmVelocityScalingFactor, reflectionType};
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

ModelParameters readModelParameters(ParameterReader* baseReader,
                                    ConfigId config,
                                    std::unordered_map<int, ConfigId> groupConfigs) {
  auto* reader = baseReader->readSubNode("equations");

  const auto boundaryFileName = reader->readPath("boundaryfilename");
  const std::string materialFileName =
      reader->readPathOrFail("materialfilename", "No material file given.");
  std::vector<std::string> plasticityFileNames(configValue(config).numSimulations);

  for (std::size_t i = 0; i < plasticityFileNames.size(); ++i) {
    const auto fieldname = "plasticityfilename" + (i == 0 ? std::string{} : std::to_string(i));
    plasticityFileNames[i] = reader->readPath(fieldname).value_or(materialFileName);
  }

  const bool hasBoundaryFile = !boundaryFileName.value_or("").empty();

  const bool plasticity = reader->readWithDefault("plasticity", false);

  const bool plasticityPointwise = reader->readWithDefault("plasticitypointwise", true);

  const auto plasticityDisabledGroupsRaw =
      reader->readWithDefault<std::string>("plasticitydisabledgroups", "");
  std::unordered_set<int> plasticityDisabledGroups;
  {
    const auto groups = utils::StringUtils::split(plasticityDisabledGroupsRaw, ',');
    for (const auto& group : groups) {
      plasticityDisabledGroups.emplace(std::stoi(group));
    }
  }

  const bool useCellHomogenizedMaterial =
      reader->readWithDefault("usecellhomogenizedmaterial", true);

  const double gravitationalAcceleration =
      reader->readWithDefault("gravitationalacceleration", 9.81);
  const double tv = reader->readWithDefault("tv", 0.1);

  bool isAnelastic = configValue(config).relaxationMechanisms > 0;
  for (const auto& [group, groupConfig] : groupConfigs) {
    isAnelastic = isAnelastic || configValue(groupConfig).relaxationMechanisms > 0;
  }

  const auto freqCentral = reader->readIfRequired<double>("freqcentral", isAnelastic);
  const auto freqRatio = reader->readIfRequired<double>("freqratio", isAnelastic);
  if (isAnelastic) {
    if (freqRatio <= 0) {
      logError() << "The freqratio parameter must be positive; but that is currently not the case.";
    }
  }

  const ITMParameters itmParameters = readITMParameters(baseReader);

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  const auto flux =
      reader->readWithDefaultStringEnum<NumericalFlux>("numflux",
                                                       "godunov",
                                                       {
                                                           {"godunov", NumericalFlux::Godunov},
                                                           {"rusanov", NumericalFlux::Rusanov},
                                                       });

  const auto fluxNearFault =
      reader->readWithDefaultStringEnum<NumericalFlux>("numfluxnearfault",
                                                       "godunov",
                                                       {
                                                           {"godunov", NumericalFlux::Godunov},
                                                           {"rusanov", NumericalFlux::Rusanov},
                                                       });

  return ModelParameters{hasBoundaryFile,
                         plasticity,
                         plasticityPointwise,
                         plasticityDisabledGroups,
                         useCellHomogenizedMaterial,
                         freqCentral,
                         freqRatio,
                         gravitationalAcceleration,
                         tv,
                         boundaryFileName.value_or(""),
                         materialFileName,
                         plasticityFileNames,
                         itmParameters,
                         flux,
                         fluxNearFault,
                         config,
                         std::move(groupConfigs)};
}

std::string fluxToString(NumericalFlux flux) {
  if (flux == NumericalFlux::Godunov) {
    return "Godunov flux";
  }
  if (flux == NumericalFlux::Rusanov) {
    return "Rusanov flux";
  }
  return "(unknown flux)";
}

} // namespace seissol::initializer::parameters
