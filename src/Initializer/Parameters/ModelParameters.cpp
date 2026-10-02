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
#include <cstddef>
#include <string>
#include <unordered_set>
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

ConfigId readConfig(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("equations");
  auto name = reader->read<std::string>("configuration");
  if (!name.has_value()) {
    return defaultConfig();
  }
  sanitize(name.value());
  const auto config = findConfig(name.value());
  if (!config.has_value()) {
    std::string names;
    for (std::size_t id = 0; id < builtConfigCount(); ++id) {
      names += (id > 0 ? ", " : "") + configName(configValue(static_cast<ConfigId>(id)));
    }
    logError() << "The configuration" << name.value()
               << "is not built into this executable. It has:" << names;
  }
  return config.value();
}

ModelParameters readModelParameters(ParameterReader* baseReader, ConfigId config) {
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

  const bool isAnelastic = configValue(config).relaxationMechanisms > 0;

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
                         {}};
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
