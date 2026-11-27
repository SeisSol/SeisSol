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
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <sstream>
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
  std::ranges::sort(ofGroups);
  std::vector<ConfigId> result{config};
  for (const auto groupConfig : ofGroups) {
    if (std::ranges::find(result, groupConfig) == result.end()) {
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

// the values numflux and numfluxnearfault take: the fluxes, and the presets of the general flux
// ('mixed'), which set its parts themselves
enum class FluxValue : std::uint8_t { Godunov, Rusanov, Centered, Mixed, GodunovShearCentered };

// godunov-shearcentered: the Godunov flux without the damping of the tangential traction jumps,
// and a tenth of the damping of the Rusanov flux, without which the orders 5 and 6 are unstable
constexpr std::array<double, 6> ShearCenteredDamping{1, 1, 1, 0, 1, 1};
constexpr double ShearCenteredRusanovShare = 0.1;

// four shares (vn vt tn tt), or six with those of a fluid (its normal velocity, its pressure),
// which take vn and tn if not given
std::array<double, 6> readDampingShares(ParameterReader* reader, const std::string& name) {
  const auto text = reader->readWithDefault(name, std::string("1 1 1 1"));
  std::istringstream stream(text);
  std::vector<double> values;
  double value = 0;
  while (stream >> value) {
    values.push_back(value);
  }
  if (!stream.eof() || (values.size() != 4 && values.size() != 6)) {
    logError() << name.c_str() << "takes four shares (vn vt tn tt) or six (vn vt tn tt vf p), not"
               << text;
  }
  std::array<double, 6> shares{};
  std::copy(values.begin(), values.end(), shares.begin());
  if (values.size() == 4) {
    shares[4] = shares[0];
    shares[5] = shares[2];
  }
  for (const auto share : shares) {
    if (!(share >= 0.0)) {
      logError() << "The shares of" << name.c_str() << "have to be at least 0.";
    }
  }
  return shares;
}

} // namespace

FluxSetting readFluxSetting(ParameterReader* reader, const std::string& suffix) {
  const auto name = "numflux" + suffix;
  const auto gpartName = "fluxgpart" + suffix;
  const auto rpartName = "fluxrpart" + suffix;
  const auto dampingName = "fluxdamping" + suffix;
  const auto value = reader->readWithDefaultStringEnum<FluxValue>(
      name,
      "godunov",
      {
          {"godunov", FluxValue::Godunov},
          {"upwind", FluxValue::Godunov},
          {"rusanov", FluxValue::Rusanov},
          {"centered", FluxValue::Centered},
          {"mixed", FluxValue::Mixed},
          {"godunov-shearcentered", FluxValue::GodunovShearCentered},
      });
  if (value == FluxValue::GodunovShearCentered) {
    for (const auto& part : {gpartName, rpartName, dampingName}) {
      if (reader->hasField(part)) {
        logError() << name.c_str() << "= 'godunov-shearcentered' sets" << part.c_str()
                   << "itself; to set the parts of the flux, take 'mixed'.";
      }
    }
    return FluxSetting{.flux = NumericalFlux::Mixed,
                       .gpart = 1.0,
                       .rpart = ShearCenteredRusanovShare,
                       .damping = ShearCenteredDamping};
  }
  const auto flux = value == FluxValue::Rusanov    ? NumericalFlux::Rusanov
                    : value == FluxValue::Centered ? NumericalFlux::Centered
                    : value == FluxValue::Mixed    ? NumericalFlux::Mixed
                                                   : NumericalFlux::Godunov;
  return FluxSetting{.flux = flux,
                     .gpart = reader->readWithDefault(gpartName, 0.0),
                     .rpart = reader->readWithDefault(rpartName, 0.0),
                     .damping = readDampingShares(reader, dampingName)};
}

bool FaceFluxRule::selects(std::uint32_t tag, int group, int neighborGroup) const {
  const auto inGroups = [&](int candidate) { return groups->at(candidate).has_value(); };
  return (!tags.has_value() || tags->at(tag).has_value()) &&
         (!groups.has_value() || inGroups(group) || inGroups(neighborGroup));
}

FluxSetting ModelParameters::fluxSetting(bool nearFault) const {
  return nearFault
             ? FluxSetting{.flux = fluxNearFault,
                           .gpart = fluxNearFaultGPart,
                           .rpart = fluxNearFaultRPart,
                           .damping = fluxNearFaultDamping}
             : FluxSetting{
                   .flux = flux, .gpart = fluxGPart, .rpart = fluxRPart, .damping = fluxDamping};
}

std::optional<std::size_t> ModelParameters::fluxRuleOf(std::uint32_t tag,
                                                       int group,
                                                       int neighborGroup,
                                                       bool nearFault) const {
  for (std::size_t i = 0; i < fluxRules.size(); ++i) {
    const auto& rule = fluxRules[i];
    const bool sets = nearFault ? rule.fluxNearFault.has_value() : rule.flux.has_value();
    if (sets && rule.selects(tag, group, neighborGroup)) {
      return i;
    }
  }
  return {};
}

FluxSetting
    ModelParameters::fluxOf(std::uint32_t tag, int group, int neighborGroup, bool nearFault) const {
  const auto rule = fluxRuleOf(tag, group, neighborGroup, nearFault);
  if (!rule.has_value()) {
    return fluxSetting(nearFault);
  }
  return nearFault ? fluxRules[*rule].fluxNearFault.value() : fluxRules[*rule].flux.value();
}

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

ModelParameters readModelParameters(ParameterReader* baseReader,
                                    ConfigId config,
                                    std::unordered_map<int, ConfigId> groupConfigs) {
  auto* reader = baseReader->readSubNode("equations");

  const auto simulationCount = readSimulationCount(baseReader, config);

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
  const auto fusedCount = configValue(config).numSimulations;
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
  if (isAnelastic && (freqRatio <= 0)) {
    logError() << "The freqratio parameter must be positive; but that is currently not the case.";
  }

  const ITMParameters itmParameters = readITMParameters(baseReader);

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  const auto flux = readFluxSetting(reader, "");
  const auto fluxNearFault = readFluxSetting(reader, "nearfault");

  return ModelParameters{.hasBoundaryFile = hasBoundaryFile,
                         .plasticity = plasticity,
                         .plasticityPointwise = plasticityPointwise,
                         .plasticityDisabledGroups = plasticityDisabledGroups,
                         .useCellHomogenizedMaterial = useCellHomogenizedMaterial,
                         .freqCentral = freqCentral,
                         .freqRatio = freqRatio,
                         .gravitationalAcceleration = gravitationalAcceleration,
                         .tv = tv,
                         .boundaryFileNames = boundaryFileNames,
                         .materialFileNames = materialFileNames,
                         .plasticityFileNames = plasticityFileNames,
                         .itmParameters = itmParameters,
                         .flux = flux.flux,
                         .fluxNearFault = fluxNearFault.flux,
                         .fluxRPart = flux.rpart,
                         .fluxNearFaultRPart = fluxNearFault.rpart,
                         .fluxGPart = flux.gpart,
                         .fluxNearFaultGPart = fluxNearFault.gpart,
                         .fluxDamping = flux.damping,
                         .fluxNearFaultDamping = fluxNearFault.damping,
                         .config = config,
                         .groupConfigs = std::move(groupConfigs),
                         .simulationCount = simulationCount};
}

std::string fluxToString(NumericalFlux flux) {
  if (flux == NumericalFlux::Godunov) {
    return "Godunov flux";
  }
  if (flux == NumericalFlux::Rusanov) {
    return "Rusanov flux";
  }
  if (flux == NumericalFlux::Centered) {
    return "Centered flux";
  }
  if (flux == NumericalFlux::Mixed) {
    return "Mixed flux (consult parameter file)";
  }
  return "(unknown flux)";
}

std::string describeFlux(NumericalFlux flux,
                         double gpart,
                         double rpart,
                         const std::array<double, 6>& damping) {
  std::ostringstream text;
  const bool upwind = std::ranges::all_of(damping, [](double share) { return share == 1.0; });
  const auto godunov = [&]() {
    text << "Godunov flux";
    if (!upwind) {
      text << " (shares of the damping";
      for (const auto share : damping) {
        text << " " << share;
      }
      text << ")";
    }
  };
  if (flux == NumericalFlux::Godunov) {
    godunov();
  } else if (flux == NumericalFlux::Rusanov) {
    text << "Rusanov flux";
  } else if (flux == NumericalFlux::Centered) {
    text << "centered flux";
  } else {
    if (gpart != 0.0) {
      text << gpart << " x ";
      godunov();
    }
    if (gpart != 1.0) {
      text << (gpart != 0.0 ? " + " : "") << 1.0 - gpart << " x centered flux";
    }
    if (rpart != 0.0) {
      text << " + " << rpart << " x damping of the Rusanov flux";
    }
  }
  return text.str();
}

} // namespace seissol::initializer::parameters
