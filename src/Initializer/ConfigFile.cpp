#include "ConfigFile.hpp"
#include <Common/cellconfig.hpp>
#include <Common/cellconfigconv.hpp>
#include <Common/configs.hpp>
#include <type_traits>
#include <yaml-cpp/node/parse.h>
#include <yaml-cpp/yaml.h>
#include <string>
#include "utils/logger.h"
#include "Initializer/InputParameters.hpp"
#include "SeisSol.h"

namespace {
template <std::size_t I>
static std::pair<std::size_t, seissol::SupportedConfigs> testConfig(const std::string& material,
                                                                    const std::string& precision,
                                                                    int order,
                                                                    bool plasticity) {
  if constexpr (I < std::variant_size_v<seissol::SupportedConfigs>) {
    using ConfigI = std::variant_alternative_t<I, seissol::SupportedConfigs>;
    bool precisionMatches = precision == seissol::PrecisionFromType<typename ConfigI::RealT>::Text;
    bool materialMatches = material == ConfigI::MaterialT::Text;
    bool orderMatches = ConfigI::ConvergenceOrder == order;
    bool plasticityMatches = ConfigI::Plasticity == plasticity;
    if (precisionMatches && materialMatches && orderMatches && plasticityMatches) {
      return std::pair(I, seissol::SupportedConfigs(ConfigI()));
    } else {
      return testConfig<I + 1>(material, precision, order, plasticity);
    }
  } else {
    // TODO(David): make more descriptive
    logError() << "Unknown cell configuration.";
    throw std::runtime_error("Unknown cell configuration.");
  }
}

constexpr SupportedConfigs defaultConfig(bool plasticity) {
  if (plasticity) {
    return SupportedConfigs(CellConfig<seissol::model::Material_t, real, ConvergenceOrder, true>());
  } else {
    return SupportedConfigs(
        CellConfig<seissol::model::Material_t, real, ConvergenceOrder, false>());
  }
}

static std::unordered_map<int, seissol::initializer::CellConfigInfo> fallbackConfig() {
  const auto& parameters = seissol::SeisSol::main.getSeisSolParameters();
  std::unordered_map<int, seissol::initializer::CellConfigInfo> configs;
  configs[0] = initializer::CellConfigInfo{0, // TODO(David): fix
                                           defaultConfig(parameters.model.plasticity),
                                           parameters.model.materialFileName};
  return configs;
}
} // namespace

namespace seissol::initializer {

std::unordered_map<int, CellConfigInfo> readConfigFile(const std::string& filename) {
  // TODO(David): remove, once we go to a new input format
  if (filename == "") {
    return fallbackConfig();
  }

  YAML::Node root = YAML::LoadFile(filename);
  std::unordered_map<int, CellConfigInfo> configs;
  for (const auto& entries : root) {
    auto group = entries["group"].as<int>();
    auto material = entries["material"].as<std::string>();
    auto precision = entries["precision"].as<std::string>();
    auto order = entries["order"].as<int>();
    auto plasticity = entries["plasticity"].as<bool>();
    auto model = entries["model"].as<std::string>();

    auto [id, config] = testConfig<0>(material, precision, order, plasticity);
    configs[group] = CellConfigInfo{id, config, model};
  }
  return configs;
}
} // namespace seissol::initializer
