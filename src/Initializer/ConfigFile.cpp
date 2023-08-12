#include "ConfigFile.hpp"
#include <Common/cellconfig.hpp>
#include <Common/configs.hpp>
#include <type_traits>
#include <yaml-cpp/node/parse.h>
#include <yaml-cpp/yaml.h>
#include <string>

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
} // namespace

namespace seissol::initializer {

std::unordered_map<int, CellConfigInfo> readConfigFile(const std::string& filename) {
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
