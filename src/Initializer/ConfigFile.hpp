#pragma once

#include "Common/configs.hpp"
#include <unordered_map>
namespace seissol::initializer {

struct CellConfigInfo {
  std::size_t id;
  SupportedConfigs config;
  std::string model;
};

using CellConfigInfoMap = std::unordered_map<int, CellConfigInfo>;

std::unordered_map<int, CellConfigInfo> readConfigFile(const std::string& filename);

} // namespace seissol::initializer
