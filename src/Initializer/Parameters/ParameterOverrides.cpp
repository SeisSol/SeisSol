// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ParameterOverrides.h"

#include "Initializer/Parameters/ParameterReader.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <exception>
#include <fty/fty.hpp>
#include <optional>
#include <string>
#include <string_view>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::initializer::parameters {

namespace {
//! @brief Whether @p name can be the name of a section or a field, which are words.
bool isName(std::string_view name) {
  return !name.empty() && std::ranges::all_of(name, [](char character) {
    return std::isalnum(static_cast<unsigned char>(character)) != 0 || character == '_';
  });
}

/**
 * @brief The value of @p entry as the YAML node reading the parameter file makes of it.
 *
 * It goes through the same converter, as the single field of a section of its own, so that
 * quotes, Fortran exponents and comments are taken care of exactly as in the file.
 */
YAML::Node convertValue(const ParameterOverride& entry) {
  const fty::StringsT lines{"&override", "value = " + entry.value, "/"};
  try {
    fty::Converter<fty::AsLowercase> converter;
    const auto converted = converter.convertToYAML(lines);
    if (converted["override"] && converted["override"]["value"]) {
      return converted["override"]["value"];
    }
  } catch (const std::exception& error) {
    logError() << "The value of the override" << entry.name()
               << "cannot be read:" << std::string(error.what());
  }
  logError() << "The value of the override" << entry.name() << "cannot be read.";
  return {};
}

void setAt(YAML::Node node,
           const std::vector<std::string>& path,
           std::size_t level,
           const YAML::Node& value) {
  if (level + 1 == path.size()) {
    node[path[level]] = value;
  } else {
    setAt(node[path[level]], path, level + 1, value);
  }
}
} // namespace

std::string ParameterOverride::name() const {
  std::string text;
  for (const auto& part : path) {
    if (!text.empty()) {
      text += '.';
    }
    text += part;
  }
  return text;
}

std::optional<ParameterOverride> parseOverride(std::string_view text) {
  const auto equals = text.find('=');
  if (equals == std::string_view::npos) {
    return {};
  }
  ParameterOverride entry;
  entry.value = std::string(text.substr(equals + 1));

  const auto name = text.substr(0, equals);
  std::size_t begin = 0;
  while (begin <= name.size()) {
    const auto dot = name.find('.', begin);
    const auto end = dot == std::string_view::npos ? name.size() : dot;
    const auto part = name.substr(begin, end - begin);
    if (!isName(part)) {
      return {};
    }
    std::string lower(part);
    for (auto& character : lower) {
      character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    }
    entry.path.push_back(lower);
    begin = end + 1;
  }
  // every field of a parameter file sits in a section
  if (entry.path.size() < 2) {
    return {};
  }
  return entry;
}

void applyOverrides(const YAML::Node& parameters, const std::vector<ParameterOverride>& overrides) {
  for (const auto& entry : overrides) {
    setAt(parameters, entry.path, 0, convertValue(entry));
  }
}

std::vector<ParameterOverride> unknownOverrides(const ParameterReader& reader,
                                                const std::vector<ParameterOverride>& overrides) {
  std::vector<ParameterOverride> unknown;
  for (const auto& entry : overrides) {
    if (!reader.visited(entry.path)) {
      unknown.push_back(entry);
    }
  }
  return unknown;
}

} // namespace seissol::initializer::parameters
