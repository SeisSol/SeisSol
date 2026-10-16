// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Binary.h"

#include "Data.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {
YAML::Node BinaryWrite::serialize() {
  YAML::Node node;
  node["file"] = filename;
  node["source"] = dataSource->serialize();
  node["writer"] = "binary";
  node["alignment"] = alignment;
  node["append"] = append;
  if (offset.has_value()) {
    node["offset"] = offset.value();
  }
  return node;
}

BinaryWrite::BinaryWrite(const std::string& filename,
                         std::shared_ptr<writer::DataSource> dataSource,
                         std::size_t alignment,
                         bool append,
                         std::optional<std::size_t> offset)
    : filename(filename), dataSource(std::move(dataSource)), alignment(alignment), append(append),
      offset(offset) {}

BinaryWrite::BinaryWrite(YAML::Node node)
    : filename(node["file"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])),
      alignment(node["alignment"].as<std::size_t>()), append(node["append"].as<bool>()),
      offset(node["offset"] ? node["offset"].as<std::size_t>() : std::optional<std::size_t>()) {}

std::vector<std::shared_ptr<DataSource>> BinaryWrite::dataSources() { return {dataSource}; }

} // namespace seissol::io::writer::instructions
