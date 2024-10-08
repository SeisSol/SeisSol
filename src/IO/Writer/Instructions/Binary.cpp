// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Binary.h"

#include "Data.h"
#include <memory>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {
YAML::Node BinaryWrite::serialize() {
  YAML::Node node;
  node["file"] = filename;
  node["source"] = dataSource->serialize();
  node["writer"] = "binary";
  return node;
}

BinaryWrite::BinaryWrite(const std::string& filename,
                         std::shared_ptr<writer::DataSource> dataSource)
    : filename(filename), dataSource(dataSource) {}

BinaryWrite::BinaryWrite(YAML::Node node)
    : filename(node["file"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])) {}

std::vector<std::shared_ptr<DataSource>> BinaryWrite::dataSources() { return {dataSource}; }

} // namespace seissol::io::writer::instructions