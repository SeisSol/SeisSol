// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Pvd.hpp"

#include "Xml.hpp"
#include <memory>
#include <vector>

namespace seissol::io::instance::metadata {

XmlFile makePvu(const std::vector<PvuEntry>& entries) {
  XmlFile file;

  auto root = XmlNode("VTKFile");
  root.addAttribute(XmlAttribute::create("type", "collection"))
      .addAttribute(XmlAttribute::create("version", "0.1"));

  auto collection = std::make_shared<XmlNode>("Collection");

  for (const auto& entry : entries) {
    auto dataset = std::make_shared<XmlData>("DataSet");
    dataset->addAttribute(XmlAttribute::create("timestep", entry.timestep))
        .addAttribute(XmlAttribute::create("group", ""))
        .addAttribute(XmlAttribute::create("part", 0))
        .addAttribute(XmlAttribute::create("file", entry.file));
    collection->addNode(dataset);
  }

  root.addNode(collection);
  file.setRoot(std::make_shared<XmlNode>(root));

  return file;
}

} // namespace seissol::io::instance::metadata
