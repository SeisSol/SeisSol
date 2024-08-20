#pragma once

#include "Xml.hpp"

namespace seissol::io::instance::metadata {

struct PvuEntry {
  std::string file;
  double timestep;
};

inline XmlFile makePvu(const std::vector<PvuEntry>& entries) {
  XmlFile file;

  auto root = XmlNode("VTKFile");
  root.addAttribute(XmlAttribute::create("type", "collection"))
      .addAttribute(XmlAttribute::create("version", "0.1"));

  auto collection = XmlNode("Collection");

  for (const auto& entry : entries) {
    auto dataset = XmlData("DataSet");
    dataset.addAttribute(XmlAttribute::create("timestep", entry.timestep))
        .addAttribute(XmlAttribute::create("group", ""))
        .addAttribute(XmlAttribute::create("part", 0))
        .addAttribute(XmlAttribute::create("file", entry.file));
    collection.addNode(std::move(dataset));
  }

  root.addNode(std::move(collection));
  file.setRoot(std::move(root));

  return file;
}

} // namespace seissol::io::instance::metadata
