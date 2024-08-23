#pragma once

#include "Xml.hpp"

namespace seissol::io::instance::metadata {

struct PvuEntry {
  std::string file;
  double timestep;
};

XmlFile makePvu(const std::vector<PvuEntry>& entries);

} // namespace seissol::io::instance::metadata
