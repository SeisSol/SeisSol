// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_HPP_
#define SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_HPP_

#include "Xml.hpp"

namespace seissol::io::instance::metadata {

struct PvuEntry {
  std::string file;
  double timestep;
};

XmlFile makePvu(const std::vector<PvuEntry>& entries);

} // namespace seissol::io::instance::metadata

#endif // SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_HPP_
