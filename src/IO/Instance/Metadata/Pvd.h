// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_H_
#define SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_H_

#include "Xml.h"

namespace seissol::io::instance::metadata {

struct PvuEntry {
  std::string file;
  double timestep;
};

XmlFile makePvu(const std::vector<PvuEntry>& entries);

} // namespace seissol::io::instance::metadata

#endif // SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_H_
