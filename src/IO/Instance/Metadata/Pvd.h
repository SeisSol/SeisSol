// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_H_
#define SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_H_

#include "Xml.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::io::instance::metadata {

struct PvuEntry {
  std::string file;
  double timestep;
};

//! The data set of a step of the collection.
std::shared_ptr<XmlEntry> makePvuDataset(const PvuEntry& entry);

//! The collection holding @p datasets.
XmlFile makePvuCollection(const std::vector<std::shared_ptr<XmlEntry>>& datasets);

XmlFile makePvu(const std::vector<PvuEntry>& entries);

} // namespace seissol::io::instance::metadata

#endif // SEISSOL_SRC_IO_INSTANCE_METADATA_PVD_H_
