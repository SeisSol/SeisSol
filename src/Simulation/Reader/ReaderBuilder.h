// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_READER_READERBUILDER_H_
#define SEISSOL_SRC_SIMULATION_READER_READERBUILDER_H_

#include "Core/Data/DataReader.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::data {

/// The reader of the model `path`, which reads the columns `defaultInArgs` (e.g. x, y and z): an
/// easi file, with or without the prefix `easi:`.
std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs);

/// The loader of the models of a run (material, fault and boundary condition files): buildReader.
ModelLoader modelLoader();

} // namespace seissol::data

#endif // SEISSOL_SRC_SIMULATION_READER_READERBUILDER_H_
