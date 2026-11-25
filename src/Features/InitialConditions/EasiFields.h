// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_EASIFIELDS_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_EASIFIELDS_H_

#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Easi/EasiContext.h"

#include <string>
#include <vector>

namespace seissol::initializer {

/// The values of the easi fields `iniFields` at the quadrature points of every element, as the
/// configuration `Cfg` projects them: by element, point, quantity and field.
template <typename Cfg>
std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime,
                                      const data::EasiContext& easi);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_EASIFIELDS_H_
