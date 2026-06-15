// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_FACEMAP_H_
#define SEISSOL_SRC_INITIALIZER_FACEMAP_H_

#include "Common/SegmentMap.h"
#include "Initializer/BasicTypedefs.h"

namespace YAML {
class Node;
} // namespace YAML

namespace seissol {

using FaceMap = SegmentMap<uint32_t, FaceType>;

FaceMap parseFaceMap(const YAML::Node& node);
FaceMap defaultFaceMap();

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_FACEMAP_H_
