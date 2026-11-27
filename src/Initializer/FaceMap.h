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
#include "Initializer/Parameters/ModelParameters.h"

#include <optional>
#include <string>
#include <vector>

namespace YAML {
class Node;
} // namespace YAML

namespace seissol {

using FaceMap = SegmentMap<uint32_t, FaceType>;

FaceMap parseFaceMap(const YAML::Node& node);

//! @brief The face map in the YAML file @p path , as parseFaceMap() reads it.
FaceMap readFaceMap(const std::string& path);

FaceMap defaultFaceMap();

/**
 * @brief The rules for the flux of regular faces in the face map @p node : its entries of regular
 * faces that are maps, select faces by their "tags", their mesh "groups", or both, and give the
 * flux keys of the parameter file (numflux, numfluxnearfault, fluxgpart, ...), in their order.
 */
std::vector<initializer::parameters::FaceFluxRule> parseFaceFluxRules(const YAML::Node& node);

//! @brief The rules of the face map in the YAML file @p path , none without a file.
std::vector<initializer::parameters::FaceFluxRule>
    readFaceFluxRules(const std::optional<std::string>& path);

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_FACEMAP_H_
