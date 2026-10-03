// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_PLASTICITYQUANTITIES_H_
#define SEISSOL_SRC_MODEL_PLASTICITYQUANTITIES_H_

#include <array>
#include <cstddef>
#include <string>

namespace seissol::model {

/// The number of quantities of the plastic strain of a cell, the same for every configuration.
inline constexpr std::size_t PlasticityQuantityCount = 7;

/// The names of the quantities of the plastic strain of a cell.
inline const std::array<std::string, PlasticityQuantityCount> PlasticityQuantities = {
    "ep_xx", "ep_yy", "ep_zz", "ep_xy", "ep_yz", "ep_xz", "eta"};

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_PLASTICITYQUANTITIES_H_
