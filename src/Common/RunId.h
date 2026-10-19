// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_RUNID_H_
#define SEISSOL_SRC_COMMON_RUNID_H_

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace seissol {

/**
 * @brief Identifies one run of SeisSol: a UUID in its canonical text form, lower case.
 *
 * Every file a run writes carries it, so that an output found later can be traced back to the
 * run, its log and its summary.
 */
constexpr std::size_t RunIdLength = 36;

//! @brief A new random (version 4) UUID.
std::string newRunId();

/**
 * @brief @p text as a run ID, or nothing if it is not a UUID.
 *
 * Accepts the canonical form, 8-4-4-4-12 hexadecimal digits, in either case, and returns it in
 * lower case.
 */
std::optional<std::string> parseRunId(std::string_view text);

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_RUNID_H_
