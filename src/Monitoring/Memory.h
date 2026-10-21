// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_MEMORY_H_
#define SEISSOL_SRC_MONITORING_MEMORY_H_

#include <cstddef>
#include <optional>

namespace seissol::monitoring {

/**
 * @brief The most host memory this process has held so far (its peak resident set), in bytes.
 *
 * Nothing if the operating system does not say.
 */
std::optional<std::size_t> peakHostMemory();

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_MEMORY_H_
