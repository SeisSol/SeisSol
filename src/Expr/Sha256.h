// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_SHA256_H_
#define SEISSOL_SRC_EXPR_SHA256_H_

// SHA-256 (FIPS 180-4), for the provenance of the core format: the hashes of the sources a core
// file was generated from, and of its body (`meta source`, `depends`, `body`).

#include <cstddef>
#include <cstdint>
#include <string>

namespace seissol::expr {

/// The SHA-256 of `data`, as 64 lowercase hexadecimal digits.
[[nodiscard]] std::string sha256(const std::string& data);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_SHA256_H_
