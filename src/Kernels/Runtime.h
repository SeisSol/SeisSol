// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_RUNTIME_H_
#define SEISSOL_SRC_KERNELS_RUNTIME_H_

#include "Config.h"
#include "GeneratedCode/runtime.h" // IWYU pragma: export
#include "GeneratedCode/variant.h"

#include <cstddef>

namespace seissol::kernels {

/// The variant of the kernels in runtime.h that computes for the configuration of this build.
///
/// Setup and output call their kernels through runtime.h, with views that carry the layout of
/// their values, in the variant they name; a build has one configuration and one variant.
constexpr std::size_t RuntimeVariant = runtime::variantOf<Config>();

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_RUNTIME_H_
