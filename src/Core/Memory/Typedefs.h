// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_CORE_MEMORY_TYPEDEFS_H_
#define SEISSOL_SRC_CORE_MEMORY_TYPEDEFS_H_

#include "Alignment.h"
#include "CellLocalInformation.h"
#include "Config.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/Constants.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/Real.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"

#include <cstddef>

namespace seissol {

namespace kernels {
template <typename Cfg>
constexpr std::size_t NumSpaceQuadraturePoints =
    (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1);
} // namespace kernels

/**
 * The generated constant matrices of the configuration `Cfg`, as one table of pointers into the
 * pool of its generator.
 *
 * There are two of these per configuration: one built on the image in this binary, one on a
 * copy of it in device memory. Which entries a table holds is decided by the code generator.
 **/
template <typename Cfg>
using GlobalData = seissol::Pool<Cfg>;

/// The global data of the configuration `Cfg` on the host and, in a GPU build, on the device.
template <typename Cfg>
struct CompoundGlobalData {
  GlobalData<Cfg>* onHost{nullptr};
  GlobalData<Cfg>* onDevice{nullptr};
};

} // namespace seissol

#endif // SEISSOL_SRC_CORE_MEMORY_TYPEDEFS_H_
