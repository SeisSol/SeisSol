// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONFIGLAYOUT_H_
#define SEISSOL_SRC_COMMON_CONFIGLAYOUT_H_

#include "Common/ConfigValue.h"

#include <cstddef>

namespace seissol {

/**
 * @brief What a configuration implies for the data of a build.
 *
 * The sizes follow from the configuration together with the target the kernels were generated
 * for (e.g. the padding to the vector width), which is why they are not part of the
 * `ConfigValue`. They are taken from the generated code, so that code which only knows the
 * configuration at runtime can size its data without including that code.
 */
struct ConfigLayout {
  /// the configuration the sizes belong to
  ConfigValue config;

  /// quantities of the material
  std::size_t numQuantities{};
  /// basis functions of a cell, per quantity and simulation
  std::size_t numBasisFunctions{};
  /// position of the basis-function index in the tensors of a cell; fused simulations come first
  std::size_t basisFunctionDimension{};
  /// stored entries of the unknowns of one cell, including the padding; a solver that keeps the
  /// anelastic unknowns apart stores them elsewhere
  std::size_t dofsSize{};

  /// quadrature points on a dynamic rupture face
  std::size_t drNumPoints{};
  /// quadrature points on a dynamic rupture face as stored: padded, and for all fused simulations
  std::size_t drNumPaddedPoints{};
  /// quantities interpolated to a dynamic rupture face
  std::size_t drNumQuantities{};
  /// time points at which a dynamic rupture face is evaluated per time step
  std::size_t drNumTimePoints{};
};

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONFIGLAYOUT_H_
