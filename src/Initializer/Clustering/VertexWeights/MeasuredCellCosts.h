// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_MEASUREDCELLCOSTS_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_MEASUREDCELLCOSTS_H_

#include <cstddef>
#include <vector>

namespace seissol::initializer {

/**
 * @brief What the cells of this rank cost in a run, as measured there: the time the loops of the
 * time stepping spent on each, in seconds per simulated second (see measureCellCosts()).
 */
struct MeasuredCellCosts {
  //! the cells, by their ids in the mesh file; in the order of the elements of the mesh reader
  std::vector<std::size_t> ids;
  std::vector<double> costs;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_MEASUREDCELLCOSTS_H_
