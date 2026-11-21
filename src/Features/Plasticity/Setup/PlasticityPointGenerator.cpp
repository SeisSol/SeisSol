// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "PlasticityPointGenerator.h"

#include "Core/Common/Constants.h"
#include "Core/Data/DataTable.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/QueryGenerator.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace seissol::initializer {

std::size_t PlasticityPointGenerator::outputPerCell() const {
  return pointwise_ ? nodes_.size() : 1;
}

data::DataTable PlasticityPointGenerator::generate() const {
  const auto pointsPerCell = outputPerCell();

  // the plasticity nodes of every element, or its barycenter only
  auto referencePoints =
      pointwise_ ? nodes_ : std::vector<std::array<double, Cell::Dim>>{{1 / 4., 1 / 4., 1 / 4.}};

  data::DataTable table(cellToVertex_.size * pointsPerCell);
  bindPointSet(
      table,
      [this, referencePoints = std::move(referencePoints)](
          std::size_t first, std::size_t count, double* xyz) {
        cellPointCoordinates(cellToVertex_, referencePoints, first, count, xyz);
      },
      [this, pointsPerCell](std::size_t first, std::size_t count, std::int32_t* out) {
        for (std::size_t i = 0; i < count; ++i) {
          out[i] = cellToVertex_.elementGroups((first + i) / pointsPerCell);
        }
      });
  return table;
}

} // namespace seissol::initializer
