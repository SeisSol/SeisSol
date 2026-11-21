// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_SETUP_PLASTICITYPOINTGENERATOR_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_SETUP_PLASTICITYPOINTGENERATOR_H_

#include "Core/Common/Constants.h"
#include "Core/Data/DataTable.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/QueryGenerator.h"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace seissol::initializer {

/// Queries the nodes of the plasticity of a configuration in each cell, given in reference
/// coordinates, or only the barycenter if not pointwise.
class PlasticityPointGenerator : public QueryGenerator {
  public:
  PlasticityPointGenerator(CellToVertexArray cellToVertex,
                           std::vector<std::array<double, Cell::Dim>> nodes,
                           bool pointwise = true)
      : cellToVertex_(std::move(cellToVertex)), nodes_(std::move(nodes)), pointwise_(pointwise) {}
  [[nodiscard]] data::DataTable generate() const override;
  [[nodiscard]] std::size_t outputPerCell() const override;

  private:
  CellToVertexArray cellToVertex_;
  std::vector<std::array<double, Cell::Dim>> nodes_;
  bool pointwise_{true};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_SETUP_PLASTICITYPOINTGENERATOR_H_
