// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Common/Real.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Model/BoundaryMappings.h"
#include "Numerical/Projection.h"
#include "Solver/MultipleSimulations.h"
#include "TestConfigs.h"
#include "TestHelper.h"

#include <algorithm>
#include <cstddef>
#include <limits>

namespace seissol::unit_test {

// The nodes of the boundary faces are the nodal points of the face basis, whichever layout the
// configuration gives the generated nodes2D.
TEST_CASE_TEMPLATE("The reference nodes of a boundary face are the nodal points of the face basis",
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  // what the code generator ships is stored in the precision of the configuration
  constexpr double Tolerance = std::max(1e-10, 10.0 * std::numeric_limits<Real<Cfg>>::epsilon());
  const auto points = numerical::projection::nodalPoints2D(Cfg::ConvergenceOrder);
  REQUIRE(points.size() == nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisDim<Cfg>]);
  for (std::size_t p = 0; p < points.size(); ++p) {
    const auto node = initializer::referenceFaceNode<Cfg>(p);
    for (std::size_t d = 0; d < node.size(); ++d) {
      REQUIRE(node[d] == AbsApprox(points[p][d]).epsilon(Tolerance).delta(Tolerance));
    }
  }
}

} // namespace seissol::unit_test
