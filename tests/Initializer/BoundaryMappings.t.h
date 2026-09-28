// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Model/BoundaryMappings.h"
#include "Kernels/Precision.h"
#include "Numerical/Projection.h"
#include "Numerical/Transformation.h"
#include "Solver/MultipleSimulations.h"
#include "TestHelper.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace seissol::unit_test {

// The C++ code reads nodes2D as [node][chi/tau] in every configuration, fused simulations included.
static_assert(nodal::tensor::nodes2D::Shape[1] == 2, "nodes2D is stored as [node][chi/tau]");
static_assert(nodal::tensor::nodes2D::Shape[0] ==
                  tensor::INodal::Shape[multisim::BasisFunctionDimension],
              "nodes2D has one row per face node");

TEST_CASE("BoundaryMappings: boundary nodes are the face nodes mapped onto the face" *
          doctest::test_suite("initializer")) {
  constexpr std::size_t NodeCount = nodal::tensor::nodes2D::Shape[0];
  constexpr double Tolerance = std::max(1e-9, 100.0 * std::numeric_limits<real>::epsilon());

  // The run-time warp&blend points, i.e. independent of the generated nodes2D table. The first one
  // is the corner (0,0), which the generated tensor of fused builds leaves out of its bounding box.
  const auto reference = numerical::projection::nodalPoints2D(ConvergenceOrder);
  REQUIRE(reference.size() == NodeCount);

  const std::array<std::array<double, Cell::Dim>, Cell::NumVertices> vertices{
      std::array<double, Cell::Dim>{0.3, -0.2, 0.1},
      std::array<double, Cell::Dim>{2.1, 0.4, -0.3},
      std::array<double, Cell::Dim>{-0.5, 1.7, 0.2},
      std::array<double, Cell::Dim>{0.1, 0.3, 1.9}};
  const std::array<const double*, Cell::NumVertices> vertexPointers{
      vertices[0].data(), vertices[1].data(), vertices[2].data(), vertices[3].data()};

  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    std::array<real, NodeCount * Cell::Dim> nodes{};
    initializer::computeBoundaryNodes(vertexPointers, side, nodes.data());

    for (std::size_t p = 0; p < NodeCount; ++p) {
      std::array<double, Cell::Dim> xiEtaZeta{};
      transformations::chiTau2XiEtaZeta(
          static_cast<std::uint32_t>(side), reference[p].data(), xiEtaZeta.data());
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        const double expected = vertices[0][d] + xiEtaZeta[0] * (vertices[1][d] - vertices[0][d]) +
                                xiEtaZeta[1] * (vertices[2][d] - vertices[0][d]) +
                                xiEtaZeta[2] * (vertices[3][d] - vertices[0][d]);
        REQUIRE(nodes[p * Cell::Dim + d] ==
                AbsApprox(expected).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }
}

} // namespace seissol::unit_test
