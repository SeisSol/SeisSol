// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Common/Constants.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "Numerical/Functions.h"
#include "Numerical/Projection.h"
#include "Numerical/Transformation.h"
#include "Solver/MultipleSimulations.h"
#include "TestHelper.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::unit_test {

namespace {

using namespace seissol::numerical;
namespace projection = seissol::numerical::projection;

constexpr double Tolerance = 1e-10;

// For fused simulations the code generator transposes everything in the `nodal` namespace
// (cf. kernels/aderdg.py); detect that from a matrix whose shape is not square.
bool nodalTransposed() {
  return static_cast<std::size_t>(nodal::tensor::nodes2D::Shape[0]) !=
         projection::modalSize(2, ConvergenceOrder);
}

// The nodal point set of the volume follows the PLASTICITY_METHOD build option: "nb" ships a
// unisolvent warp&blend set, "ip" the conical-product quadrature points.
constexpr auto VolumeNodalSet =
    static_cast<std::size_t>(tensor::vNodes::Shape[0]) == projection::modalSize(3, ConvergenceOrder)
        ? projection::NodalSet::WarpBlend
        : projection::NodalSet::Stroud;

using Tet = std::array<std::array<double, 3>, 4>;

//! @brief The reference tetrahedron, i.e. element.vertices[0..3] in reference coordinates.
constexpr Tet ReferenceTet{std::array<double, 3>{0, 0, 0},
                           std::array<double, 3>{1, 0, 0},
                           std::array<double, 3>{0, 1, 0},
                           std::array<double, 3>{0, 0, 1}};

std::array<double, 3> midpoint(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
}

// Transcriptions of the legacy refinement::DivideTetrahedronBy4 and ...By8. The legacy
// refinement::MeshRefiner passed element.vertices[0..3] as a, b, c, d, so evaluating these on
// ReferenceTet gives the subcells the refined output mesh used to consist of.
std::vector<Tet> legacyRefine4(const Tet& in) {
  const auto& [a, b, c, d] = in;
  const auto center = midpoint(midpoint(a, b), midpoint(c, d));
  return {Tet{a, b, c, center}, Tet{a, b, d, center}, Tet{a, c, d, center}, Tet{b, c, d, center}};
}

std::vector<Tet> legacyRefine8(const Tet& in) {
  const auto& [a, b, c, d] = in;
  const auto ab = midpoint(a, b);
  const auto ac = midpoint(a, c);
  const auto ad = midpoint(a, d);
  const auto bc = midpoint(b, c);
  const auto bd = midpoint(b, d);
  const auto cd = midpoint(c, d);
  return {Tet{a, ab, ac, ad},
          Tet{b, ab, bc, bd},
          Tet{c, ac, bc, cd},
          Tet{d, ad, bd, cd},
          Tet{ab, ac, ad, bd},
          Tet{ab, ac, bc, bd},
          Tet{ac, ad, bd, cd},
          Tet{ac, bc, bd, cd}};
}

//! @brief Swaps the last two vertices of a negatively oriented tetrahedron, as the tables do.
Tet orientedTet(const Tet& in) {
  const std::array<double, 3> u{in[1][0] - in[0][0], in[1][1] - in[0][1], in[1][2] - in[0][2]};
  const std::array<double, 3> v{in[2][0] - in[0][0], in[2][1] - in[0][1], in[2][2] - in[0][2]};
  const std::array<double, 3> w{in[3][0] - in[0][0], in[3][1] - in[0][1], in[3][2] - in[0][2]};
  const auto determinant = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
                           u[2] * (v[0] * w[1] - v[1] * w[0]);
  return determinant > 0 ? in : Tet{in[0], in[1], in[3], in[2]};
}

Tet subcellVertices(const AffineMap<3, 3>& map) {
  Tet vertices{};
  for (std::size_t d = 0; d < 3; ++d) {
    vertices[0][d] = map.offset[d];
    for (std::size_t j = 0; j < 3; ++j) {
      vertices[j + 1][d] = map.offset[d] + map.matrix[d][j];
    }
  }
  return vertices;
}

double signedVolume(const AffineMap<3, 3>& map) {
  const auto& m = map.matrix;
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Compares two tetrahedra as vertex sets, i.e. up to the vertex order inside the cell. All
// coordinates involved are dyadic rationals, so sorting them exactly is safe.
void requireSameCell(const Tet& actual, const Tet& expected) {
  auto lhs = actual;
  auto rhs = expected;
  std::sort(lhs.begin(), lhs.end());
  std::sort(rhs.begin(), rhs.end());
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t d = 0; d < 3; ++d) {
      REQUIRE(lhs[i][d] == AbsApprox(rhs[i][d]).epsilon(Tolerance).delta(Tolerance));
    }
  }
}

seissol::numerical::AffineMap<2, 3> faceEmbedding(std::size_t side) {
  const std::array<std::array<double, 2>, 3> corners = {
      std::array<double, 2>{0, 0}, std::array<double, 2>{1, 0}, std::array<double, 2>{0, 1}};
  std::vector<std::array<double, 3>> vertices;
  for (const auto& chiTau : corners) {
    std::array<double, 3> xez{};
    seissol::transformations::chiTau2XiEtaZeta(
        static_cast<std::uint32_t>(side), chiTau.data(), xez.data());
    vertices.push_back(xez);
  }
  return seissol::numerical::AffineMap<2, 3>::fromVertices(vertices);
}

} // namespace

// The projection module re-derives, at run time, what the code generator ships as constants. These
// tests pin that agreement down; if they fail, the run-time output matrices and the generated
// kernels disagree about the basis or nodal conventions.
TEST_CASE("Numerical/Projection: nodal point sets match the generated ones") {
  SUBCASE("2D face nodes vs. nodes2D") {
    const auto points = projection::nodalPoints2D(ConvergenceOrder);
    const auto nodes = nodal::init::nodes2D::view::create(nodal::init::nodes2D::Values);
    const auto transposed = nodalTransposed();
    for (std::size_t p = 0; p < points.size(); ++p) {
      for (std::size_t d = 0; d < 2; ++d) {
        const auto i = transposed ? d : p;
        const auto j = transposed ? p : d;
        const auto reference = nodes.isInRange(i, j) ? nodes(i, j) : 0.0;
        REQUIRE(points[p][d] == AbsApprox(reference).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }

  SUBCASE("3D volume nodes vs. vNodes") {
    const auto points = projection::nodalPoints3D(ConvergenceOrder, VolumeNodalSet);
    const auto nodes = init::vNodes::view::create(init::vNodes::Values);
    REQUIRE(points.size() == static_cast<std::size_t>(tensor::vNodes::Shape[0]));
    for (std::size_t p = 0; p < points.size(); ++p) {
      for (std::size_t d = 0; d < 3; ++d) {
        const auto reference = nodes.isInRange(p, d) ? nodes(p, d) : 0.0;
        REQUIRE(points[p][d] == AbsApprox(reference).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }
}

TEST_CASE("Numerical/Projection: nodal-to-modal transforms match the generated ones") {
  SUBCASE("2D vs. MV2nTo2m") {
    // MV2nTo2m is stored as [modalBasis][node]
    const auto matrix =
        projection::nodalToModal<2>(ConvergenceOrder, projection::NodalSet::WarpBlend);
    const auto reference = nodal::init::MV2nTo2m::view::create(nodal::init::MV2nTo2m::Values);
    const auto transposed = nodalTransposed();
    for (std::size_t b = 0; b < matrix.rows(); ++b) {
      for (std::size_t n = 0; n < matrix.cols(); ++n) {
        const auto i = transposed ? n : b;
        const auto j = transposed ? b : n;
        const auto expected = reference.isInRange(i, j) ? reference(i, j) : 0.0;
        REQUIRE(matrix(b, n) == AbsApprox(expected).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }

  SUBCASE("3D vs. vInv") {
    const auto matrix = projection::nodalToModal<3>(ConvergenceOrder, VolumeNodalSet);
    const auto reference = init::vInv::view::create(init::vInv::Values);
    for (std::size_t b = 0; b < matrix.rows(); ++b) {
      for (std::size_t n = 0; n < matrix.cols(); ++n) {
        const auto expected = reference.isInRange(b, n) ? reference(b, n) : 0.0;
        REQUIRE(matrix(b, n) == AbsApprox(expected).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }
}

TEST_CASE("Numerical/Projection: both volume nodal sets round-trip") {
  // Only one of the two is shipped by the code generator in a given build, so exercise both
  // against the definition instead: sampling a modal function at the nodes and transforming back
  // has to reproduce the coefficients.
  const auto indices = projection::modalIndices<3>(ConvergenceOrder);

  for (const auto set : {projection::NodalSet::WarpBlend, projection::NodalSet::Stroud}) {
    const auto points = projection::nodalPoints3D(ConvergenceOrder, set);
    const auto matrix = projection::nodalToModal<3>(ConvergenceOrder, set);

    REQUIRE(points.size() == projection::nodalSize(3, ConvergenceOrder, set));
    REQUIRE(matrix.rows() == indices.size());
    REQUIRE(matrix.cols() == points.size());

    for (const auto& point : points) {
      const auto sum = point[0] + point[1] + point[2];
      REQUIRE(point[0] >= -Tolerance);
      REQUIRE(point[1] >= -Tolerance);
      REQUIRE(point[2] >= -Tolerance);
      REQUIRE(sum <= 1 + Tolerance);
    }

    for (std::size_t coefficient = 0; coefficient < indices.size(); ++coefficient) {
      std::vector<double> nodalValues(points.size());
      for (std::size_t n = 0; n < points.size(); ++n) {
        nodalValues[n] = functions::DubinerP<3>(indices[coefficient], points[n]);
      }
      for (std::size_t b = 0; b < indices.size(); ++b) {
        double value = 0;
        for (std::size_t n = 0; n < points.size(); ++n) {
          value += matrix(b, n) * nodalValues[n];
        }
        const auto expected = b == coefficient ? 1.0 : 0.0;
        REQUIRE(value == AbsApprox(expected).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }
}

TEST_CASE("Numerical/Projection: L2 projection reproduces interpolation when exact") {
  // If the target space contains the source space, the L2 projection and the pointwise evaluation
  // have to agree -- this exercises the quadrature, the mass matrix and the Vandermonde inversion
  // in one go.
  const auto degree = ConvergenceOrder - 1;
  const auto points = io::instance::geometry::pointsTetrahedron(static_cast<int>(degree));

  AffineMap<3, 3> subcell;
  subcell.offset = {0.25, 0.125, 0.5};
  subcell.matrix = {{{0.5, 0, 0}, {0, 0.5, 0}, {0, 0, 0.5}}};

  projection::Spec interpolate;
  interpolate.order = ConvergenceOrder;
  auto project = interpolate;
  project.target = projection::Target::Project;

  const auto a = projection::build<3, 3>(points, degree, subcell, interpolate);
  const auto b = projection::build<3, 3>(points, degree, subcell, project);

  REQUIRE(a.rows() == b.rows());
  REQUIRE(a.cols() == b.cols());
  for (std::size_t r = 0; r < a.rows(); ++r) {
    for (std::size_t c = 0; c < a.cols(); ++c) {
      REQUIRE(b(r, c) == AbsApprox(a(r, c)).epsilon(1e-9).delta(1e-9));
    }
  }
}

TEST_CASE("Numerical/Projection: derivatives match finite differences") {
  constexpr std::size_t Degree = 2;
  const auto points = io::instance::geometry::pointsTetrahedron(Degree);

  AffineMap<3, 3> subcell;
  subcell.offset = {0.1, 0.2, 0.05};
  subcell.matrix = {{{0.25, 0.125, 0}, {0, 0.25, 0.0625}, {0, 0, 0.25}}};

  projection::Spec spec;
  spec.order = ConvergenceOrder;

  const double h = 1e-5;
  for (std::size_t direction = 0; direction < Cell::Dim; ++direction) {
    auto derivativeSpec = spec;
    derivativeSpec.derivative = direction;
    const auto derivative = projection::build<3, 3>(points, Degree, subcell, derivativeSpec);

    auto plus = subcell;
    auto minus = subcell;
    plus.offset[direction] += h;
    minus.offset[direction] -= h;
    const auto valuePlus = projection::build<3, 3>(points, Degree, plus, spec);
    const auto valueMinus = projection::build<3, 3>(points, Degree, minus, spec);

    for (std::size_t r = 0; r < derivative.rows(); ++r) {
      for (std::size_t c = 0; c < derivative.cols(); ++c) {
        const auto finiteDifference = (valuePlus(r, c) - valueMinus(r, c)) / (2 * h);
        REQUIRE(derivative(r, c) == AbsApprox(finiteDifference).epsilon(1e-5).delta(1e-5));
      }
    }
  }
}

TEST_CASE("Numerical/Projection: the table matches build() in the generated layout") {
  constexpr std::size_t Degree = 3;
  const auto points = io::instance::geometry::pointsTetrahedron(Degree);

  auto subcells = io::instance::geometry::unrefined<3>();
  subcells =
      io::instance::geometry::subdivideMaps(subcells, io::instance::geometry::TetrahedronRefine4);

  const auto index = tensor::collvv::index(ConvergenceOrder, Degree);
  const std::size_t stride = tensor::collvv::Size[index] / tensor::collvv::Shape[index][1];
  REQUIRE(stride >= points.size());

  projection::Spec spec;
  const projection::Table<3, 3> table(subcells, points, Degree, stride, spec, 1, ConvergenceOrder);
  REQUIRE(table.subcellCount() == subcells.size());

  for (std::size_t subcell = 0; subcell < subcells.size(); ++subcell) {
    for (std::size_t order = 1; order <= ConvergenceOrder; ++order) {
      auto local = spec;
      local.order = order;
      const auto reference = projection::build<3, 3>(points, Degree, subcells[subcell], local);
      const auto* stored = table(subcell, order);
      for (std::size_t p = 0; p < reference.rows(); ++p) {
        for (std::size_t b = 0; b < reference.cols(); ++b) {
          REQUIRE(
              stored[p + b * stride] ==
              AbsApprox(static_cast<real>(reference(p, b))).epsilon(Tolerance).delta(Tolerance));
        }
      }
    }
  }
}

TEST_CASE("Numerical/Projection: subcells tile the reference cell without inverting") {
  // The union of the subcell volumes has to be the reference cell. Note that the volumes are
  // summed signed on purpose: an inverted subcell would still tile the reference cell, but it
  // would be written out as a tetrahedron with a negative Jacobian.
  for (const auto& refine :
       {io::instance::geometry::TetrahedronRefine4, io::instance::geometry::TetrahedronRefine8}) {
    const auto subcells =
        io::instance::geometry::subdivideMaps(io::instance::geometry::unrefined<3>(), refine);
    double total = 0;
    for (const auto& subcell : subcells) {
      REQUIRE(signedVolume(subcell) > Tolerance);
      total += signedVolume(subcell);
    }
    REQUIRE(total == AbsApprox(1.0).epsilon(Tolerance).delta(Tolerance));
  }

  const auto faceRefine = io::instance::geometry::subdivideMaps(
      io::instance::geometry::unrefined<2>(), io::instance::geometry::TriangleRefine4);
  double area = 0;
  for (const auto& subcell : faceRefine) {
    const auto& m = subcell.matrix;
    const auto subarea = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    REQUIRE(subarea > Tolerance);
    area += subarea;
  }
  REQUIRE(area == AbsApprox(1.0).epsilon(Tolerance).delta(Tolerance));
}

// The subdivision tables have to be read in the vertex order of the reference simplex that
// AffineMap::fromVertices reconstructs from -- not in the (rotated) vertex order of the legacy
// refinement::Tetrahedron::unitTetrahedron(). Getting that wrong is invisible in the volumes:
// the relabelling is an even permutation, so the subcells still tile the reference cell, they
// just pick a different diagonal of the inner octahedron.
TEST_CASE("Numerical/Projection: tetrahedron subdivisions match the legacy refiner") {
  const auto check = [](const std::vector<AffineMap<3, 3>>& subcells,
                        const std::vector<Tet>& legacy) {
    REQUIRE(subcells.size() == legacy.size());
    for (std::size_t i = 0; i < subcells.size(); ++i) {
      requireSameCell(subcellVertices(subcells[i]), legacy[i]);
    }
  };

  const auto refine4 = io::instance::geometry::subdivideMaps(
      io::instance::geometry::unrefined<3>(), io::instance::geometry::TetrahedronRefine4);
  check(refine4, legacyRefine4(ReferenceTet));

  const auto refine8 = io::instance::geometry::subdivideMaps(
      io::instance::geometry::unrefined<3>(), io::instance::geometry::TetrahedronRefine8);
  check(refine8, legacyRefine8(ReferenceTet));

  // refinement::DivideTetrahedronBy32 divided by 8 first and split each of those subcells by
  // its center point afterwards, in the ordering 4*i + j. The center split enumerates the four
  // faces of the cell it subdivides, so it depends on that cell's vertex order: for the four
  // outer subcells the legacy By8 table had inverted, reorienting them swaps the first two of
  // the four cells the split produces. That is unavoidable -- no odd permutation leaves the
  // enumeration fixed -- and it only renumbers cells within a block of four; the resulting set
  // of 32 cells is the legacy one either way.
  std::vector<Tet> legacy32;
  for (const auto& outer : legacyRefine8(ReferenceTet)) {
    for (const auto& inner : legacyRefine4(orientedTet(outer))) {
      legacy32.push_back(inner);
    }
  }
  const auto refine32 =
      io::instance::geometry::subdivideMaps(refine8, io::instance::geometry::TetrahedronRefine4);
  check(refine32, legacy32);
}

TEST_CASE("Numerical/Projection: face embedding is consistent with chiTau2XiEtaZeta") {
  const auto points = io::instance::geometry::pointsTriangle(4);
  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    const auto embedding = faceEmbedding(side);
    for (const auto& point : points) {
      std::array<double, 3> expected{};
      seissol::transformations::chiTau2XiEtaZeta(
          static_cast<std::uint32_t>(side), point.data(), expected.data());
      const auto actual = embedding(point);
      for (std::size_t d = 0; d < 3; ++d) {
        REQUIRE(actual[d] == AbsApprox(expected[d]).epsilon(Tolerance).delta(Tolerance));
      }
    }
  }
}

} // namespace seissol::unit_test
