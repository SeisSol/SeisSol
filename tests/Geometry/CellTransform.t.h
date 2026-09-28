// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <random>
#include <vector>

namespace seissol::unit_test {

namespace {
using seissol::geometry::AffineTransform;
using seissol::geometry::CellTransform;
using VectorT = CellTransform::VectorEigenT;

/**
 * A quadratic cell transform. It is not a transform SeisSol uses; it exists so that the generic
 * parts of CellTransform are exercised on something whose Jacobian actually varies, which an
 * AffineTransform cannot do.
 */
class QuadraticTransform : public CellTransform {
  public:
  [[nodiscard]] auto refToSpace(const VectorT& input) const -> VectorT override {
    return {1.0 + 2.0 * input(0) + 0.3 * input(1) * input(1),
            2.0 + 1.5 * input(1) + 0.2 * input(2) * input(2),
            3.0 + 1.0 * input(2) + 0.1 * input(0) * input(0)};
  }

  [[nodiscard]] auto refToSpaceJacobian(const VectorT& input) const -> MatrixEigenT override {
    MatrixEigenT jacobian;
    jacobian << 2.0, 0.6 * input(1), 0.0, 0.0, 1.5, 0.4 * input(2), 0.2 * input(0), 0.0, 1.0;
    return jacobian;
  }
};

auto randomTetrahedron(std::mt19937& rng, double scale, const VectorT& shift)
    -> std::array<VectorT, Cell::NumVertices> {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::array<VectorT, Cell::NumVertices> vertices{};
  // rejection sampling: a random tetrahedron may be arbitrarily flat, which would make the round
  // trip below fail for reasons of conditioning rather than correctness
  while (true) {
    for (auto& vertex : vertices) {
      vertex = VectorT(
          shift(0) + scale * dist(rng), shift(1) + scale * dist(rng), shift(2) + scale * dist(rng));
    }
    const auto transform = AffineTransform(vertices);
    if (std::abs(transform.determinant()) > 1e-3 * scale * scale * scale) {
      return vertices;
    }
  }
}
} // namespace

TEST_CASE("Affine cell transform" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-12;
  std::mt19937 rng(20260928);

  SUBCASE("Maps the reference vertices onto the cell vertices") {
    const auto vertices = randomTetrahedron(rng, 1.0, VectorT::Zero());
    const auto transform = AffineTransform(vertices);

    const std::array<VectorT, Cell::NumVertices> reference{
        VectorT(0, 0, 0), VectorT(1, 0, 0), VectorT(0, 1, 0), VectorT(0, 0, 1)};

    for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
      const auto mapped = transform.refToSpace(reference[i]);
      REQUIRE((mapped - vertices[i]).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("Round trips reference to space and back") {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int sample = 0; sample < 100; ++sample) {
      const auto vertices = randomTetrahedron(rng, 1.0, VectorT::Zero());
      const auto transform = AffineTransform(vertices);

      const auto point = VectorT(0.2 * dist(rng), 0.2 * dist(rng), 0.2 * dist(rng));
      const auto back = transform.spaceToRef(transform.refToSpace(point));
      REQUIRE((back - point).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("The determinant is six times the cell volume") {
    for (int sample = 0; sample < 100; ++sample) {
      const auto vertices = randomTetrahedron(rng, 1.0, VectorT::Zero());
      const auto transform = AffineTransform(vertices);

      const auto volume = std::abs((vertices[1] - vertices[0])
                                       .cross(vertices[2] - vertices[0])
                                       .dot(vertices[3] - vertices[0])) /
                          6.0;
      REQUIRE(std::abs(transform.determinant()) / 6.0 == AbsApprox(volume).epsilon(Epsilon));
    }
  }

  SUBCASE("The inverse Jacobian inverts the Jacobian") {
    const auto vertices = randomTetrahedron(rng, 1.0, VectorT::Zero());
    const auto transform = AffineTransform(vertices);

    const auto point = VectorT(0.1, 0.2, 0.3);
    const auto product =
        transform.refToSpaceJacobian(point) * transform.refToSpaceJacobianInverse(point);
    REQUIRE((product - CellTransform::MatrixEigenT::Identity()).norm() ==
            AbsApprox(0.0).epsilon(Epsilon));
  }

  SUBCASE("Stays accurate for cells far away from the origin") {
    // production meshes are commonly given in metres with a geographic offset; a cell is then
    // small compared to its own coordinates, and every intermediate that is not differenced first
    // loses most of its significant digits
    const auto shift = VectorT(5e5, 5e5, 0.0);
    for (int sample = 0; sample < 100; ++sample) {
      const auto vertices = randomTetrahedron(rng, 1e2, shift);
      const auto transform = AffineTransform(vertices);

      const auto point = VectorT(0.1, 0.2, 0.3);
      const auto back = transform.spaceToRef(transform.refToSpace(point));
      REQUIRE((back - point).norm() == AbsApprox(0.0).epsilon(1e-10));

      const auto product =
          transform.refToSpaceJacobian(point) * transform.refToSpaceJacobianInverse(point);
      REQUIRE((product - CellTransform::MatrixEigenT::Identity()).norm() ==
              AbsApprox(0.0).epsilon(1e-10));
    }
  }

  SUBCASE("The batched overload agrees with the scalar one") {
    const auto vertices = randomTetrahedron(rng, 1.0, VectorT::Zero());
    const auto transform = AffineTransform(vertices);

    const std::vector<CellTransform::VectorT> points{
        {0.1, 0.2, 0.3}, {0.0, 0.0, 0.0}, {0.25, 0.25, 0.25}};
    const auto batched = transform.refToSpace(points);

    REQUIRE(batched.size() == points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      const auto single = transform.refToSpace(points[i]);
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        REQUIRE(batched[i][d] == AbsApprox(single[d]).epsilon(Epsilon));
      }
    }
  }
}

TEST_CASE("Generic cell transform inversion" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-9;
  const QuadraticTransform transform;

  SUBCASE("The Newton fallback inverts a non-affine transform") {
    const std::array<VectorT, 3> points{
        VectorT(0.2, 0.3, 0.4), VectorT(0.0, 0.0, 0.0), VectorT(0.25, 0.25, 0.25)};
    for (const auto& point : points) {
      const auto target = transform.refToSpace(point);
      const auto recovered = transform.spaceToRef(target);
      REQUIRE((recovered - point).norm() == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE((transform.refToSpace(recovered) - target).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("The space-coordinate Jacobian agrees with the reference-coordinate one") {
    // the two differ for a non-affine transform, and must agree once the same point is named in
    // both coordinate systems
    const auto point = VectorT(0.2, 0.3, 0.4);
    const auto space = transform.refToSpace(point);
    const auto fromRef = transform.refToSpaceJacobianInverse(point);
    const auto fromSpace = transform.spaceToRefJacobian(space);
    REQUIRE((fromRef - fromSpace).norm() == AbsApprox(0.0).epsilon(Epsilon));
  }
}

} // namespace seissol::unit_test
