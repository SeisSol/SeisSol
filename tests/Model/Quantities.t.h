// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Model/Common.h"
#include "Model/Quantities.h"

#include <cmath>
#include <array>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace seissol::unit_test {

namespace {

/// A right-handed orthonormal frame, built the way SeisSol's face normals are.
struct Frame {
  VrtxCoords normal{};
  VrtxCoords tangent1{};
  VrtxCoords tangent2{};
};

Frame randomFrame(std::mt19937& rng) {
  std::normal_distribution<double> gauss(0.0, 1.0);
  Frame frame{};
  double norm = 0.0;
  for (auto& component : frame.normal) {
    component = gauss(rng);
    norm += component * component;
  }
  norm = std::sqrt(norm);
  for (auto& component : frame.normal) {
    component /= norm;
  }

  const double aux[3] = {gauss(rng), gauss(rng), gauss(rng)};
  const double projection =
      aux[0] * frame.normal[0] + aux[1] * frame.normal[1] + aux[2] * frame.normal[2];
  double tangentNorm = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    frame.tangent1[i] = aux[i] - projection * frame.normal[i];
    tangentNorm += frame.tangent1[i] * frame.tangent1[i];
  }
  tangentNorm = std::sqrt(tangentNorm);
  for (auto& component : frame.tangent1) {
    component /= tangentNorm;
  }

  frame.tangent2[0] = frame.normal[1] * frame.tangent1[2] - frame.normal[2] * frame.tangent1[1];
  frame.tangent2[1] = frame.normal[2] * frame.tangent1[0] - frame.normal[0] * frame.tangent1[2];
  frame.tangent2[2] = frame.normal[0] * frame.tangent1[1] - frame.normal[1] * frame.tangent1[0];
  return frame;
}

/// True if (row, column) falls inside one of the declared diagonal blocks.
template <std::size_t N>
bool insideAGroup(const std::array<model::QuantityGroup, N>& groups,
                  std::size_t row,
                  std::size_t column) {
  std::size_t offset = 0;
  for (const auto& group : groups) {
    const auto end = offset + group.extent();
    if (row >= offset && row < end && column >= offset && column < end) {
      return true;
    }
    offset = end;
  }
  return false;
}

} // namespace

TEST_CASE("Quantity groups describe the configured material" * doctest::test_suite("model")) {
  constexpr auto Groups = model::MaterialT::RotationGroups;

  // Each declaration has to account for every quantity its matrix covers, and
  // may nominate at most one traction and one velocity group. Getting this
  // wrong is the failure mode the declaration exists to prevent, so it is a
  // compile-time error rather than a test assertion.
  //
  // T and Tinv are checked separately because they need not span the same
  // quantities: a solver keeping the mechanism index in its own tensor
  // dimension rotates one anelastic block forwards and none back.
  static_assert(model::quantitiesWellFormed(Groups, tensor::T::Shape[0]),
                "the quantity groups do not describe the rotation matrix");
  static_assert(model::quantitiesWellFormed(model::MaterialT::InverseRotationGroups,
                                            tensor::Tinv::Shape[0]),
                "the quantity groups do not describe the inverse rotation matrix");
  static_assert(tensor::T::Shape[0] == tensor::T::Shape[1]);
  static_assert(tensor::Tinv::Shape[0] == tensor::Tinv::Shape[1]);
  static_assert(tensor::Tinv::Shape[0] <= tensor::T::Shape[0]);

  // The velocity components start where the declaration says they do; the rest
  // of the code reaches for them through this offset.
  CHECK(model::roleOffset(Groups, model::FaceRole::Velocity) ==
        model::MaterialT::TractionQuantities);
}

TEST_CASE("Face rotation follows the quantity groups" * doctest::test_suite("model")) {
  constexpr double Epsilon = 1e4 * std::numeric_limits<real>::epsilon();
  constexpr std::size_t Size = tensor::T::Shape[0];
  constexpr std::size_t InverseSize = tensor::Tinv::Shape[0];

  std::vector<real> matTData(tensor::T::size());
  std::vector<real> matTinvData(tensor::Tinv::size());
  auto matT = init::T::view::create(matTData.data());
  auto matTinv = init::Tinv::view::create(matTinvData.data());

  std::mt19937 rng(20260904);
  for (int sample = 0; sample < 32; ++sample) {
    const auto frame = randomFrame(rng);
    model::getFaceRotationMatrix(frame.normal, frame.tangent1, frame.tangent2, matT, matTinv);

    SUBCASE("the inverse inverts the quantities it spans") {
      // This is what the two writers per kind are for: for a symmetric
      // second-order tensor the inverse is not the transpose, because the
      // Voigt weights differ. A block that only gets its forward writer
      // therefore shows up here as a missing identity.
      for (std::size_t i = 0; i < InverseSize; ++i) {
        for (std::size_t j = 0; j < InverseSize; ++j) {
          double accumulator = 0.0;
          for (std::size_t k = 0; k < InverseSize; ++k) {
            accumulator += static_cast<double>(matTinv(i, k)) * static_cast<double>(matT(k, j));
          }
          CHECK(std::abs(accumulator - (i == j ? 1.0 : 0.0)) < Epsilon);
        }
      }
    }

    SUBCASE("nothing is written outside the declared blocks") {
      for (std::size_t i = 0; i < Size; ++i) {
        for (std::size_t j = 0; j < Size; ++j) {
          if (!insideAGroup(model::MaterialT::RotationGroups, i, j)) {
            CHECK(matT(i, j) == static_cast<real>(0.0));
          }
        }
      }
      for (std::size_t i = 0; i < InverseSize; ++i) {
        for (std::size_t j = 0; j < InverseSize; ++j) {
          if (!insideAGroup(model::MaterialT::InverseRotationGroups, i, j)) {
            CHECK(matTinv(i, j) == static_cast<real>(0.0));
          }
        }
      }
    }

    SUBCASE("scalar groups do not rotate") {
      std::size_t offset = 0;
      for (const auto& group : model::MaterialT::RotationGroups) {
        if (group.kind == model::QuantityKind::Scalar) {
          CHECK(std::abs(static_cast<double>(matT(offset, offset)) - 1.0) < Epsilon);
          if (offset < InverseSize) {
            CHECK(std::abs(static_cast<double>(matTinv(offset, offset)) - 1.0) < Epsilon);
          }
        }
        offset += group.extent();
      }
    }
  }
}

} // namespace seissol::unit_test
