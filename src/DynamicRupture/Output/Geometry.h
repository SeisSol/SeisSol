// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_GEOMETRY_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_GEOMETRY_H_

#include "Common/CompactOptional.h"
#include "Geometry/MeshDefinition.h"

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace seissol::dr {
struct ExtTriangle {
  ExtTriangle() = default;

  explicit ExtTriangle(const CoordinateT& p1, const CoordinateT& p2, const CoordinateT& p3) {
    points_[0] = p1;
    points_[1] = p2;
    points_[2] = p3;
  }

  CoordinateT& point(size_t index) {
    assert((index < points_.size()) && "ExtTriangle index must be less than 3");
    return points_[index];
  }

  [[nodiscard]] const CoordinateT& point(size_t index) const {
    assert((index < points_.size()) && "ExtTriangle index must be less than 3");
    return points_[index];
  }

  static constexpr std::size_t size() { return Size; }

  private:
  static constexpr std::size_t Size = 3;
  std::array<CoordinateT, Size> points_{};
};

struct Receiver {
  // physical coords of a receiver
  CoordinateT global{};

  // reference coords of a receiver
  CoordinateT reference{};

  // a surrounding triangle of a receiver
  ExtTriangle globalTriangle;

  // Face Fault index which the receiver belongs to
  OptionalSize faultFaceIndex;

  // Side ID of a reference element
  OptionalSide localFaceSideId;

  // Side ID (minus) of a reference element
  OptionalSide localNeighborFaceSideId;

  // Rank-local element ID to which the receiver belongs
  OptionalSize elementIndex;

  // Global element ID to which the receiver belongs
  OptionalSize elementGlobalIndex;

  // Global element ID (minus) to which the receiver belongs
  OptionalSize elementNeighborGlobalIndex;

  // receiver index of global list
  OptionalSize globalReceiverIndex;

  // If a point is inside the mesh or not
  bool isInside{false};

  int nearestGpIndex{-1};

  int faultTag{-1};

  // Simulation index for multisim
  int simIndex{0};

  // Index of the nearest quadrature point considering fused simulations
  int gpIndex{-1};

  // Internal points are required because computed gradients
  // are inaccurate near triangle edges,
  // specifically for low-order elements
  int nearestInternalGpIndex{-1};

  // Index of the nearest internal quadrature point considering fused simulations
  int internalGpIndexFused{-1};

  [[nodiscard]] constexpr std::size_t globalFaultFaceId() const {
    return elementGlobalIndex.value() * Cell::NumFaces +
           static_cast<std::size_t>(localFaceSideId.value());
  }
};
using Receivers = std::vector<Receiver>;

struct FaultDirections {
  std::array<double, 3> faceNormal{};
  std::array<double, 3> tangent1{};
  std::array<double, 3> tangent2{};
  std::array<double, 3> strike{};
  std::array<double, 3> dip{};
};
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_GEOMETRY_H_
