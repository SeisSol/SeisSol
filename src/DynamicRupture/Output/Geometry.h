// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_GEOMETRY_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_GEOMETRY_H_

#include "Geometry/MeshDefinition.h"
#include "Kernels/Precision.h"

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <limits>

namespace seissol::dr {
struct ExtTriangle {
  ExtTriangle() = default;

  explicit ExtTriangle(const CoordinateT& p1, const CoordinateT& p2, const CoordinateT& p3) {
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
  }

  CoordinateT& point(size_t index) {
    assert((index < points.size()) && "ExtTriangle index must be less than 3");
    return points[index];
  }

  [[nodiscard]] const CoordinateT& point(size_t index) const {
    assert((index < points.size()) && "ExtTriangle index must be less than 3");
    return points[index];
  }

  static constexpr std::size_t size() { return Size; }

  private:
  static constexpr std::size_t Size = 3;
  std::array<CoordinateT, Size> points{};
};

struct ReceiverPoint {
  CoordinateT global{};       // physical coords of a receiver
  CoordinateT reference{};    // reference coords of a receiver
  ExtTriangle globalTriangle; // a surrounding triangle of a receiver
  size_t faultFaceIndex{
      std::numeric_limits<std::size_t>::max()}; // Face Fault index which the receiver belongs to
  int8_t localFaceSideId{-1};                   // Side ID of a reference element
  size_t elementIndex{
      std::numeric_limits<std::size_t>::max()}; // Element which the receiver belongs to
  std::size_t elementGlobalIndex{
      std::numeric_limits<std::size_t>::max()}; // Element which the receiver belongs to
  size_t globalReceiverIndex{
      std::numeric_limits<std::size_t>::max()}; // receiver index of global list
  bool isInside{false};                         // If a point is inside the mesh or not
  int nearestGpIndex{-1};
  int faultTag{-1};
  int simIndex{0}; // Simulation index for multisim
  int gpIndex{-1}; // Index of the nearest gaussian point considering fused simulations

  // Internal points are required because computed gradients
  // are inaccurate near triangle edges,
  // specifically for low-order elements
  int nearestInternalGpIndex{-1};
  int internalGpIndexFused{
      -1}; // Index of the nearest internal gaussian point considering fused simulations
};
using ReceiverPoints = std::vector<ReceiverPoint>;

struct FaultDirections {
  std::array<double, 3> faceNormal{};
  std::array<double, 3> tangent1{};
  std::array<double, 3> tangent2{};
  std::array<double, 3> strike{};
  std::array<double, 3> dip{};
};
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_GEOMETRY_H_
