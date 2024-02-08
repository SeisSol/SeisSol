#ifndef SEISSOL_DR_OUTPUT_GEOMETRY_HPP
#define SEISSOL_DR_OUTPUT_GEOMETRY_HPP

#include "Geometry/MeshDefinition.h"
#include "Kernels/precision.hpp"
#include <Eigen/Dense>
#include <array>
#include <cassert>

namespace seissol::dr {
struct ExtVrtxCoords {
  ExtVrtxCoords() = default;
  ~ExtVrtxCoords() = default;

  template <typename T>
  ExtVrtxCoords(const T& other) {
    for (int i = 0; i < 3; ++i)
      coords[i] = other[i];
  }

  template <typename T>
  ExtVrtxCoords& operator=(const T& other) {
    for (int i = 0; i < 3; ++i)
      coords[i] = other[i];
    return *this;
  }

  ExtVrtxCoords(std::initializer_list<double> inputCoords) {
    assert(inputCoords.size() == 3 && "ExtVrtxCoords must get initialized with 3 values");
    const auto* begin = inputCoords.begin();
    for (int i = 0; i < 3; ++i, ++begin)
      coords[i] = *begin;
  }

  double& operator[](size_t index) {
    assert((index < 3) && "ExtVrtxCoords index must be less than 3");
    return coords[index];
  }

  double operator[](size_t index) const {
    assert((index < 3) && "ExtVrtxCoords index must be less than 3");
    return coords[index];
  }

  Eigen::Vector3d getAsEigen3LibVector() const {
    return Eigen::Vector3d(coords[0], coords[1], coords[2]);
  }

  constexpr static int size() { return 3; }

  VrtxCoords coords = {0.0, 0.0, 0.0};
};

struct ExtTriangle {
  ExtTriangle() = default;
  ~ExtTriangle() = default;
  explicit ExtTriangle(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2, const ExtVrtxCoords& p3) {
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
  }

  ExtTriangle(const ExtTriangle& other) {
    for (int i = 0; i < 3; ++i)
      points[i] = other.points[i];
  }

  ExtTriangle& operator=(const ExtTriangle& other) {
    for (int i = 0; i < 3; ++i)
      points[i] = other.points[i];
    return *this;
  }

  ExtVrtxCoords& point(size_t index) {
    assert((index < 3) && "ExtTriangle index must be less than 3");
    return points[index];
  }

  ExtVrtxCoords point(size_t index) const {
    assert((index < 3) && "ExtTriangle index must be less than 3");
    return points[index];
  }

  static int size() { return 3; }

  private:
  std::array<ExtVrtxCoords, 3> points{};
};

struct ReceiverPoint {
  ExtVrtxCoords global{};       // physical coords of a receiver
  ExtVrtxCoords reference{};    // reference coords of a receiver
  ExtTriangle globalTriangle{}; // a surrounding triangle of a receiver
  int faultFaceIndex{-1};       // Face Fault index which the receiver belongs to
  int localFaceSideId{-1};      // Side ID of a reference element
  int elementIndex{-1};         // Element which the receiver belongs to
  int globalReceiverIndex{-1};  // receiver index of global list
  bool isInside{false};         // If a point is inside the mesh or not
  int nearestGpIndex{-1};
  int faultTag{-1};

  // Internal points are required because computed gradients
  // are inaccurate near triangle edges,
  // specifically for low-order elements
  int nearestInternalGpIndex{-1};
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

#endif // SEISSOL_DR_OUTPUT_GEOMETRY_HPP
