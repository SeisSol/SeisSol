// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_TYPEDEFS_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_TYPEDEFS_H_

#include <cstddef>
#include <string>
#include <utils/logger.h>
namespace seissol::io::instance::geometry {

enum class Shape {
  Line,
  Triangle,
  Tetrahedron,
  Quadrangle,
  Hexahedron,

  // NYI:
  // Prism,
  // Pyramid
};

constexpr std::size_t numPoints(std::size_t order, Shape shape) {
  order = std::max(static_cast<std::size_t>(1), order);
  switch (shape) {
  case Shape::Line:
    return order + 1;
  case Shape::Triangle:
    return ((order + 1) * (order + 2)) / 2;
  case Shape::Tetrahedron:
    return ((order + 1) * (order + 2) * (order + 3)) / 6;
  case Shape::Quadrangle:
    return (order + 1) * (order + 1);
  case Shape::Hexahedron:
    return (order + 1) * (order + 1) * (order + 1);
  }
}

constexpr std::size_t vtkType(Shape shape) {
  switch (shape) {
  case Shape::Line:
    return 1;
  case Shape::Triangle:
    return 69;
  case Shape::Tetrahedron:
    return 71;
  case Shape::Quadrangle:
    return 1;
  case Shape::Hexahedron:
    return 1;
  }
}

inline std::string xdmfType(Shape shape, std::size_t order) {
  if (order > 1) {
    logError() << "High-order output is not supported for Xdmf.";
  }
  switch (shape) {
  case Shape::Line:
    return "Polyline";
  case Shape::Triangle:
    return "Triangle";
  case Shape::Tetrahedron:
    return "Tetrahedron";
  case Shape::Quadrangle:
    return "Quadrilateral";
  case Shape::Hexahedron:
    return "Hexahedron";
  }
}

constexpr std::size_t dimension(Shape shape) {
  switch (shape) {
  case Shape::Line:
    return 1;
  case Shape::Triangle:
    return 2;
  case Shape::Tetrahedron:
    return 3;
  case Shape::Quadrangle:
    return 2;
  case Shape::Hexahedron:
    return 3;
  }
}

} // namespace seissol::io::instance::geometry
#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_TYPEDEFS_H_
