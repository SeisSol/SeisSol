// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Points.h"

#include <array>
#include <vector>

namespace {

std::vector<std::array<double, 1>> innerLine(int order) {
  std::vector<std::array<double, 1>> points;

  for (int i = 1; i < order; ++i) {
    points.emplace_back(std::array<double, 1>{i / static_cast<double>(order)});
  }

  return points;
}

template <std::size_t TargetDim, std::size_t SourceDim>
std::vector<std::array<double, TargetDim>>
    transform(const std::vector<std::array<double, SourceDim>>& points,
              const std::vector<std::vector<double>>& map,
              const std::vector<double>& offset) {
  std::vector<std::array<double, TargetDim>> output(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = 0; j < TargetDim; ++j) {
      for (std::size_t k = 0; k < SourceDim; ++k) {
        output[i][j] += map[j][k] * points[i][k];
      }
      output[i][j] += offset[j];
    }
  }
  return output;
}

template <std::size_t Dim, typename... Args>
std::vector<std::array<double, Dim>> concat(const Args&... args) {
  std::vector<std::array<double, Dim>> output;
  const auto addArg = [&](const auto& arg) { output.insert(output.end(), arg.begin(), arg.end()); };
  (addArg(args), ...);
  return output;
}

template <std::size_t Count, std::size_t Dim>
std::vector<std::array<double, Dim * Count>>
    cartesian(const std::vector<std::array<double, Dim>>& points) {
  if constexpr (Count == 0) {
    return {};
  } else if constexpr (Count == 1) {
    return points;
  } else {
    const auto subProduct = cartesian<Count - 1>(points);
    std::vector<std::array<double, Dim * Count>> product(subProduct.size() * points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      for (std::size_t j = 0; j < subProduct.size(); ++j) {
        for (std::size_t k = 0; k < Dim * (Count - 1); ++k) {
        }
      }
    }
    return product;
  }
}

template <std::size_t Dim>
std::vector<std::array<double, Dim>> simplex() {
  std::vector<std::array<double, Dim>> output(Dim + 1);

  // output[0] is implicitly all zero

  for (std::size_t i = 0; i < Dim; ++i) {
    output[i + 1][i] = 1;
  }

  return output;
}
} // namespace

namespace seissol::io::instance::geometry {
std::vector<std::array<double, 1>> pointsLine(int order) {
  if (order < 0) {
    return {};
  }
  if (order == 0) {
    return {{0.5}};
  }

  return concat<1>(simplex<1>(), innerLine(order));
}

std::vector<std::array<double, 2>> pointsTriangle(int order) {
  if (order < 0) {
    return {};
  }
  if (order == 0) {
    return {std::array<double, 2>{1. / 3, 1. / 3}};
  }

  const auto line = innerLine(order);
  const auto inner = pointsTriangle(order - 2);
  const auto offset = 1.0 / static_cast<double>(order);
  const auto scale = (order - 2) / static_cast<double>(order);

  return concat<2>(simplex<2>(),
                   transform<2>(line, {{1}, {0}}, {0, 0}),
                   transform<2>(line, {{0}, {1}}, {0, 0}),
                   transform<2>(line, {{-1}, {-1}}, {1, 1}),
                   transform<2>(inner, {{scale, 0}, {0, scale}}, {offset, offset}));
}

std::vector<std::array<double, 2>> pointsQuadrangle(int order) {
  return cartesian<2>(pointsLine(order));
}

std::vector<std::array<double, 3>> pointsHexahedron(int order) {
  return cartesian<3>(pointsLine(order));
}

std::vector<std::array<double, 3>> pointsTetrahedron(int order) {
  if (order < 0) {
    return {};
  }
  if (order == 0) {
    return {std::array<double, 3>{1. / 4, 1. / 4, 1. / 4}};
  }

  const auto line = innerLine(order);
  const auto face = pointsTriangle(order - 2);
  const auto inner = pointsTetrahedron(order - 3);
  const auto offsetFace = 1.0 / static_cast<double>(order);
  const auto scaleFace = (order - 2) / static_cast<double>(order);
  const auto offsetInner = 1.0 / static_cast<double>(order);
  const auto scaleInner = (order - 3) / static_cast<double>(order);

  return concat<3>(simplex<3>(),
                   transform<3>(line, {{1}, {0}, {0}}, {0, 0, 0}),
                   transform<3>(line, {{-1}, {1}, {0}}, {1, 0, 0}),
                   transform<3>(line, {{0}, {-1}, {0}}, {0, 1, 0}),
                   transform<3>(line, {{0}, {0}, {1}}, {0, 0, 0}),
                   transform<3>(line, {{-1}, {0}, {1}}, {1, 0, 0}),
                   transform<3>(line, {{0}, {-1}, {1}}, {0, 1, 0}),

                   transform<3>(face, {{1, 0}, {0, 0}, {0, 1}}, {0, 0, 0}),
                   transform<3>(face, {{0, 1}, {-1, -1}, {1, 0}}, {0, 1, 0}),
                   transform<3>(face, {{0, 0}, {0, 1}, {1, 0}}, {0, 0, 0}),
                   transform<3>(face, {{0, 1}, {1, 0}, {0, 0}}, {0, 0, 0}),

                   transform<3>(inner,
                                {{scaleInner, 0, 0}, {0, scaleInner, 0}, {0, 0, scaleInner}},
                                {offsetInner, offsetInner, offsetInner}));
}
} // namespace seissol::io::instance::geometry
