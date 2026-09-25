// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Projection.h"

#include "Kernels/Precision.h"
#include "Numerical/Functions.h"
#include "Numerical/Quadrature.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace {

using namespace seissol;
using namespace seissol::numerical;

/**
 * Quadrature on the Dim-dimensional reference simplex, exact for polynomial degree 2n - 1.
 */
template <std::size_t Dim>
std::pair<std::vector<std::array<double, Dim>>, std::vector<double>>
    simplexQuadrature(std::size_t n) {
  std::vector<std::array<double, Dim>> points;
  std::vector<double> weights;
  if constexpr (Dim == 1) {
    auto line = quadrature::ShiftedGaussLegendre(n, 0, 1);
    points.reserve(n);
    for (const auto& point : line.first) {
      points.push_back(std::array<double, 1>{point});
    }
    weights = std::move(line.second);
  } else {
    std::size_t count = 1;
    for (std::size_t i = 0; i < Dim; ++i) {
      count *= n;
    }
    std::vector<double> flat(Dim * count);
    weights.resize(count);
    if constexpr (Dim == 2) {
      quadrature::TriangleQuadrature(
          reinterpret_cast<double (*)[2]>(flat.data()), weights.data(), n);
    } else {
      quadrature::TetrahedronQuadrature(
          reinterpret_cast<double (*)[3]>(flat.data()), weights.data(), n);
    }
    points.resize(count);
    for (std::size_t q = 0; q < count; ++q) {
      for (std::size_t i = 0; i < Dim; ++i) {
        points[q][i] = flat[q * Dim + i];
      }
    }
  }
  return {points, weights};
}

constexpr double Pi = 3.14159265358979323846;

// Hesthaven/Warburton warp&blend blending parameters, indexed by polynomial degree - 1.
constexpr std::array<double, 15> WarpBlendAlpha2D = {0.0000,
                                                     0.0000,
                                                     1.4152,
                                                     0.1001,
                                                     0.2751,
                                                     0.9800,
                                                     1.0999,
                                                     1.2832,
                                                     1.3648,
                                                     1.4773,
                                                     1.4959,
                                                     1.5743,
                                                     1.5770,
                                                     1.6223,
                                                     1.6258};

// The same, for the tetrahedron.
constexpr std::array<double, 15> WarpBlendAlpha3D = {0.0000,
                                                     0.0000,
                                                     0.0000,
                                                     0.1002,
                                                     1.1332,
                                                     1.5608,
                                                     1.3413,
                                                     1.2577,
                                                     1.1603,
                                                     1.10153,
                                                     0.6080,
                                                     0.4523,
                                                     0.8856,
                                                     0.8717,
                                                     0.9655};

//! Gauss-Lobatto points on [-1, 1]; returns degree + 1 points.
std::vector<double> gaussLobatto(std::size_t degree) {
  std::vector<double> points;
  points.reserve(degree + 1);
  points.push_back(-1.0);
  if (degree >= 2) {
    std::vector<double> inner(degree - 1);
    std::vector<double> innerWeights(degree - 1);
    quadrature::GaussJacobi(inner.data(), innerWeights.data(), degree - 1, 1, 1);
    std::sort(inner.begin(), inner.end());
    points.insert(points.end(), inner.begin(), inner.end());
  }
  points.push_back(1.0);
  return points;
}

/**
 * Interpolates the displacement between the equidistant and the Gauss-Lobatto nodes of the given
 * degree onto @p positions , divided by @f$ 1 - x^2 @f$ .
 *
 * The division is carried out analytically: the two endpoint factors of the Lagrange polynomial
 * are exactly @f$ (x \mp 1) @f$ and are simply left out, which keeps the result finite (and
 * accurate) at and near the interval ends. The endpoint terms themselves vanish because the
 * Gauss-Lobatto and the equidistant nodes coincide there.
 */
std::vector<double> warpInterpolation(std::size_t degree, const std::vector<double>& positions) {
  const auto lobatto = gaussLobatto(degree);
  std::vector<double> equidistant(degree + 1);
  for (std::size_t i = 0; i <= degree; ++i) {
    equidistant[i] = -1.0 + (2.0 * static_cast<double>(i)) / static_cast<double>(degree);
  }

  std::vector<double> warp(positions.size(), 0.0);
  for (std::size_t k = 0; k < positions.size(); ++k) {
    for (std::size_t i = 0; i <= degree; ++i) {
      double term = lobatto[i] - equidistant[i];
      for (std::size_t j = 1; j < degree; ++j) {
        if (i != j) {
          term *= (positions[k] - equidistant[j]) / (equidistant[i] - equidistant[j]);
        }
      }
      if (i != 0) {
        term /= -(equidistant[i] - equidistant[0]);
      }
      if (i != degree) {
        term /= equidistant[i] - equidistant[degree];
      }
      warp[k] += term;
    }
  }
  return warp;
}

/**
 * The same interpolant on the reversed node set, which is what Hesthaven/Warburton use inside the
 * tetrahedron. Both node sets are symmetric about zero, so reversing them just mirrors the
 * interpolant.
 */
std::vector<double> warpInterpolationReversed(std::size_t degree, std::vector<double> positions) {
  for (auto& position : positions) {
    position = -position;
  }
  auto warp = warpInterpolation(degree, positions);
  for (auto& value : warp) {
    value = -value;
  }
  return warp;
}

} // namespace

namespace seissol::numerical {

DenseMatrix DenseMatrix::multiply(const DenseMatrix& other) const {
  assert(cols_ == other.rows());
  DenseMatrix result(rows_, other.cols());
  for (std::size_t i = 0; i < rows_; ++i) {
    for (std::size_t k = 0; k < cols_; ++k) {
      const auto factor = (*this)(i, k);
      for (std::size_t j = 0; j < other.cols(); ++j) {
        result(i, j) += factor * other(k, j);
      }
    }
  }
  return result;
}

DenseMatrix DenseMatrix::inverse() const {
  assert(rows_ == cols_);
  const auto n = rows_;
  DenseMatrix work = *this;
  DenseMatrix result(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    result(i, i) = 1.0;
  }

  for (std::size_t column = 0; column < n; ++column) {
    // partial pivoting
    std::size_t pivot = column;
    for (std::size_t row = column + 1; row < n; ++row) {
      if (std::abs(work(row, column)) > std::abs(work(pivot, column))) {
        pivot = row;
      }
    }
    if (work(pivot, column) == 0.0) {
      logError() << "Tried to invert a singular matrix of size" << n << "in the projection setup.";
    }
    if (pivot != column) {
      for (std::size_t j = 0; j < n; ++j) {
        std::swap(work(pivot, j), work(column, j));
        std::swap(result(pivot, j), result(column, j));
      }
    }

    const auto scale = 1.0 / work(column, column);
    for (std::size_t j = 0; j < n; ++j) {
      work(column, j) *= scale;
      result(column, j) *= scale;
    }
    for (std::size_t row = 0; row < n; ++row) {
      if (row != column) {
        const auto factor = work(row, column);
        if (factor != 0.0) {
          for (std::size_t j = 0; j < n; ++j) {
            work(row, j) -= factor * work(column, j);
            result(row, j) -= factor * result(column, j);
          }
        }
      }
    }
  }
  return result;
}

namespace projection {

template <std::size_t Dim>
std::vector<std::array<unsigned, Dim>> modalIndices(std::size_t order) {
  std::vector<std::array<unsigned, Dim>> indices;
  indices.reserve(modalSize(Dim, order));

  std::array<unsigned, Dim> current{};
  const auto enumerate = [&](auto&& self, std::size_t dimension, std::size_t remaining) -> void {
    if (dimension == 0) {
      current[0] = static_cast<unsigned>(remaining);
      indices.push_back(current);
      return;
    }
    for (std::size_t value = 0; value <= remaining; ++value) {
      current[dimension] = static_cast<unsigned>(value);
      self(self, dimension - 1, remaining - value);
    }
  };

  for (std::size_t degree = 0; degree < order; ++degree) {
    enumerate(enumerate, Dim - 1, degree);
  }

  assert(indices.size() == modalSize(Dim, order));
  return indices;
}

template std::vector<std::array<unsigned, 1>> modalIndices<1>(std::size_t);
template std::vector<std::array<unsigned, 2>> modalIndices<2>(std::size_t);
template std::vector<std::array<unsigned, 3>> modalIndices<3>(std::size_t);

std::vector<std::array<double, 2>> nodalPoints2D(std::size_t order) {
  assert(order >= 1);
  if (order == 1) {
    // the Lagrange space is spanned by the constant; any interior point does
    return {{1.0 / 3.0, 1.0 / 3.0}};
  }

  const auto degree = order - 1;
  const auto alpha = degree <= WarpBlendAlpha2D.size() ? WarpBlendAlpha2D[degree - 1] : 5.0 / 3.0;
  const auto scale = 1.0 / static_cast<double>(degree);

  // equidistant lattice in barycentric coordinates; eta (= lambda2) slowest
  std::vector<double> lambda1;
  std::vector<double> lambda2;
  lambda1.reserve(modalSize(2, order));
  lambda2.reserve(modalSize(2, order));
  for (std::size_t j = 0; j <= degree; ++j) {
    for (std::size_t i = 0; i + j <= degree; ++i) {
      lambda1.push_back(static_cast<double>(i) * scale);
      lambda2.push_back(static_cast<double>(j) * scale);
    }
  }
  const auto count = lambda1.size();
  std::vector<double> lambda3(count);
  for (std::size_t p = 0; p < count; ++p) {
    lambda3[p] = 1.0 - lambda1[p] - lambda2[p];
  }

  // equilateral reference triangle with the vertices (-1,-1/sqrt(3)), (1,-1/sqrt(3)), (0,2/sqrt(3))
  std::vector<double> x(count);
  std::vector<double> y(count);
  for (std::size_t p = 0; p < count; ++p) {
    x[p] = lambda3[p] - lambda2[p];
    y[p] = (2.0 * lambda1[p] - lambda2[p] - lambda3[p]) / std::sqrt(3.0);
  }

  const auto difference = [&](const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(count);
    for (std::size_t p = 0; p < count; ++p) {
      result[p] = a[p] - b[p];
    }
    return result;
  };

  const std::array<std::vector<double>, 3> warp = {
      warpInterpolation(degree, difference(lambda3, lambda2)),
      warpInterpolation(degree, difference(lambda1, lambda3)),
      warpInterpolation(degree, difference(lambda2, lambda1))};
  const std::array<const std::vector<double>*, 3> blendOuter = {&lambda2, &lambda1, &lambda1};
  const std::array<const std::vector<double>*, 3> blendInner = {&lambda3, &lambda3, &lambda2};
  const std::array<const std::vector<double>*, 3> attenuate = {&lambda1, &lambda2, &lambda3};

  for (std::size_t direction = 0; direction < 3; ++direction) {
    const auto angle = 2.0 * Pi * static_cast<double>(direction) / 3.0;
    for (std::size_t p = 0; p < count; ++p) {
      const auto blend = 4.0 * (*blendOuter[direction])[p] * (*blendInner[direction])[p];
      const auto damping = alpha * (*attenuate[direction])[p];
      const auto value = warp[direction][p] * blend * (1.0 + damping * damping);
      x[p] += std::cos(angle) * value;
      y[p] += std::sin(angle) * value;
    }
  }

  std::vector<std::array<double, 2>> points(count);
  for (std::size_t p = 0; p < count; ++p) {
    points[p][0] = (std::sqrt(3.0) * y[p] + 1.0) / 3.0;
    points[p][1] = (2.0 - 3.0 * x[p] - std::sqrt(3.0) * y[p]) / 6.0;
  }
  return points;
}

std::vector<std::array<double, 3>> nodalPoints3D(std::size_t order, NodalSet set) {
  if (set == NodalSet::Stroud) {
    auto quadrature = simplexQuadrature<3>(order + 1);
    assert(quadrature.first.size() == nodalSize(3, order, set));
    return std::move(quadrature.first);
  }

  assert(order >= 1);
  if (order == 1) {
    return {{0.25, 0.25, 0.25}};
  }

  const auto degree = order - 1;
  const auto alpha = degree <= WarpBlendAlpha3D.size() ? WarpBlendAlpha3D[degree - 1] : 1.0;
  const auto scale = 1.0 / static_cast<double>(degree);
  constexpr double Tolerance = 1e-10;

  // equidistant lattice in barycentric coordinates; zeta slowest, xi fastest. lambda[0..2] are
  // the coordinates of the reference tetrahedron, lambda[3] the one belonging to its origin.
  std::array<std::vector<double>, 4> lambda;
  for (std::size_t k = 0; k <= degree; ++k) {
    for (std::size_t j = 0; j + k <= degree; ++j) {
      for (std::size_t i = 0; i + j + k <= degree; ++i) {
        lambda[0].push_back(static_cast<double>(i) * scale);
        lambda[1].push_back(static_cast<double>(j) * scale);
        lambda[2].push_back(static_cast<double>(k) * scale);
        lambda[3].push_back(1.0 - static_cast<double>(i + j + k) * scale);
      }
    }
  }
  const auto count = lambda[0].size();
  assert(count == nodalSize(3, order, set));

  // the equilateral reference tetrahedron the warp is applied in; vertex v[i] carries lambda[i]
  const std::array<std::array<double, 3>, 4> vertices = {
      std::array<double, 3>{1.0, -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(6.0)},
      std::array<double, 3>{0.0, 2.0 / std::sqrt(3.0), -1.0 / std::sqrt(6.0)},
      std::array<double, 3>{0.0, 0.0, 3.0 / std::sqrt(6.0)},
      std::array<double, 3>{-1.0, -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(6.0)}};

  // the two orthogonal in-plane directions of each face, indexed as the faces below
  const auto subtract = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  };
  const auto middle = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::array<double, 3>{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
  };
  const auto normalize = [](std::array<double, 3> a) {
    const auto norm = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    for (auto& value : a) {
      value /= norm;
    }
    return a;
  };
  const std::array<std::array<double, 3>, 4> tangent1 = {
      normalize(subtract(vertices[0], vertices[3])),
      normalize(subtract(vertices[0], vertices[3])),
      normalize(subtract(vertices[1], vertices[0])),
      normalize(subtract(vertices[1], vertices[3]))};
  const std::array<std::array<double, 3>, 4> tangent2 = {
      normalize(subtract(vertices[1], middle(vertices[3], vertices[0]))),
      normalize(subtract(vertices[2], middle(vertices[3], vertices[0]))),
      normalize(subtract(vertices[2], middle(vertices[0], vertices[1]))),
      normalize(subtract(vertices[2], middle(vertices[3], vertices[1])))};

  // the face opposite to lambda[opposite[f]], warped in the plane spanned by the other three
  const std::array<std::array<std::size_t, 4>, 4> faceLambda = {
      std::array<std::size_t, 4>{2, 1, 3, 0},
      std::array<std::size_t, 4>{1, 2, 3, 0},
      std::array<std::size_t, 4>{3, 2, 0, 1},
      std::array<std::size_t, 4>{0, 2, 3, 1}};

  std::vector<std::array<double, 3>> position(count, std::array<double, 3>{});
  for (std::size_t p = 0; p < count; ++p) {
    for (std::size_t l = 0; l < 4; ++l) {
      for (std::size_t d = 0; d < 3; ++d) {
        position[p][d] += lambda[l][p] * vertices[l][d];
      }
    }
  }

  const auto difference = [&](const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(count);
    for (std::size_t p = 0; p < count; ++p) {
      result[p] = a[p] - b[p];
    }
    return result;
  };

  std::vector<std::array<double, 3>> shift(count, std::array<double, 3>{});
  for (std::size_t face = 0; face < 4; ++face) {
    const auto& la = lambda[faceLambda[face][0]];
    const auto& lb = lambda[faceLambda[face][1]];
    const auto& lc = lambda[faceLambda[face][2]];
    const auto& ld = lambda[faceLambda[face][3]];

    // the 2D warp of the face, in the frame spanned by tangent1/tangent2
    const std::array<std::vector<double>, 3> warp = {
        warpInterpolationReversed(degree, difference(ld, lc)),
        warpInterpolationReversed(degree, difference(lb, ld)),
        warpInterpolationReversed(degree, difference(lc, lb))};
    const std::array<const std::vector<double>*, 3> blendOuter = {&lc, &lb, &lb};
    const std::array<const std::vector<double>*, 3> blendInner = {&ld, &ld, &lc};
    const std::array<const std::vector<double>*, 3> attenuate = {&lb, &lc, &ld};

    std::vector<double> shift1(count, 0.0);
    std::vector<double> shift2(count, 0.0);
    for (std::size_t direction = 0; direction < 3; ++direction) {
      const auto angle = 2.0 * Pi * static_cast<double>(direction) / 3.0;
      for (std::size_t p = 0; p < count; ++p) {
        const auto blend = 4.0 * (*blendOuter[direction])[p] * (*blendInner[direction])[p];
        const auto damping = alpha * (*attenuate[direction])[p];
        const auto value = warp[direction][p] * blend * (1.0 + damping * damping);
        shift1[p] += std::cos(angle) * value;
        shift2[p] += std::sin(angle) * value;
      }
    }

    for (std::size_t p = 0; p < count; ++p) {
      const auto contribution = [&](std::size_t d) {
        return shift1[p] * tangent1[face][d] + shift2[p] * tangent2[face][d];
      };

      // On the boundary of this face the warp is prescribed exactly and replaces whatever the
      // previously visited faces contributed; in the interior it is blended into the volume.
      const bool onBoundary = la[p] < Tolerance && (static_cast<int>(lb[p] > Tolerance) +
                                                    static_cast<int>(lc[p] > Tolerance) +
                                                    static_cast<int>(ld[p] > Tolerance)) < 3;
      if (onBoundary) {
        for (std::size_t d = 0; d < 3; ++d) {
          shift[p][d] = contribution(d);
        }
      } else {
        const auto denominator =
            (lb[p] + 0.5 * la[p]) * (lc[p] + 0.5 * la[p]) * (ld[p] + 0.5 * la[p]);
        auto blend = lb[p] * lc[p] * ld[p];
        if (denominator > Tolerance) {
          blend *= (1.0 + (alpha * la[p]) * (alpha * la[p])) / denominator;
        }
        for (std::size_t d = 0; d < 3; ++d) {
          shift[p][d] += blend * contribution(d);
        }
      }
    }
  }

  for (std::size_t p = 0; p < count; ++p) {
    for (std::size_t d = 0; d < 3; ++d) {
      position[p][d] += shift[p][d];
    }
  }

  // back to the reference tetrahedron: undo the barycentric map given by `vertices`
  DenseMatrix basis(3, 3);
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t l = 0; l < 3; ++l) {
      basis(d, l) = vertices[l][d] - vertices[3][d];
    }
  }
  const auto inverse = basis.inverse();

  std::vector<std::array<double, 3>> points(count);
  for (std::size_t p = 0; p < count; ++p) {
    for (std::size_t d = 0; d < 3; ++d) {
      double value = 0.0;
      for (std::size_t e = 0; e < 3; ++e) {
        value += inverse(d, e) * (position[p][e] - vertices[3][e]);
      }
      points[p][d] = value;
    }
  }
  return points;
}

std::vector<double> stroudWeights3D(std::size_t order) {
  return std::move(simplexQuadrature<3>(order + 1).second);
}

template <std::size_t Dim>
DenseMatrix nodalToModal(std::size_t order, NodalSet set) {
  const auto indices = modalIndices<Dim>(order);
  if constexpr (Dim == 3) {
    if (set == NodalSet::Stroud) {
      // over-determined nodal set: the transform is the (exact) L2 projection using the
      // conical-product quadrature the nodes stem from
      const auto points = nodalPoints3D(order, set);
      const auto weights = stroudWeights3D(order);
      DenseMatrix vandermonde(points.size(), indices.size());
      for (std::size_t n = 0; n < points.size(); ++n) {
        for (std::size_t b = 0; b < indices.size(); ++b) {
          vandermonde(n, b) = functions::DubinerP<3>(indices[b], points[n]);
        }
      }
      std::vector<double> mass(indices.size(), 0.0);
      for (std::size_t b = 0; b < indices.size(); ++b) {
        for (std::size_t n = 0; n < points.size(); ++n) {
          mass[b] += weights[n] * vandermonde(n, b) * vandermonde(n, b);
        }
      }
      DenseMatrix result(indices.size(), points.size());
      for (std::size_t b = 0; b < indices.size(); ++b) {
        for (std::size_t n = 0; n < points.size(); ++n) {
          result(b, n) = weights[n] * vandermonde(n, b) / mass[b];
        }
      }
      return result;
    }
  } else {
    static_assert(Dim == 2, "Only 2D and 3D nodal bases are supported.");
    assert(set == NodalSet::WarpBlend);
  }

  // unisolvent nodal set: the transform is the inverse Vandermonde matrix
  const auto points = [&]() {
    if constexpr (Dim == 3) {
      return nodalPoints3D(order, set);
    } else {
      return nodalPoints2D(order);
    }
  }();
  assert(points.size() == indices.size());
  DenseMatrix vandermonde(points.size(), indices.size());
  for (std::size_t n = 0; n < points.size(); ++n) {
    for (std::size_t b = 0; b < indices.size(); ++b) {
      vandermonde(n, b) = functions::DubinerP<Dim>(indices[b], points[n]);
    }
  }
  return vandermonde.inverse();
}

template DenseMatrix nodalToModal<2>(std::size_t, NodalSet);
template DenseMatrix nodalToModal<3>(std::size_t, NodalSet);

template <std::size_t From, std::size_t To>
DenseMatrix build(const std::vector<std::array<double, From>>& referenceTargetPoints,
                  std::size_t targetDegree,
                  const AffineMap<From, To>& map,
                  const Spec& spec) {
  const auto sourceIndices = modalIndices<To>(spec.order);
  const auto sourceCount = sourceIndices.size();
  const auto pointCount = referenceTargetPoints.size();

  const auto evaluate = [&](const std::array<double, To>& position, std::size_t basis) {
    if (spec.derivative.has_value()) {
      assert(spec.derivative.value() < To);
      return functions::gradDubinerP<To>(sourceIndices[basis], position)[spec.derivative.value()];
    }
    return functions::DubinerP<To>(sourceIndices[basis], position);
  };

  DenseMatrix matrix(pointCount, sourceCount);

  if (spec.target == Target::Interpolate) {
    for (std::size_t p = 0; p < pointCount; ++p) {
      const auto position = map(referenceTargetPoints[p]);
      for (std::size_t b = 0; b < sourceCount; ++b) {
        matrix(p, b) = evaluate(position, b);
      }
    }
  } else {
    // L2 projection onto the Lagrange space of the target points. With psi_i the modal basis of
    // the target space and V(p,i) = psi_i(x_p), the Lagrange mass matrix is
    //   M = |det A| V^-T diag(m_i) V^-1  ,  m_i = int psi_i^2 ,
    // and the right hand side is |det A| V^-T C with C(i,b) = int psi_i * phi_b(A .) ; hence
    //   P = M^-1 B = V diag(1/m_i) C ,
    // in particular the subcell volume cancels.
    const auto targetIndices = modalIndices<From>(targetDegree + 1);
    const auto targetCount = targetIndices.size();
    if (targetCount != pointCount) {
      logError() << "The target point set of size" << pointCount
                 << "is not unisolvent for a polynomial space of degree" << targetDegree
                 << "and dimension" << From << "(expected" << targetCount << "points).";
    }

    // exact for degree targetDegree + (order - 1)
    const auto quadratureSize = (targetDegree + spec.order) / 2 + 1;
    const auto quadrature = simplexQuadrature<From>(quadratureSize);
    const auto& quadraturePoints = quadrature.first;
    const auto& quadratureWeights = quadrature.second;

    DenseMatrix coupling(targetCount, sourceCount);
    std::vector<double> mass(targetCount, 0.0);
    std::vector<double> targetValues(targetCount);
    for (std::size_t q = 0; q < quadraturePoints.size(); ++q) {
      const auto position = map(quadraturePoints[q]);
      for (std::size_t i = 0; i < targetCount; ++i) {
        targetValues[i] = functions::DubinerP<From>(targetIndices[i], quadraturePoints[q]);
        mass[i] += quadratureWeights[q] * targetValues[i] * targetValues[i];
      }
      for (std::size_t b = 0; b < sourceCount; ++b) {
        const auto sourceValue = evaluate(position, b);
        for (std::size_t i = 0; i < targetCount; ++i) {
          coupling(i, b) += quadratureWeights[q] * targetValues[i] * sourceValue;
        }
      }
    }
    for (std::size_t i = 0; i < targetCount; ++i) {
      for (std::size_t b = 0; b < sourceCount; ++b) {
        coupling(i, b) /= mass[i];
      }
    }

    DenseMatrix vandermonde(pointCount, targetCount);
    for (std::size_t p = 0; p < pointCount; ++p) {
      for (std::size_t i = 0; i < targetCount; ++i) {
        vandermonde(p, i) = functions::DubinerP<From>(targetIndices[i], referenceTargetPoints[p]);
      }
    }

    matrix = vandermonde.multiply(coupling);
  }

  if (spec.source == Source::Nodal) {
    matrix = matrix.multiply(nodalToModal<To>(spec.order, spec.nodalSet));
    assert(matrix.cols() == nodalSize(To, spec.order, spec.nodalSet));
  }

  return matrix;
}

template DenseMatrix build<2, 2>(const std::vector<std::array<double, 2>>&,
                                 std::size_t,
                                 const AffineMap<2, 2>&,
                                 const Spec&);
template DenseMatrix build<2, 3>(const std::vector<std::array<double, 2>>&,
                                 std::size_t,
                                 const AffineMap<2, 3>&,
                                 const Spec&);
template DenseMatrix build<3, 3>(const std::vector<std::array<double, 3>>&,
                                 std::size_t,
                                 const AffineMap<3, 3>&,
                                 const Spec&);

template <std::size_t From, std::size_t To>
Table<From, To>::Table(const std::vector<AffineMap<From, To>>& subcells,
                       const std::vector<std::array<double, From>>& referenceTargetPoints,
                       std::size_t targetDegree,
                       std::size_t leadingDimension,
                       const Spec& spec,
                       std::size_t minOrder,
                       std::size_t maxOrder)
    : subcellCount_(subcells.size()), minOrder_(minOrder), maxOrder_(maxOrder),
      leadingDimension_(leadingDimension) {
  assert(minOrder >= 1 && minOrder <= maxOrder);
  assert(leadingDimension >= referenceTargetPoints.size());

  orderOffsets_.resize(maxOrder - minOrder + 1);
  std::size_t offset = 0;
  for (std::size_t order = minOrder; order <= maxOrder; ++order) {
    orderOffsets_[order - minOrder] = offset;
    const auto columns =
        spec.source == Source::Nodal ? nodalSize(To, order, spec.nodalSet) : modalSize(To, order);
    offset += columns * leadingDimension;
  }
  subcellStride_ = offset;

  if (subcellCount_ * subcellStride_ == 0) {
    return;
  }

  storage_.resize(subcellCount_ * subcellStride_);
  std::fill(storage_.begin(), storage_.end(), static_cast<real>(0));

  for (std::size_t subcell = 0; subcell < subcellCount_; ++subcell) {
    for (std::size_t order = minOrder; order <= maxOrder; ++order) {
      auto localSpec = spec;
      localSpec.order = order;
      const auto matrix =
          build<From, To>(referenceTargetPoints, targetDegree, subcells[subcell], localSpec);
      auto* target = storage_.data() + this->offset(subcell, order);
      for (std::size_t p = 0; p < matrix.rows(); ++p) {
        for (std::size_t b = 0; b < matrix.cols(); ++b) {
          target[p + b * leadingDimension] = static_cast<real>(matrix(p, b));
        }
      }
    }
  }
}

template class Table<2, 2>;
template class Table<2, 3>;
template class Table<3, 3>;

} // namespace projection

} // namespace seissol::numerical
