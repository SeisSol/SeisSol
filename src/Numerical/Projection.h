// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_PROJECTION_H_
#define SEISSOL_SRC_NUMERICAL_PROJECTION_H_

#include "Kernels/Precision.h"
#include "Memory/MemoryAllocator.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::numerical {

/**
 * @brief An affine map from the @c From -dimensional reference simplex into @c To -dimensional
 * space, i.e. @f$ y = o + M x @f$ .
 *
 * Used both for cell subdivision (From == To) and for embedding a face into the volume
 * (From == 2, To == 3).
 */
template <std::size_t From, std::size_t To = From>
struct AffineMap {
  std::array<double, To> offset{};
  //! matrix[i][j] = dy_i / dx_j
  std::array<std::array<double, From>, To> matrix{};

  static AffineMap identity() {
    static_assert(From == To, "The identity map requires From == To.");
    AffineMap map;
    for (std::size_t i = 0; i < To; ++i) {
      map.matrix[i][i] = 1;
    }
    return map;
  }

  /**
   * @brief Reconstructs the map from the images of the reference simplex vertices,
   * i.e. from @f$ f(0), f(e_1), \dots, f(e_{From}) @f$ .
   */
  static AffineMap fromVertices(const std::vector<std::array<double, To>>& vertices) {
    assert(vertices.size() == From + 1);
    AffineMap map;
    map.offset = vertices[0];
    for (std::size_t i = 0; i < To; ++i) {
      for (std::size_t j = 0; j < From; ++j) {
        map.matrix[i][j] = vertices[j + 1][i] - vertices[0][i];
      }
    }
    return map;
  }

  [[nodiscard]] std::array<double, To> operator()(const std::array<double, From>& point) const {
    std::array<double, To> result = offset;
    for (std::size_t i = 0; i < To; ++i) {
      for (std::size_t j = 0; j < From; ++j) {
        result[i] += matrix[i][j] * point[j];
      }
    }
    return result;
  }

  //! @brief Returns this ∘ inner.
  template <std::size_t Inner>
  [[nodiscard]] AffineMap<Inner, To> compose(const AffineMap<Inner, From>& inner) const {
    AffineMap<Inner, To> result;
    result.offset = (*this)(inner.offset);
    for (std::size_t i = 0; i < To; ++i) {
      for (std::size_t j = 0; j < Inner; ++j) {
        double sum = 0;
        for (std::size_t k = 0; k < From; ++k) {
          sum += matrix[i][k] * inner.matrix[k][j];
        }
        result.matrix[i][j] = sum;
      }
    }
    return result;
  }

  [[nodiscard]] std::vector<std::array<double, To>>
      apply(const std::vector<std::array<double, From>>& points) const {
    std::vector<std::array<double, To>> result;
    result.reserve(points.size());
    for (const auto& point : points) {
      result.emplace_back((*this)(point));
    }
    return result;
  }
};

/**
 * @brief A minimal dense row-major matrix, for setup code only.
 */
class DenseMatrix {
  public:
  DenseMatrix() = default;
  DenseMatrix(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), data_(rows * cols) {}

  double& operator()(std::size_t row, std::size_t col) { return data_[row * cols_ + col]; }
  [[nodiscard]] double operator()(std::size_t row, std::size_t col) const {
    return data_[row * cols_ + col];
  }

  [[nodiscard]] std::size_t rows() const { return rows_; }
  [[nodiscard]] std::size_t cols() const { return cols_; }
  [[nodiscard]] const std::vector<double>& data() const { return data_; }

  //! @brief Returns this * other.
  [[nodiscard]] DenseMatrix multiply(const DenseMatrix& other) const;

  //! @brief Returns the inverse; the matrix needs to be square and regular.
  [[nodiscard]] DenseMatrix inverse() const;

  private:
  std::size_t rows_{};
  std::size_t cols_{};
  std::vector<double> data_;
};

namespace projection {

//! @brief How the source data is represented.
enum class Source {
  //! Dubiner (modal) coefficients.
  Modal,
  //! Point values at the canonical nodal points of the respective dimension.
  Nodal
};

//! @brief How the source function is transferred onto the target points.
enum class Target {
  //! Evaluate the source function at the target points.
  Interpolate,
  //! L2-project the source function onto the Lagrange space spanned by the target points.
  Project
};

/**
 * @brief Which nodal point set a Source::Nodal representation lives on.
 *
 * Both are in use in the volume: the generated plasticity matrices (@c vNodes , @c v , @c vInv )
 * are built for exactly one of them, selected at configure time by @c PLASTICITY_METHOD . The 2D
 * face nodes (@c nodes2D ) are always WarpBlend.
 */
enum class NodalSet {
  //! Hesthaven/Warburton warp&blend points: unisolvent, one point per basis function
  //! (@c PLASTICITY_METHOD=nb , the default).
  WarpBlend,
  //! Conical-product (Stroud) quadrature points: over-determined, the nodal-to-modal transform is
  //! the corresponding L2 projection (@c PLASTICITY_METHOD=ip ). 3D only.
  Stroud
};

//! @brief Number of Dubiner basis functions of convergence order @p order .
constexpr std::size_t modalSize(std::size_t dim, std::size_t order) {
  std::size_t count = order;
  for (std::size_t i = 1; i < dim; ++i) {
    count *= order + i;
    count /= i + 1;
  }
  return count;
}

//! @brief Number of points of the nodal point set @p set of convergence order @p order .
constexpr std::size_t nodalSize(std::size_t dim, std::size_t order, NodalSet set) {
  if (set == NodalSet::Stroud) {
    // cf. quadrature::TetrahedronQuadrature(order + 1); defined for dim == 3 only
    return (order + 1) * (order + 1) * (order + 1);
  }
  return modalSize(dim, order);
}

/**
 * @brief Enumerates the Dubiner multi-indices of convergence order @p order .
 *
 * The enumeration is by total degree; within one degree the last coordinate varies slowest and
 * the first coordinate takes up the remainder. This matches
 * basisFunction::SampledBasisFunctions (Dim == 3) and basisFunction::tri_dubiner (Dim == 2), and
 * therefore also the modal index order of the generated V2mTo2n / MV2nTo2m / v / vInv matrices.
 */
template <std::size_t Dim>
std::vector<std::array<unsigned, Dim>> modalIndices(std::size_t order);

/**
 * @brief The nodal points of the 2D face basis (the generated @c nodes2D ).
 *
 * Hesthaven/Warburton warp&blend nodes; the enumeration runs over the equidistant lattice with
 * eta slowest.
 */
std::vector<std::array<double, 2>> nodalPoints2D(std::size_t order);

/**
 * @brief The nodal points of the 3D volume basis (the generated @c vNodes ).
 *
 * NodalSet::Stroud gives the conical-product points, identical to
 * quadrature::TetrahedronQuadrature(order + 1); NodalSet::WarpBlend gives the
 * Hesthaven/Warburton warp&blend points, enumerated over the equidistant lattice with zeta
 * slowest and xi fastest.
 */
std::vector<std::array<double, 3>> nodalPoints3D(std::size_t order, NodalSet set);

//! @brief The quadrature weights belonging to nodalPoints3D(order, NodalSet::Stroud).
std::vector<double> stroudWeights3D(std::size_t order);

/**
 * @brief The nodal-to-modal transform of a nodal set, i.e. the generated @c MV2nTo2m (Dim == 2)
 * resp. @c vInv (Dim == 3). Shape: [modalSize][nodalSize].
 *
 * For a unisolvent set this is the inverse Vandermonde matrix; for the over-determined Stroud set
 * it is the L2 projection using the quadrature the points stem from.
 */
template <std::size_t Dim>
DenseMatrix nodalToModal(std::size_t order, NodalSet set);

struct Spec {
  //! Convergence order of the source representation.
  std::size_t order{1};
  Source source{Source::Modal};
  Target target{Target::Interpolate};
  //! The point set a Source::Nodal representation lives on; ignored for Source::Modal.
  NodalSet nodalSet{NodalSet::WarpBlend};
  //! If set, the derivative direction (in the coordinates of the source reference cell).
  std::optional<std::size_t> derivative;
};

/**
 * @brief Builds a single projection matrix, laid out row-major as [targetPoint][sourceCoefficient].
 *
 * @param referenceTargetPoints The target points on the @c From -dimensional reference simplex,
 *        i.e. the output of io::instance::geometry::points* for @p targetDegree .
 * @param targetDegree The polynomial degree of the target (Lagrange) space; only relevant for
 *        Target::Project.
 * @param map Maps the reference simplex of the target onto the source reference cell. For volume
 *        output this is the subcell map; for face output the subcell map composed with the
 *        face embedding.
 *
 * The number of columns is projection::modalSize(To, order) for Source::Modal and
 * projection::nodalSize(To, order, spec.nodalSet) for Source::Nodal, i.e. the nodal-to-modal
 * transform is folded into the matrix.
 */
template <std::size_t From, std::size_t To>
DenseMatrix build(const std::vector<std::array<double, From>>& referenceTargetPoints,
                  std::size_t targetDegree,
                  const AffineMap<From, To>& map,
                  const Spec& spec);

/**
 * @brief A set of projection matrices for all convergence orders and all subcells, stored in the
 * memory layout expected by the generated (Yateto) kernels.
 *
 * The generated tensors are declared as [points][coefficients] with an aligned stride on the point
 * dimension; i.e. the address of entry (p, c) is p + c * leadingDimension .
 */
template <std::size_t From, std::size_t To>
class Table {
  public:
  Table() = default;

  /**
   * @param subcells One map per subcell (cf. io::instance::geometry::subdivideMaps).
   * @param leadingDimension The padded point count of the target tensor; obtain it from the
   *        generated metadata as tensor::X::size(...) / tensor::X::Shape[...][1] .
   * @param minOrder, maxOrder Inclusive range of convergence orders to generate.
   */
  Table(const std::vector<AffineMap<From, To>>& subcells,
        const std::vector<std::array<double, From>>& referenceTargetPoints,
        std::size_t targetDegree,
        std::size_t leadingDimension,
        const Spec& spec,
        std::size_t minOrder,
        std::size_t maxOrder);

  //! @brief The matrix for the given subcell and convergence order.
  [[nodiscard]] const real* operator()(std::size_t subcell, std::size_t order) const {
    return storage_.data() + offset(subcell, order);
  }

  [[nodiscard]] std::size_t subcellCount() const { return subcellCount_; }
  [[nodiscard]] std::size_t leadingDimension() const { return leadingDimension_; }
  //! @brief Total number of stored reals; mostly for diagnostics and tests.
  [[nodiscard]] std::size_t size() const { return storage_.size(); }

  private:
  [[nodiscard]] std::size_t offset(std::size_t subcell, std::size_t order) const {
    assert(subcell < subcellCount_);
    assert(order >= minOrder_ && order <= maxOrder_);
    return subcell * subcellStride_ + orderOffsets_[order - minOrder_];
  }

  std::size_t subcellCount_{};
  std::size_t minOrder_{};
  std::size_t maxOrder_{};
  std::size_t leadingDimension_{};
  //! Offset of a given order inside one subcell block.
  std::vector<std::size_t> orderOffsets_;
  std::size_t subcellStride_{};
  memory::MemkindArray<real> storage_{memory::Memkind::Standard};
};

} // namespace projection

} // namespace seissol::numerical

#endif // SEISSOL_SRC_NUMERICAL_PROJECTION_H_
