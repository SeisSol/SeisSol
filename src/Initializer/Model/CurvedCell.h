// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_CURVEDCELL_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_CURVEDCELL_H_

#include "Alignment.h"
#include "Common/Constants.h"
#include "Equations/Setup.h" // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace seissol::initializer {

/**
 * What a cell that may be curved carries of its geometry.
 *
 * The predictor of every cell applies the strong form, whose metric is the inverse Jacobian at a
 * point (CurvedCell::setGradients). A straight-sided cell runs the corrector of the one-per-cell
 * form: tested with the basis, the face terms take the surface Jacobian over the Jacobian
 * determinant at every point of a face its flux is integrated at, together with the rotation into
 * the frame the face has there (CurvedCell::setFace), and all of these are the constants the
 * one-per-cell form carries. A curved cell is weight-adjusted with the weight J W instead, J its
 * Jacobian determinant: its faces take the surface Jacobian alone, its inverse weight takes 1 / J
 * at every point (CurvedCell::setWeightScale), and it forms its volume term itself, in the split
 * form (see kernels::CurvedSplit), from J times the coefficients of the operator without material
 * (CurvedCell::setSplitCoefficients). A weight-adjusted cell that is straight-sided takes its
 * scales the same way, with a J that is constant.
 */
struct CurvedCell {
  /// The reference coordinates of the point `point` the operator is formed at.
  static auto operatorPoint(std::size_t point) -> geometry::CellTransform::VectorEigenT {
    const auto nodes =
        init::operatorNodes::view::create(const_cast<real*>(init::operatorNodes::Values));
    geometry::CellTransform::VectorEigenT result = geometry::CellTransform::VectorEigenT::Zero();
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      if (nodes.isInRange(point, d)) {
        result(d) = nodes(point, d);
      }
    }
    return result;
  }

  /// The reference coordinates, on the reference face, of the point `node` the flux of a face is
  /// integrated at.
  static auto facePoint(std::size_t node) -> geometry::FaceTransform::FaceVectorT {
    const auto points =
        init::faceQuadPoints::view::create(const_cast<real*>(init::faceQuadPoints::Values));
    geometry::FaceTransform::FaceVectorT result = geometry::FaceTransform::FaceVectorT::Zero();
    for (std::size_t d = 0; d < Face::Dim; ++d) {
      if (points.isInRange(node, d)) {
        result(d) = points(node, d);
      }
    }
    return result;
  }

  /// The reference coordinates, in the cell, of the point `node` of the face `side`.
  static auto faceNode(std::size_t side, std::size_t node)
      -> geometry::CellTransform::VectorEigenT {
    return geometry::ReferenceFaceMap(side).faceToCell(facePoint(node));
  }

  /**
   * The rows of the inverse Jacobian at every point the operator is formed at, in the layout the
   * kernels read referenceGradients(dim) in; `target` is LocalIntegrationData::referenceGradients.
   */
  template <typename TargetT>
  static void setGradients(const geometry::CellTransform& transform, TargetT& target) {
    for (std::size_t point = 0; point < OperatorPointCount; ++point) {
      const auto inverse = transform.refToSpaceJacobianInverse(operatorPoint(point));
      setGradientsOfDirection<0>(target[0], point, inverse);
      setGradientsOfDirection<1>(target[1], point, inverse);
      setGradientsOfDirection<2>(target[2], point, inverse);
    }
  }

  /// The Jacobian determinant of a cell at the point `point` the operator is formed at; positive
  /// for a cell that is not turned inside out, whichever way round its vertices are numbered.
  static auto jacobianDeterminant(const geometry::CellTransform& transform, std::size_t point)
      -> double {
    return std::abs(transform.refToSpaceJacobian(operatorPoint(point)).determinant());
  }

  /// 1 / J at every point the operator is formed at, which the inverse weight of a weight-adjusted
  /// cell carries (see setFace); `target` is LocalIntegrationData::weightScale.
  template <typename TargetT>
  static void setWeightScale(const geometry::CellTransform& transform, TargetT& target) {
    for (std::size_t point = 0; point < OperatorPointCount; ++point) {
      target[point] = 1.0 / jacobianDeterminant(transform, point);
    }
  }

  /**
   * J times the coefficients of the operator without material at every point the operator is
   * formed at, from `samples`, the material there; what the volume term of a curved cell reads
   * (see kernels::CurvedSplit). `target` is LocalIntegrationData::splitCoefficients.
   */
  template <typename MaterialT, typename TargetT>
  static void setSplitCoefficients(const geometry::CellTransform& transform,
                                   const MaterialT* samples,
                                   TargetT& target) {
    static_assert(MaterialSampleCount == OperatorPointCount,
                  "the material is sampled where the operator is formed");
    for (std::size_t point = 0; point < OperatorPointCount; ++point) {
      const auto coefficients = seissol::model::getUnweightedCoefficients(samples[point]);
      const double determinant = jacobianDeterminant(transform, point);
      for (std::size_t a = 0; a < coefficients.size(); ++a) {
        target[a][point] = determinant * coefficients[a];
      }
    }
  }

  /**
   * The scalars the normal flux of a curved cell's own trace on a face reads at the point of the
   * face with the material `material` and the scale `scale` (setFace) there, into the slot `node`
   * of `target`, LocalIntegrationData::ownTractionCoefficients of the face: the scale times the
   * coefficients of the operator without material named by OwnTractionCoefficientIndices.
   */
  template <typename MaterialT, typename TargetT>
  static void
      setOwnTraction(const MaterialT& material, double scale, std::size_t node, TargetT& target) {
    const auto coefficients = seissol::model::getUnweightedCoefficients(material);
    for (std::size_t c = 0; c < OwnTractionCoefficientCount; ++c) {
      target[c][node] = scale * coefficients[OwnTractionCoefficientIndices[c]];
    }
  }

  /**
   * The rotation of the face `side` at each of the points its flux is integrated at, into
   * `rotation` in the layout the kernels read TNodes in, and the scale its flux takes there: the
   * surface Jacobian, over the Jacobian determinant of the cell unless the cell is weight-adjusted
   * (its weight carries that, see setWeightScale), with the sign of a flux that is subtracted. The
   * normal is the outward one; the rotation's frame is the one faceAlignedBasis gives at the point,
   * so a straight face gets the rotation it would get as one face.
   */
  static void setFace(const geometry::CellTransform& cell,
                      const geometry::FaceTransform& face,
                      std::size_t side,
                      bool weightAdjusted,
                      real* rotation,
                      std::array<double, FluxFaceNodes>& fluxScale) {
    auto target = init::TNodes::view::create(rotation);
    for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
      const auto point = facePoint(node);
      const double determinant = cell.refToSpaceJacobian(face.refToCell(point)).determinant();
      const auto basis = face.faceAlignedBasis(point);
      // the face normal is the outward one of a cell that is not turned inside out, which the mesh
      // checks; its norm is the surface Jacobian
      fluxScale[node] = -basis[0].norm() / (weightAdjusted ? 1.0 : std::abs(determinant));

      alignas(Alignment) std::array<real, tensor::T::size()> matT{};
      alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
      auto viewT = init::T::view::create(matT.data());
      auto viewTinv = init::Tinv::view::create(matTinv.data());
      seissol::model::getFaceRotationMatrix(
          basis[0].normalized(), basis[1].normalized(), basis[2].normalized(), viewT, viewTinv);
      setRotationOfNode(target, node, viewT);
    }
  }

  /**
   * The rows of a constant inverse Jacobian -- the one of a straight-sided cell -- written into
   * LocalIntegrationData::referenceGradients in whichever layout the build carries them: once, or
   * at every point the operator is formed at.
   */
  template <typename MatrixT, typename TargetT>
  static void setConstantGradients(const MatrixT& inverse, TargetT& target) {
    if constexpr (Curvilinear) {
      for (std::size_t point = 0; point < OperatorPointCount; ++point) {
        setGradientsOfDirection<0>(target[0], point, inverse);
        setGradientsOfDirection<1>(target[1], point, inverse);
        setGradientsOfDirection<2>(target[2], point, inverse);
      }
    } else {
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
        for (std::size_t component = 0; component < Cell::Dim; ++component) {
          target[dim][component] = inverse(dim, component);
        }
      }
    }
  }

  /**
   * The rotation of a straight face, given in the layout of T, written into
   * LocalIntegrationData::faceRotation in whichever layout the build carries it: once, or at every
   * point of the face.
   */
  static void setConstantRotation(const real* rotation, real* target) {
    const auto source = init::T::view::create(const_cast<real*>(rotation));
    if constexpr (Curvilinear) {
      auto view = init::TNodes::view::create(target);
      for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
        setRotationOfNode(view, node, source);
      }
    } else {
      auto view = init::T::view::create(target);
      view.setZero();
      for (std::size_t row = 0; row < tensor::T::Shape[0]; ++row) {
        for (std::size_t column = 0; column < tensor::T::Shape[1]; ++column) {
          if (view.isInRange(row, column)) {
            view(row, column) = source(row, column);
          }
        }
      }
    }
  }

  private:
  template <typename TargetT, typename SourceT>
  static void setRotationOfNode(TargetT& target, std::size_t node, const SourceT& source) {
    for (std::size_t row = 0; row < tensor::T::Shape[0]; ++row) {
      for (std::size_t column = 0; column < tensor::T::Shape[1]; ++column) {
        if (target.isInRange(node, row, column)) {
          target(node, row, column) = source.isInRange(row, column) ? source(row, column) : 0;
        }
      }
    }
  }

  template <std::size_t Dim, typename TargetT, typename MatrixT>
  static void setGradientsOfDirection(TargetT& target, std::size_t point, const MatrixT& inverse) {
    auto view = init::referenceGradients::view<Dim>::create(target);
    for (std::size_t component = 0; component < Cell::Dim; ++component) {
      view(point, component) = inverse(Dim, component);
    }
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_CURVEDCELL_H_
