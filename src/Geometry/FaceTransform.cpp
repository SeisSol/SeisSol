// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FaceTransform.h"

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/MeshReader.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <tuple>
#include <utility>
#include <utils/logger.h>

namespace {

using seissol::Cell;
using seissol::Face;
using seissol::geometry::FaceOrientation;

using OrientationMatrixT = Eigen::Matrix<double, Face::Dim, Face::Dim>;
using OrientationVectorT = Eigen::Vector<double, Face::Dim>;
using EmbeddingT = Eigen::Matrix<double, Cell::Dim, Face::Dim>;
using ProjectionT = Eigen::Matrix<double, Face::Dim, Cell::Dim>;
using CellVectorT = Eigen::Vector<double, Cell::Dim>;

/**
 * Permutes the face coordinates according to the orientation.
 *
 * All three rotations are involutions, hence each matrix is its own inverse; the entries are 0 and
 * +-1 throughout, so applying them is exact.
 */
auto orientationMap(FaceOrientation orientation)
    -> std::pair<OrientationMatrixT, OrientationVectorT> {
  auto matrix = OrientationMatrixT::Identity().eval();
  auto offset = OrientationVectorT::Zero().eval();
  switch (orientation) {
  case FaceOrientation::Rotate0:
    matrix << 0, 1, 1, 0;
    break;
  case FaceOrientation::Rotate1:
    matrix << -1, -1, 0, 1;
    offset << 1, 0;
    break;
  case FaceOrientation::Rotate2:
    matrix << 1, 0, -1, -1;
    offset << 0, 1;
    break;
  case FaceOrientation::Local:
    break;
  }
  return {matrix, offset};
}

/// Embeds the reference face into the given side of the reference cell, plus a left inverse.
auto sideMap(std::size_t side) -> std::tuple<EmbeddingT, CellVectorT, ProjectionT> {
  auto embedding = EmbeddingT::Zero().eval();
  auto offset = CellVectorT::Zero().eval();
  auto projection = ProjectionT::Zero().eval();
  switch (side) {
  case 0:
    embedding << 0, 1, 1, 0, 0, 0;
    projection << 0, 1, 0, 1, 0, 0;
    break;
  case 1:
    embedding << 1, 0, 0, 0, 0, 1;
    projection << 1, 0, 0, 0, 0, 1;
    break;
  case 2:
    embedding << 0, 0, 0, 1, 1, 0;
    projection << 0, 0, 1, 0, 1, 0;
    break;
  case 3:
    embedding << -1, -1, 1, 0, 0, 1;
    offset << 1, 0, 0;
    projection << 0, 1, 0, 0, 0, 1;
    break;
  default:
    logError() << "There is no side" << side << "on a cell; expected 0 <=side <" << Cell::NumFaces
               << ".";
  }
  return {embedding, offset, projection};
}

} // namespace

namespace seissol::geometry {

ReferenceFaceMap::ReferenceFaceMap(std::size_t side, FaceOrientation orientation)
    : side_(side), orientation_(orientation) {
  const auto [orientationMatrix, orientationOffset] = orientationMap(orientation);
  const auto [embedding, embeddingOffset, sideProjection] = sideMap(side);

  // the two affine maps are composed once here, rather than applied one after the other on every
  // evaluation; this also collapses the cancellation each of them would contribute separately
  embedding_ = embedding * orientationMatrix;
  offset_ = embedding * orientationOffset + embeddingOffset;

  projection_ = orientationMatrix * sideProjection;
  projectionOffset_ = -orientationMatrix * (sideProjection * embeddingOffset + orientationOffset);
}

auto ReferenceFaceMap::faceToCell(const FaceVectorT& input) const -> CellVectorT {
  return embedding_ * input + offset_;
}

auto ReferenceFaceMap::cellToFace(const CellVectorT& input) const -> FaceVectorT {
  return projection_ * input + projectionOffset_;
}

FaceTransform::~FaceTransform() = default;

auto FaceTransform::normal(const FaceVectorT& input) const -> VectorT {
  const auto jacobian = refToSpaceJacobian(input);
  return VectorT(jacobian.col(0)).cross(VectorT(jacobian.col(1)));
}

auto FaceTransform::surfaceJacobian(const FaceVectorT& input) const -> double {
  return normal(input).norm();
}

auto FaceTransform::area() const -> double {
  return 0.5 * surfaceJacobian(FaceVectorT(Face::ReferenceBarycenter.data()));
}

auto FaceTransform::center() const -> VectorT {
  return refToSpace(FaceVectorT(Face::ReferenceBarycenter.data()));
}

auto FaceTransform::vertices() const -> std::array<VectorT, Face::NumVertices> {
  std::array<VectorT, Face::NumVertices> result{};
  result[0] = refToSpace(FaceVectorT(0.0, 0.0));
  result[1] = refToSpace(FaceVectorT(1.0, 0.0));
  result[2] = refToSpace(FaceVectorT(0.0, 1.0));
  return result;
}

auto FaceTransform::faceAlignedBasis() const -> std::array<VectorT, Cell::Dim> {
  const auto corners = vertices();
  const VectorT faceNormal = normal(FaceVectorT(Face::ReferenceBarycenter.data()));
  const VectorT tangent1 = corners[1] - corners[0];
  return {faceNormal, tangent1, VectorT(faceNormal.cross(tangent1))};
}

AffineFaceTransform::AffineFaceTransform(const AffineTransform& cell,
                                         const ReferenceFaceMap& embedding)
    : map_(embedding), cell_(cell) {
  transform_ = cell.jacobian() * embedding.embedding();
  offset_ = cell.refToSpace(embedding.faceToCell(FaceVectorT::Zero()));
}

auto AffineFaceTransform::fromMeshCell(std::size_t id,
                                       std::size_t side,
                                       const MeshReader& mesh,
                                       FaceOrientation orientation) -> AffineFaceTransform {
  return {AffineTransform::fromMeshCell(id, mesh), ReferenceFaceMap(side, orientation)};
}

auto AffineFaceTransform::fromMeshFault(std::size_t faultId,
                                        const MeshReader& mesh) -> AffineFaceTransform {
  const auto& fault = mesh.getFault()[faultId];
  if (fault.element.hasValue()) {
    return fromMeshCell(fault.element.value(), fault.side, mesh);
  }
  if (!fault.neighborElement.hasValue()) {
    logError() << "Fault face" << faultId << "has no adjacent cell on this rank.";
  }
  const auto neighbor = fault.neighborElement.value();
  const auto orientation = static_cast<FaceOrientation>(
      mesh.getElements()[neighbor].sideOrientations[fault.neighborSide]);
  return fromMeshCell(neighbor, fault.neighborSide, mesh, orientation);
}

auto AffineFaceTransform::refToSpace(const FaceVectorT& input) const -> VectorT {
  return transform_ * input + offset_;
}

auto AffineFaceTransform::refToCell(const FaceVectorT& input) const -> VectorT {
  return map_.faceToCell(input);
}

auto AffineFaceTransform::refToSpaceJacobian(const FaceVectorT& /*input*/) const -> JacobianT {
  // since we're linear—no dependency on the input vector here
  return transform_;
}

} // namespace seissol::geometry
