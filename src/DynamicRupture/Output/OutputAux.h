// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTAUX_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTAUX_H_

#include "DataTypes.h"
#include "Geometry/MeshReader.h"

#include <array>
#include <cstddef>
#include <memory>

namespace seissol {
template <int N, typename T>
auto unsafe_reshape(T* ptr) -> T (*)[N] {
  return reinterpret_cast<T(*)[N]>(ptr);
}
} // namespace seissol

namespace seissol::dr {
int getElementVertexId(int localSideId, int localFaceVertexId);

ExtTriangle getReferenceTriangle(int sideIdx);

ExtTriangle getGlobalTriangle(int localSideId,
                              const Element& element,
                              const std::vector<Vertex>& verticesInfo);

ExtVrtxCoords getMidPointTriangle(const ExtTriangle& triangle);

ExtVrtxCoords getTrianglePointByCoords(const ExtTriangle& triangle,
                                       const std::array<double, 2>& point);

ExtVrtxCoords getMidPoint(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2);

struct TriangleQuadratureData {
  static constexpr size_t Size{tensor::quadweights::Shape[0]};
  std::array<double, 2 * Size> points{};
  std::array<double, Size> weights{};
};

TriangleQuadratureData generateTriangleQuadrature();

void assignNearestGaussianPoints(ReceiverPoints& geoPoints);

int getClosestInternalStroudGp(int nearestGpIndex, int nPoly);

std::pair<int, double> getNearestFacePoint(const double targetPoint[2],
                                           const double (*facePoints)[2],
                                           std::size_t numFacePoints);

double
    isInsideFace(const ExtVrtxCoords& point, const ExtTriangle& face, const VrtxCoords faceNormal);

void projectPointToFace(ExtVrtxCoords& point, const ExtTriangle& face, const VrtxCoords faceNormal);

double getDistanceFromPointToFace(const ExtVrtxCoords& point,
                                  const ExtTriangle& face,
                                  const VrtxCoords faceNormal);

PlusMinusBasisFunctions getPlusMinusBasisFunctions(const VrtxCoords point,
                                                   const VrtxCoords* plusElementCoords[4],
                                                   const VrtxCoords* minusElementCoords[4]);

real computeTriangleArea(ExtTriangle& triangle);

/**
 * @brief The index of the first receiver an output cell owns.
 *
 * The refiner emits the receivers of a cell point-major and simulation-minor, so a cell owns
 * @p pointsPerCell * @p simulationCount consecutive entries. Properties that every receiver of
 * the cell shares are read off the first of them.
 */
std::size_t
    firstReceiverOfCell(std::size_t cell, std::size_t pointsPerCell, std::size_t simulationCount);

/**
 * @brief The fault tag of an output cell, as the mesh assigned it to the face.
 *
 * This is the group the face was tagged with in the mesh file, not an identifier: several faces
 * carry the same tag, and a mesh that tags nothing leaves it at its default.
 */
int faultTagOfCell(const ReceiverPoints& receiverPoints,
                   std::size_t cell,
                   std::size_t pointsPerCell,
                   std::size_t simulationCount);

/**
 * @brief The global identifier of the face an output cell sits on.
 *
 * Unique across the mesh, since it is built from the global element index and the side.
 */
std::size_t globalFaceIdOfCell(const ReceiverPoints& receiverPoints,
                               std::size_t cell,
                               std::size_t pointsPerCell,
                               std::size_t simulationCount);
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTAUX_H_
