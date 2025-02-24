// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
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

ExtVrtxCoords getMidPoint(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2);

struct TriangleQuadratureData {
  static constexpr size_t Size{tensor::quadweights::Shape[0]};
  std::array<double, 2 * Size> points{};
  std::array<double, Size> weights{};
};

TriangleQuadratureData generateTriangleQuadrature(unsigned polyDegree);

void assignNearestGaussianPoints(ReceiverPoints& geoPoints);

int getClosestInternalStroudGp(int nearestGpIndex, int nPoly);

std::pair<int, double> getNearestFacePoint(const double targetPoint[2],
                                           const double (*facePoints)[2],
                                           unsigned numFacePoints);

void projectPointToFace(ExtVrtxCoords& point, const ExtTriangle& face, const VrtxCoords faceNormal);

double getDistanceFromPointToFace(const ExtVrtxCoords& point,
                                  const ExtTriangle& face,
                                  const VrtxCoords faceNormal);

PlusMinusBasisFunctions getPlusMinusBasisFunctions(const VrtxCoords point,
                                                   const VrtxCoords* plusElementCoords[4],
                                                   const VrtxCoords* minusElementCoords[4]);

std::vector<double> getAllVertices(const seissol::dr::ReceiverPoints& receiverPoints);

std::vector<unsigned int> getCellConnectivity(const seissol::dr::ReceiverPoints& receiverPoints);
std::vector<unsigned int> getFaultTags(const seissol::dr::ReceiverPoints& receiverPoints);

real computeTriangleArea(ExtTriangle& triangle);

template <int Size>
std::unique_ptr<int[]> convertMaskFromBoolToInt(const std::array<bool, Size>& boolMask) {
  auto intMask = std::unique_ptr<int[]>(new int[boolMask.size()]);

  for (size_t i = 0; i < boolMask.size(); ++i) {
    intMask[i] = static_cast<int>(boolMask[i]);
  }

  return intMask;
}
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_OUTPUTAUX_H_
