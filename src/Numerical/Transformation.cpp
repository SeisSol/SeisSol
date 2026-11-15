// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Transformation.h"

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>
#include <cstdint>
#include <yateto.h>

namespace seissol::transformations {

template <typename RealT>
void inverseTensor1RotationMatrix(const CoordinateT& iNormal,
                                  const CoordinateT& iTangent1,
                                  const CoordinateT& iTangent2,
                                  yateto::DenseTensorView<2, RealT, unsigned>& oTinv,
                                  std::uint32_t row,
                                  std::uint32_t col) {
  for (std::uint32_t i = 0; i < Cell::Dim; ++i) {
    oTinv(row + 0, col + i) = iNormal[i];
    oTinv(row + 1, col + i) = iTangent1[i];
    oTinv(row + 2, col + i) = iTangent2[i];
  }
}

template void inverseTensor1RotationMatrix(const CoordinateT& iNormal,
                                           const CoordinateT& iTangent1,
                                           const CoordinateT& iTangent2,
                                           yateto::DenseTensorView<2, float, unsigned>& oTinv,
                                           std::uint32_t row,
                                           std::uint32_t col);

template void inverseTensor1RotationMatrix(const CoordinateT& iNormal,
                                           const CoordinateT& iTangent1,
                                           const CoordinateT& iTangent2,
                                           yateto::DenseTensorView<2, double, unsigned>& oTinv,
                                           std::uint32_t row,
                                           std::uint32_t col);

template <typename RealT>
void tensor1RotationMatrix(const CoordinateT& iNormal,
                           const CoordinateT& iTangent1,
                           const CoordinateT& iTangent2,
                           yateto::DenseTensorView<2, RealT, unsigned>& oT,
                           std::uint32_t row,
                           std::uint32_t col) {
  for (std::uint32_t i = 0; i < Cell::Dim; ++i) {
    oT(row + i, col + 0) = iNormal[i];
    oT(row + i, col + 1) = iTangent1[i];
    oT(row + i, col + 2) = iTangent2[i];
  }
}

template void tensor1RotationMatrix(const CoordinateT& iNormal,
                                    const CoordinateT& iTangent1,
                                    const CoordinateT& iTangent2,
                                    yateto::DenseTensorView<2, float, unsigned>& oT,
                                    std::uint32_t row,
                                    std::uint32_t col);

template void tensor1RotationMatrix(const CoordinateT& iNormal,
                                    const CoordinateT& iTangent1,
                                    const CoordinateT& iTangent2,
                                    yateto::DenseTensorView<2, double, unsigned>& oT,
                                    std::uint32_t row,
                                    std::uint32_t col);

template <typename RealT>
void symmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                    const CoordinateT& iTangent1,
                                    const CoordinateT& iTangent2,
                                    yateto::DenseTensorView<2, RealT, unsigned>& oT,
                                    std::uint32_t row,
                                    std::uint32_t col) {
  const auto nx = iNormal[0];
  const auto ny = iNormal[1];
  const auto nz = iNormal[2];
  const auto sx = iTangent1[0];
  const auto sy = iTangent1[1];
  const auto sz = iTangent1[2];
  const auto tx = iTangent2[0];
  const auto ty = iTangent2[1];
  const auto tz = iTangent2[2];

  oT(row + 0, col + 0) = nx * nx;
  oT(row + 1, col + 0) = ny * ny;
  oT(row + 2, col + 0) = nz * nz;
  oT(row + 3, col + 0) = ny * nx;
  oT(row + 4, col + 0) = nz * ny;
  oT(row + 5, col + 0) = nz * nx;
  oT(row + 0, col + 1) = sx * sx;
  oT(row + 1, col + 1) = sy * sy;
  oT(row + 2, col + 1) = sz * sz;
  oT(row + 3, col + 1) = sy * sx;
  oT(row + 4, col + 1) = sz * sy;
  oT(row + 5, col + 1) = sz * sx;
  oT(row + 0, col + 2) = tx * tx;
  oT(row + 1, col + 2) = ty * ty;
  oT(row + 2, col + 2) = tz * tz;
  oT(row + 3, col + 2) = ty * tx;
  oT(row + 4, col + 2) = tz * ty;
  oT(row + 5, col + 2) = tz * tx;
  oT(row + 0, col + 3) = 2.0 * nx * sx;
  oT(row + 1, col + 3) = 2.0 * ny * sy;
  oT(row + 2, col + 3) = 2.0 * nz * sz;
  oT(row + 3, col + 3) = ny * sx + nx * sy;
  oT(row + 4, col + 3) = nz * sy + ny * sz;
  oT(row + 5, col + 3) = nz * sx + nx * sz;
  oT(row + 0, col + 4) = 2.0 * sx * tx;
  oT(row + 1, col + 4) = 2.0 * sy * ty;
  oT(row + 2, col + 4) = 2.0 * sz * tz;
  oT(row + 3, col + 4) = sy * tx + sx * ty;
  oT(row + 4, col + 4) = sz * ty + sy * tz;
  oT(row + 5, col + 4) = sz * tx + sx * tz;
  oT(row + 0, col + 5) = 2.0 * nx * tx;
  oT(row + 1, col + 5) = 2.0 * ny * ty;
  oT(row + 2, col + 5) = 2.0 * nz * tz;
  oT(row + 3, col + 5) = ny * tx + nx * ty;
  oT(row + 4, col + 5) = nz * ty + ny * tz;
  oT(row + 5, col + 5) = nz * tx + nx * tz;
}

template void symmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                             const CoordinateT& iTangent1,
                                             const CoordinateT& iTangent2,
                                             yateto::DenseTensorView<2, float, unsigned>& oT,
                                             std::uint32_t row,
                                             std::uint32_t col);

template void symmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                             const CoordinateT& iTangent1,
                                             const CoordinateT& iTangent2,
                                             yateto::DenseTensorView<2, double, unsigned>& oT,
                                             std::uint32_t row,
                                             std::uint32_t col);

template <typename RealT>
void inverseSymmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                           const CoordinateT& iTangent1,
                                           const CoordinateT& iTangent2,
                                           yateto::DenseTensorView<2, RealT, unsigned>& oTinv,
                                           std::uint32_t row,
                                           std::uint32_t col) {
  const auto nx = iNormal[0];
  const auto ny = iNormal[1];
  const auto nz = iNormal[2];
  const auto sx = iTangent1[0];
  const auto sy = iTangent1[1];
  const auto sz = iTangent1[2];
  const auto tx = iTangent2[0];
  const auto ty = iTangent2[1];
  const auto tz = iTangent2[2];

  oTinv(row + 0, col + 0) = nx * nx;
  oTinv(row + 1, col + 0) = sx * sx;
  oTinv(row + 2, col + 0) = tx * tx;
  oTinv(row + 3, col + 0) = nx * sx;
  oTinv(row + 4, col + 0) = sx * tx;
  oTinv(row + 5, col + 0) = nx * tx;
  oTinv(row + 0, col + 1) = ny * ny;
  oTinv(row + 1, col + 1) = sy * sy;
  oTinv(row + 2, col + 1) = ty * ty;
  oTinv(row + 3, col + 1) = ny * sy;
  oTinv(row + 4, col + 1) = sy * ty;
  oTinv(row + 5, col + 1) = ny * ty;
  oTinv(row + 0, col + 2) = nz * nz;
  oTinv(row + 1, col + 2) = sz * sz;
  oTinv(row + 2, col + 2) = tz * tz;
  oTinv(row + 3, col + 2) = nz * sz;
  oTinv(row + 4, col + 2) = sz * tz;
  oTinv(row + 5, col + 2) = nz * tz;
  oTinv(row + 0, col + 3) = 2.0 * ny * nx;
  oTinv(row + 1, col + 3) = 2.0 * sy * sx;
  oTinv(row + 2, col + 3) = 2.0 * ty * tx;
  oTinv(row + 3, col + 3) = ny * sx + nx * sy;
  oTinv(row + 4, col + 3) = sy * tx + sx * ty;
  oTinv(row + 5, col + 3) = ny * tx + nx * ty;
  oTinv(row + 0, col + 4) = 2.0 * nz * ny;
  oTinv(row + 1, col + 4) = 2.0 * sz * sy;
  oTinv(row + 2, col + 4) = 2.0 * tz * ty;
  oTinv(row + 3, col + 4) = nz * sy + ny * sz;
  oTinv(row + 4, col + 4) = sz * ty + sy * tz;
  oTinv(row + 5, col + 4) = nz * ty + ny * tz;
  oTinv(row + 0, col + 5) = 2.0 * nz * nx;
  oTinv(row + 1, col + 5) = 2.0 * sz * sx;
  oTinv(row + 2, col + 5) = 2.0 * tz * tx;
  oTinv(row + 3, col + 5) = nz * sx + nx * sz;
  oTinv(row + 4, col + 5) = sz * tx + sx * tz;
  oTinv(row + 5, col + 5) = nz * tx + nx * tz;
}

template void
    inverseSymmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                          const CoordinateT& iTangent1,
                                          const CoordinateT& iTangent2,
                                          yateto::DenseTensorView<2, float, unsigned>& oTinv,
                                          std::uint32_t row,
                                          std::uint32_t col);

template void
    inverseSymmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                          const CoordinateT& iTangent1,
                                          const CoordinateT& iTangent2,
                                          yateto::DenseTensorView<2, double, unsigned>& oTinv,
                                          std::uint32_t row,
                                          std::uint32_t col);

} // namespace seissol::transformations
