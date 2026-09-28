// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_NUMERICAL_TRANSFORMATION_H_
#define SEISSOL_SRC_NUMERICAL_TRANSFORMATION_H_

#include "Geometry/MeshDefinition.h"
#include "Initializer/Typedefs.h"

#include <Eigen/Dense>
#include <yateto.h>

namespace seissol::transformations {
/**
 * Inverse of Tensor1RotationMatrix().
 **/
template <typename RealT>
void inverseTensor1RotationMatrix(const CoordinateT& iNormal,
                                  const CoordinateT& iTangent1,
                                  const CoordinateT& iTangent2,
                                  yateto::DenseTensorView<2, RealT, unsigned>& oTinv,
                                  std::uint32_t row = 0,
                                  std::uint32_t col = 0);

/**
 * Returns a column-major matrix that rotates a first-order tensor
 * (u_x, u_y, u_y) into a new coordinate system aligned with normal
 * and tangents.
 * u' = T*u
 **/
template <typename RealT>
void tensor1RotationMatrix(const CoordinateT& iNormal,
                           const CoordinateT& iTangent1,
                           const CoordinateT& iTangent2,
                           yateto::DenseTensorView<2, RealT, unsigned>& oT,
                           std::uint32_t row = 0,
                           std::uint32_t col = 0);

/**
 * Inverse of SymmetricTensor2RotationMatrix().
 **/
template <typename RealT>
void inverseSymmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                           const CoordinateT& iTangent1,
                                           const CoordinateT& iTangent2,
                                           yateto::DenseTensorView<2, RealT, unsigned>& oTinv,
                                           std::uint32_t row = 0,
                                           std::uint32_t col = 0);

/**
 * Returns a column-major matrix that rotates a symmetric second-order
 * tensor, given as vector (u_xx, u_yy, u_zz, u_xy, u_yz, u_xz),
 * into a new coordinate system aligned with normal and tangents.
 * u' = T*u
 **/
template <typename RealT>
void symmetricTensor2RotationMatrix(const CoordinateT& iNormal,
                                    const CoordinateT& iTangent1,
                                    const CoordinateT& iTangent2,
                                    yateto::DenseTensorView<2, RealT, unsigned>& oTinv,
                                    std::uint32_t row = 0,
                                    std::uint32_t col = 0);

} // namespace seissol::transformations

#endif // SEISSOL_SRC_NUMERICAL_TRANSFORMATION_H_
