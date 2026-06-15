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
 * Calculates the global coordinates from
 * reference tetrahedron coordinates.
 */
void tetrahedronReferenceToGlobal(const double v0[3],
                                  const double v1[3],
                                  const double v2[3],
                                  const double v3[3],
                                  const double xiEtaZeta[3],
                                  double xyz[3]);

Eigen::Vector3d tetrahedronReferenceToGlobal(const Eigen::Vector3d& v0,
                                             const Eigen::Vector3d& v1,
                                             const Eigen::Vector3d& v2,
                                             const Eigen::Vector3d& v3,
                                             const double xiEtaZeta[3]);

/**
 * Calculates the reference tetrahedron coordinates from
 * global tetrahedron coordinates.
 */
Eigen::Vector3d tetrahedronGlobalToReference(const double v0[3],
                                             const double v1[3],
                                             const double v2[3],
                                             const double v3[3],
                                             const Eigen::Vector3d& xyz);

/**
 * Calculates the Jacobian for the coordinate transformation
 * xi(x, y, z), eta(x, y, z), zeta(x, y, z)
 * from a global tetrahedron to the reference tetrahedron.
 **/
void tetrahedronGlobalToReferenceJacobian(const double iX[4],
                                          const double iY[4],
                                          const double iZ[4],
                                          double oGradXi[3],
                                          double oGradEta[3],
                                          double oGradZeta[3]);

/**
 * Inverse of Tensor1RotationMatrix().
 **/
template <typename RealT>
void inverseTensor1RotationMatrix(const VrtxCoords iNormal,
                                  const VrtxCoords iTangent1,
                                  const VrtxCoords iTangent2,
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
void tensor1RotationMatrix(const VrtxCoords iNormal,
                           const VrtxCoords iTangent1,
                           const VrtxCoords iTangent2,
                           yateto::DenseTensorView<2, RealT, unsigned>& oT,
                           std::uint32_t row = 0,
                           std::uint32_t col = 0);

/**
 * Inverse of SymmetricTensor2RotationMatrix().
 **/
template <typename RealT>
void inverseSymmetricTensor2RotationMatrix(const VrtxCoords iNormal,
                                           const VrtxCoords iTangent1,
                                           const VrtxCoords iTangent2,
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
void symmetricTensor2RotationMatrix(const VrtxCoords iNormal,
                                    const VrtxCoords iTangent1,
                                    const VrtxCoords iTangent2,
                                    yateto::DenseTensorView<2, RealT, unsigned>& oTinv,
                                    std::uint32_t row = 0,
                                    std::uint32_t col = 0);

void chiTau2XiEtaZeta(std::uint32_t face,
                      const double chiTau[2],
                      double xiEtaZeta[3],
                      std::int32_t sideOrientation = -1);
void XiEtaZeta2chiTau(std::uint32_t face, const double xiEtaZeta[3], double chiTau[2]);
} // namespace seissol::transformations

#endif // SEISSOL_SRC_NUMERICAL_TRANSFORMATION_H_
