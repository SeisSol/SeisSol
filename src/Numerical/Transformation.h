// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
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
void tetrahedronGlobalToReferenceJacobian(const real iX[4],
                                          const real iY[4],
                                          const real iZ[4],
                                          real oGradXi[3],
                                          real oGradEta[3],
                                          real oGradZeta[3]);

/**
 * Inverse of Tensor1RotationMatrix().
 **/
void inverseTensor1RotationMatrix(const VrtxCoords iNormal,
                                  const VrtxCoords iTangent1,
                                  const VrtxCoords iTangent2,
                                  yateto::DenseTensorView<2, real, unsigned>& oTinv,
                                  unsigned row = 0,
                                  unsigned col = 0);

/**
 * Returns a column-major matrix that rotates a first-order tensor
 * (u_x, u_y, u_y) into a new coordinate system aligned with normal
 * and tangents.
 * u' = T*u
 **/
void tensor1RotationMatrix(const VrtxCoords iNormal,
                           const VrtxCoords iTangent1,
                           const VrtxCoords iTangent2,
                           yateto::DenseTensorView<2, real, unsigned>& oT,
                           unsigned row = 0,
                           unsigned col = 0);

/**
 * Inverse of SymmetricTensor2RotationMatrix().
 **/
template <typename Tmatrix>
void inverseSymmetricTensor2RotationMatrix(const VrtxCoords iNormal,
                                           const VrtxCoords iTangent1,
                                           const VrtxCoords iTangent2,
                                           Tmatrix& oTinv,
                                           unsigned row = 0,
                                           unsigned col = 0) {
  const real nx = iNormal[0];
  const real ny = iNormal[1];
  const real nz = iNormal[2];
  const real sx = iTangent1[0];
  const real sy = iTangent1[1];
  const real sz = iTangent1[2];
  const real tx = iTangent2[0];
  const real ty = iTangent2[1];
  const real tz = iTangent2[2];

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

/**
 * Returns a column-major matrix that rotates a symmetric second-order
 * tensor, given as vector (u_xx, u_yy, u_zz, u_xy, u_yz, u_xz),
 * into a new coordinate system aligned with normal and tangents.
 * u' = T*u
 **/
void symmetricTensor2RotationMatrix(const VrtxCoords iNormal,
                                    const VrtxCoords iTangent1,
                                    const VrtxCoords iTangent2,
                                    yateto::DenseTensorView<2, real, unsigned>& oTinv,
                                    unsigned row = 0,
                                    unsigned col = 0);

void chiTau2XiEtaZeta(unsigned face,
                      const double chiTau[2],
                      double xiEtaZeta[3],
                      int sideOrientation = -1);
void XiEtaZeta2chiTau(unsigned face, const double xiEtaZeta[3], double chiTau[2]);
} // namespace seissol::transformations

#endif // SEISSOL_SRC_NUMERICAL_TRANSFORMATION_H_
