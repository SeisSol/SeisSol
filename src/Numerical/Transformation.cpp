// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Transformation.h"
#include <Eigen/Dense>
#include <Geometry/MeshDefinition.h>
#include <Kernels/Precision.h>
#include <cassert>
#include <cmath>
#include <utils/logger.h>
#include <yateto.h>

void seissol::transformations::tetrahedronReferenceToGlobal(const double v0[3],
                                                            const double v1[3],
                                                            const double v2[3],
                                                            const double v3[3],
                                                            const double xiEtaZeta[3],
                                                            double xyz[3]) {
  for (unsigned i = 0; i < 3; ++i) {
    xyz[i] = v0[i] + (v1[i] - v0[i]) * xiEtaZeta[0] + (v2[i] - v0[i]) * xiEtaZeta[1] +
             (v3[i] - v0[i]) * xiEtaZeta[2];
  }
}

Eigen::Vector3d seissol::transformations::tetrahedronReferenceToGlobal(const Eigen::Vector3d& v0,
                                                                       const Eigen::Vector3d& v1,
                                                                       const Eigen::Vector3d& v2,
                                                                       const Eigen::Vector3d& v3,
                                                                       const double xiEtaZeta[3]) {
  return v0 + (v1 - v0) * xiEtaZeta[0] + (v2 - v0) * xiEtaZeta[1] + (v3 - v0) * xiEtaZeta[2];
}

Eigen::Vector3d seissol::transformations::tetrahedronGlobalToReference(const double v0[3],
                                                                       const double v1[3],
                                                                       const double v2[3],
                                                                       const double v3[3],
                                                                       const Eigen::Vector3d& xyz) {
  // Forward transformation
  Eigen::Matrix4d a;
  a << v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2], 0.0, v2[0] - v0[0], v2[1] - v0[1],
      v2[2] - v0[2], 0.0, v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2], 0.0, v0[0], v0[1], v0[2],
      1.0;
  a = a.transpose().eval();

  const Eigen::Vector4d rhs(xyz[0], xyz[1], xyz[2], 1.0);
  Eigen::Vector4d sol = a.partialPivLu().solve(rhs);
  Eigen::Vector3d xiEtaZeta(sol[0], sol[1], sol[2]);

  return xiEtaZeta;
}

void seissol::transformations::tetrahedronGlobalToReferenceJacobian(const real iX[4],
                                                                    const real iY[4],
                                                                    const real iZ[4],
                                                                    real oGradXi[3],
                                                                    real oGradEta[3],
                                                                    real oGradZeta[3]) {
  const real determinant =
      iX[0] * (iY[1] * (iZ[3] - iZ[2]) + iY[2] * (iZ[1] - iZ[3]) + iY[3] * (iZ[2] - iZ[1])) +
      iX[1] * (iY[0] * (iZ[2] - iZ[3]) + iY[2] * (iZ[3] - iZ[0]) + iY[3] * (iZ[0] - iZ[2])) +
      iX[2] * (iY[0] * (iZ[3] - iZ[1]) + iY[1] * (iZ[0] - iZ[3]) + iY[3] * (iZ[1] - iZ[0])) +
      iX[3] * (iY[0] * (iZ[1] - iZ[2]) + iY[1] * (iZ[2] - iZ[0]) + iY[2] * (iZ[0] - iZ[1]));
  const real inverseDeterminant = 1.0 / determinant;

  // dxi_dx, dxi_dy, dxi_dz
  oGradXi[0] = inverseDeterminant *
               (iY[0] * (iZ[2] - iZ[3]) + iY[2] * (iZ[3] - iZ[0]) + iY[3] * (iZ[0] - iZ[2]));
  oGradXi[1] = inverseDeterminant *
               (iX[0] * (iZ[3] - iZ[2]) + iX[2] * (iZ[0] - iZ[3]) + iX[3] * (iZ[2] - iZ[0]));
  oGradXi[2] = inverseDeterminant *
               (iX[0] * (iY[2] - iY[3]) + iX[2] * (iY[3] - iY[0]) + iX[3] * (iY[0] - iY[2]));

  // deta_dx, deta_dy, deta_dz
  oGradEta[0] = inverseDeterminant *
                (iY[0] * (iZ[3] - iZ[1]) + iY[1] * (iZ[0] - iZ[3]) + iY[3] * (iZ[1] - iZ[0]));
  oGradEta[1] = inverseDeterminant *
                (iX[0] * (iZ[1] - iZ[3]) + iX[1] * (iZ[3] - iZ[0]) + iX[3] * (iZ[0] - iZ[1]));
  oGradEta[2] = inverseDeterminant *
                (iX[0] * (iY[3] - iY[1]) + iX[1] * (iY[0] - iY[3]) + iX[3] * (iY[1] - iY[0]));

  // dzeta_dx, dzeta_dx, dzeta_dx
  oGradZeta[0] = inverseDeterminant *
                 (iY[0] * (iZ[1] - iZ[2]) + iY[1] * (iZ[2] - iZ[0]) + iY[2] * (iZ[0] - iZ[1]));
  oGradZeta[1] = inverseDeterminant *
                 (iX[0] * (iZ[2] - iZ[1]) + iX[1] * (iZ[0] - iZ[2]) + iX[2] * (iZ[1] - iZ[0]));
  oGradZeta[2] = inverseDeterminant *
                 (iX[0] * (iY[1] - iY[2]) + iX[1] * (iY[2] - iY[0]) + iX[2] * (iY[0] - iY[1]));
}

void seissol::transformations::inverseTensor1RotationMatrix(
    const VrtxCoords iNormal,
    const VrtxCoords iTangent1,
    const VrtxCoords iTangent2,
    yateto::DenseTensorView<2, real, unsigned>& oTinv,
    unsigned row,
    unsigned col) {
  for (unsigned i = 0; i < 3; ++i) {
    oTinv(row + 0, col + i) = iNormal[i];
    oTinv(row + 1, col + i) = iTangent1[i];
    oTinv(row + 2, col + i) = iTangent2[i];
  }
}

void seissol::transformations::tensor1RotationMatrix(const VrtxCoords iNormal,
                                                     const VrtxCoords iTangent1,
                                                     const VrtxCoords iTangent2,
                                                     yateto::DenseTensorView<2, real, unsigned>& oT,
                                                     unsigned row,
                                                     unsigned col) {
  for (unsigned i = 0; i < 3; ++i) {
    oT(row + i, col + 0) = iNormal[i];
    oT(row + i, col + 1) = iTangent1[i];
    oT(row + i, col + 2) = iTangent2[i];
  }
}

void seissol::transformations::symmetricTensor2RotationMatrix(
    const VrtxCoords iNormal,
    const VrtxCoords iTangent1,
    const VrtxCoords iTangent2,
    yateto::DenseTensorView<2, real, unsigned>& oT,
    unsigned row,
    unsigned col) {
  const real nx = iNormal[0];
  const real ny = iNormal[1];
  const real nz = iNormal[2];
  const real sx = iTangent1[0];
  const real sy = iTangent1[1];
  const real sz = iTangent1[2];
  const real tx = iTangent2[0];
  const real ty = iTangent2[1];
  const real tz = iTangent2[2];

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

void seissol::transformations::chiTau2XiEtaZeta(unsigned face,
                                                const double chiTau[2],
                                                double xiEtaZeta[3],
                                                int sideOrientation) {
  double chiTauTilde[2];

  switch (sideOrientation) {
  case 0:
    chiTauTilde[0] = chiTau[1];
    chiTauTilde[1] = chiTau[0];
    break;
  case 1:
    chiTauTilde[0] = 1.0 - chiTau[0] - chiTau[1];
    chiTauTilde[1] = chiTau[1];
    break;
  case 2:
    chiTauTilde[0] = chiTau[0];
    chiTauTilde[1] = 1.0 - chiTau[0] - chiTau[1];
    break;
  default:
    chiTauTilde[0] = chiTau[0];
    chiTauTilde[1] = chiTau[1];
    break;
  }

  switch (face) {
  case 0:
    xiEtaZeta[0] = chiTauTilde[1];
    xiEtaZeta[1] = chiTauTilde[0];
    xiEtaZeta[2] = 0.0;
    break;
  case 1:
    xiEtaZeta[0] = chiTauTilde[0];
    xiEtaZeta[1] = 0.0;
    xiEtaZeta[2] = chiTauTilde[1];
    break;
  case 2:
    xiEtaZeta[0] = 0.0;
    xiEtaZeta[1] = chiTauTilde[1];
    xiEtaZeta[2] = chiTauTilde[0];
    break;
  case 3:
    xiEtaZeta[0] = 1.0 - chiTauTilde[0] - chiTauTilde[1];
    xiEtaZeta[1] = chiTauTilde[0];
    xiEtaZeta[2] = chiTauTilde[1];
    break;
  default:
    break;
  }
}

void seissol::transformations::XiEtaZeta2chiTau(unsigned face,
                                                const double xiEtaZeta[3],
                                                double chiTau[2]) {
  [[maybe_unused]] constexpr double Eps = 1e-6;

  switch (face) {
  case 0: {
    chiTau[1] = xiEtaZeta[0];
    chiTau[0] = xiEtaZeta[1];
    assert((std::abs(xiEtaZeta[2]) < Eps) && "reference coord is not on the 1st face");
    break;
  }
  case 1: {
    chiTau[0] = xiEtaZeta[0];
    chiTau[1] = xiEtaZeta[2];
    assert((std::abs(xiEtaZeta[1]) < Eps) && "reference coord is not on the 2nd face");
    break;
  }
  case 2: {
    chiTau[1] = xiEtaZeta[1];
    chiTau[0] = xiEtaZeta[2];
    assert((std::abs(xiEtaZeta[0]) < Eps) && "reference coord is not on the 3rd face");
    break;
  }
  case 3: {
    chiTau[0] = xiEtaZeta[1];
    chiTau[1] = xiEtaZeta[2];
    assert((std::abs(xiEtaZeta[0] + xiEtaZeta[1] + xiEtaZeta[2] - 1.0) < Eps) &&
           "reference coord is not on the 4th face");
    break;
  }
  default:
    logError() << "Tried to get the XiEtaZeta2chiTau transformation for face" << face
               << ", which is not possible. Provide 0 <= face <= 3.";
  }
}
