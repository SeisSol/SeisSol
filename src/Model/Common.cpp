// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "Model/Common.h"
#include <Geometry/MeshDefinition.h>
#include <Kernels/Precision.h>
#include <Numerical/Transformation.h>
#include <cmath>
#include <generated_code/init.h>
#include <limits>

bool seissol::model::testIfAcoustic(real mu) {
  return std::abs(mu) <= std::numeric_limits<real>::epsilon();
}

void seissol::model::getBondMatrix(const VrtxCoords normal,
                                   const VrtxCoords tangent1,
                                   const VrtxCoords tangent2,
                                   real* matN) {
  matN[0 * 6 + 0] = normal[0] * normal[0];
  matN[0 * 6 + 1] = normal[1] * normal[1];
  matN[0 * 6 + 2] = normal[2] * normal[2];
  matN[0 * 6 + 3] = 2 * normal[2] * normal[1];
  matN[0 * 6 + 4] = 2 * normal[2] * normal[0];
  matN[0 * 6 + 5] = 2 * normal[1] * normal[0];
  matN[1 * 6 + 0] = tangent1[0] * tangent1[0];
  matN[1 * 6 + 1] = tangent1[1] * tangent1[1];
  matN[1 * 6 + 2] = tangent1[2] * tangent1[2];
  matN[1 * 6 + 3] = 2 * tangent1[2] * tangent1[1];
  matN[1 * 6 + 4] = 2 * tangent1[2] * tangent1[0];
  matN[1 * 6 + 5] = 2 * tangent1[1] * tangent1[0];
  matN[2 * 6 + 0] = tangent2[0] * tangent2[0];
  matN[2 * 6 + 1] = tangent2[1] * tangent2[1];
  matN[2 * 6 + 2] = tangent2[2] * tangent2[2];
  matN[2 * 6 + 3] = 2 * tangent2[2] * tangent2[1];
  matN[2 * 6 + 4] = 2 * tangent2[2] * tangent2[0];
  matN[2 * 6 + 5] = 2 * tangent2[1] * tangent2[0];

  matN[3 * 6 + 0] = tangent1[0] * tangent2[0];
  matN[3 * 6 + 1] = tangent1[1] * tangent2[1];
  matN[3 * 6 + 2] = tangent1[2] * tangent2[2];
  matN[3 * 6 + 3] = tangent1[1] * tangent2[2] + tangent1[2] * tangent2[1];
  matN[3 * 6 + 4] = tangent1[0] * tangent2[2] + tangent1[2] * tangent2[0];
  matN[3 * 6 + 5] = tangent1[1] * tangent2[0] + tangent1[0] * tangent2[1];
  matN[4 * 6 + 0] = normal[0] * tangent2[0];
  matN[4 * 6 + 1] = normal[1] * tangent2[1];
  matN[4 * 6 + 2] = normal[2] * tangent2[2];
  matN[4 * 6 + 3] = normal[1] * tangent2[2] + normal[2] * tangent2[1];
  matN[4 * 6 + 4] = normal[0] * tangent2[2] + normal[2] * tangent2[0];
  matN[4 * 6 + 5] = normal[1] * tangent2[0] + normal[0] * tangent2[1];
  matN[5 * 6 + 0] = normal[0] * tangent1[0];
  matN[5 * 6 + 1] = normal[1] * tangent1[1];
  matN[5 * 6 + 2] = normal[2] * tangent1[2];
  matN[5 * 6 + 3] = normal[1] * tangent1[2] + normal[2] * tangent1[1];
  matN[5 * 6 + 4] = normal[0] * tangent1[2] + normal[2] * tangent1[0];
  matN[5 * 6 + 5] = normal[1] * tangent1[0] + normal[0] * tangent1[1];
}
void seissol::model::getFaceRotationMatrix(const Eigen::Vector3d& normal,
                                           const Eigen::Vector3d& tangent1,
                                           const Eigen::Vector3d& tangent2,
                                           init::T::view::type& matT,
                                           init::Tinv::view::type& matTinv) {
  const VrtxCoords n = {normal(0), normal(1), normal(2)};
  const VrtxCoords s = {tangent1(0), tangent1(1), tangent1(2)};
  const VrtxCoords t = {tangent2(0), tangent2(1), tangent2(2)};
  getFaceRotationMatrix(n, s, t, matT, matTinv);
}

void seissol::model::getFaceRotationMatrix(const VrtxCoords normal,
                                           const VrtxCoords tangent1,
                                           const VrtxCoords tangent2,
                                           init::T::view::type& matT,
                                           init::Tinv::view::type& matTinv) {
  matT.setZero();
  matTinv.setZero();

  seissol::transformations::symmetricTensor2RotationMatrix(normal, tangent1, tangent2, matT, 0, 0);
  seissol::transformations::tensor1RotationMatrix(normal, tangent1, tangent2, matT, 6, 6);

  seissol::transformations::inverseSymmetricTensor2RotationMatrix(
      normal, tangent1, tangent2, matTinv, 0, 0);
  seissol::transformations::inverseTensor1RotationMatrix(normal, tangent1, tangent2, matTinv, 6, 6);

#ifdef USE_VISCOELASTIC
  for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
    const unsigned origin = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
    seissol::transformations::symmetricTensor2RotationMatrix(
        normal, tangent1, tangent2, matT, origin, origin);
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        normal, tangent1, tangent2, matTinv, origin, origin);
  }
#elif USE_VISCOELASTIC2
  seissol::transformations::symmetricTensor2RotationMatrix(normal,
                                                           tangent1,
                                                           tangent2,
                                                           matT,
                                                           MaterialT::NumElasticQuantities,
                                                           MaterialT::NumElasticQuantities);
#elif USE_POROELASTIC
  // pressure
  matT(9, 9) = 1;
  matTinv(9, 9) = 1;
  // fluid velocities
  unsigned origin = 10;
  seissol::transformations::tensor1RotationMatrix(normal, tangent1, tangent2, matT, origin, origin);
  seissol::transformations::inverseTensor1RotationMatrix(
      normal, tangent1, tangent2, matTinv, origin, origin);
#endif
}
