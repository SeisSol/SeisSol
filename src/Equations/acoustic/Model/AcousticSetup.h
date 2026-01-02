// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf
// SPDX-FileContributor: Jinwen Pan

#ifndef SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ACOUSTICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ACOUSTICSETUP_H_

#include "Equations/acoustic/Model/Datastructures.h"
#include "Equations/acoustic/Model/IntegrationData.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"

namespace seissol::model {
using Matrix44 = Eigen::Matrix<double, 4, 4>;

template <>
struct MaterialSetup<AcousticMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const AcousticMaterial& material, unsigned dim, T& matM) {
    matM.setZero();

    const double rhoInv = 1.0 / material.rho;

    switch (dim) {
    case 0:
      matM(1, 0) = material.lambda;
      matM(0, 1) = rhoInv;
      break;

    case 1:
      matM(2, 0) = material.lambda;
      matM(0, 2) = rhoInv;
      break;

    case 2:
      matM(3, 0) = material.lambda;
      matM(0, 3) = rhoInv;
      break;

    default:
      break;
    }
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const AcousticMaterial& local,
                                        const AcousticMaterial& neighbor,
                                        FaceType faceType,
                                        Tloc& qGodLocal,
                                        Tneigh& qGodNeighbor) {
    qGodNeighbor.setZero();

    // matR == eigenvector matrix

    // Eigenvectors are precomputed
    Matrix44 matR = Matrix44::Zero();
    // scale for better condition number of R
    matR(0, 0) = std::sqrt(local.lambda * local.rho);
    matR(1, 0) = -local.lambda;
    matR(0, 1) = std::sqrt(neighbor.lambda * neighbor.rho);
    matR(1, 1) = neighbor.lambda;
    matR(2, 2) = local.lambda;
    matR(3, 3) = local.lambda;

    if (faceType == FaceType::FreeSurface) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          qGodNeighbor(i, j) = std::numeric_limits<double>::signaling_NaN();
        }
      }
      qGodLocal.setZero();
      qGodLocal(0, 1) = -1 * matR(1, 0) * 1 / matR(0, 0);
      qGodLocal(1, 1) = 1.0;
    } else {
      Matrix44 chi = Matrix44::Zero();
      chi(0, 0) = 1.0;

      const auto godunov = ((matR * chi) * matR.inverse()).eval();

      // qGodLocal = I - qGodNeighbor
      for (unsigned i = 0; i < godunov.cols(); ++i) {
        for (unsigned j = 0; j < godunov.rows(); ++j) {
          qGodLocal(i, j) = -godunov(j, i);
          qGodNeighbor(i, j) = godunov(j, i);
        }
      }
      for (unsigned idx = 0; idx < 4; ++idx) {
        qGodLocal(idx, idx) += 1.0;
      }
    }
  }

  static AcousticMaterial
      getRotatedMaterialCoefficients(const std::array<double, 36>& /*rotationParameters*/,
                                     AcousticMaterial& material) {
    return material;
  }
  static void initializeSpecificLocalData(const AcousticMaterial& material,
                                          double timeStepWidth,
                                          AcousticLocalData* localData) {}

  static void initializeSpecificNeighborData(const AcousticMaterial& material,
                                             AcousticNeighborData* localData) {}

  static void getPlaneWaveOperator(const AcousticMaterial& material,
                                   const double n[3],
                                   std::complex<double> mdata[AcousticMaterial::NumQuantities *
                                                              AcousticMaterial::NumQuantities]) {
    getElasticPlaneWaveOperator(material, n, mdata);
  }
  template <typename T>
  static void getTransposedSourceCoefficientTensor(const AcousticMaterial& material,
                                                   T& sourceMatrix) {}

  static void getFaceRotationMatrix(const CoordinateT& normal,
                                    const CoordinateT& tangent1,
                                    const CoordinateT& tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    matT.setZero();
    matTinv.setZero();

    seissol::transformations::tensor1RotationMatrix(normal, tangent1, tangent2, matT, 1, 1);

    seissol::transformations::inverseTensor1RotationMatrix(
        normal, tangent1, tangent2, matTinv, 1, 1);
  }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ACOUSTICSETUP_H_
