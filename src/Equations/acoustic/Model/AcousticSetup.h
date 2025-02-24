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

#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"
#include <Equations/acoustic/Model/Datastructures.h>
#include <Equations/acoustic/Model/IntegrationData.h>

namespace seissol::model {
using Matrix44 = Eigen::Matrix<double, 4, 4>;

template <>
struct MaterialSetup<AcousticMaterial> {
  template <typename T>
  static void getTransposedCoefficientMatrix(const AcousticMaterial& material, unsigned dim, T& M) {
    M.setZero();

    real rhoInv = 1.0 / material.rho;

    switch (dim) {
    case 0:
      M(1, 0) = material.lambda;
      M(0, 1) = rhoInv;
      break;

    case 1:
      M(2, 0) = material.lambda;
      M(0, 2) = rhoInv;
      break;

    case 2:
      M(3, 0) = material.lambda;
      M(0, 3) = rhoInv;
      break;

    default:
      break;
    }
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const AcousticMaterial& local,
                                        const AcousticMaterial& neighbor,
                                        FaceType faceType,
                                        Tloc& QgodLocal,
                                        Tneigh& QgodNeighbor) {
    QgodNeighbor.setZero();

    // Eigenvectors are precomputed
    Matrix44 R = Matrix44::Zero();
    // scale for better condition number of R
    R(0, 0) = std::sqrt(local.lambda * local.rho);
    R(1, 0) = -local.lambda;
    R(0, 1) = std::sqrt(neighbor.lambda * neighbor.rho);
    R(1, 1) = neighbor.lambda;
    R(2, 2) = local.lambda;
    R(3, 3) = local.lambda;

    if (faceType == FaceType::FreeSurface) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          QgodNeighbor(i, j) = std::numeric_limits<double>::signaling_NaN();
        }
      }
      QgodLocal.setZero();
      QgodLocal(0, 1) = -1 * R(1, 0) * 1 / R(0, 0);
      QgodLocal(1, 1) = 1.0;
    } else {
      Matrix44 chi = Matrix44::Zero();
      chi(0, 0) = 1.0;

      const auto godunov = ((R * chi) * R.inverse()).eval();

      // QgodLocal = I - QgodNeighbor
      for (unsigned i = 0; i < godunov.cols(); ++i) {
        for (unsigned j = 0; j < godunov.rows(); ++j) {
          QgodLocal(i, j) = -godunov(j, i);
          QgodNeighbor(i, j) = godunov(j, i);
        }
      }
      for (unsigned idx = 0; idx < 4; ++idx) {
        QgodLocal(idx, idx) += 1.0;
      }
    }
  }

  static AcousticMaterial getRotatedMaterialCoefficients(real rotationParameters[36],
                                                         AcousticMaterial& material) {
    return material;
  }
  static void initializeSpecificLocalData(const AcousticMaterial& material,
                                          real timeStepWidth,
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

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
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
