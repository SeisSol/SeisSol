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
      matM(1, 0) = -material.lambda;
      matM(0, 1) = -rhoInv;
      break;

    case 1:
      matM(2, 0) = -material.lambda;
      matM(0, 2) = -rhoInv;
      break;

    case 2:
      matM(3, 0) = -material.lambda;
      matM(0, 3) = -rhoInv;
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

    // Eigenvector matrix for acoustic wave equation
    // The acoustic system has 4 quantities: pressure p and velocities (v1, v2, v3)
    // Coefficient matrices use negative sign convention to match elastic: A = [[0, -λ], [-1/ρ, 0]]
    // This gives eigenvalues ±√(λ/ρ) = ±c (acoustic wave speed)
    Matrix44 matR = Matrix44::Zero();
    
    // Column 0: outgoing acoustic wave (eigenvalue +c) - ALWAYS uses local material
    matR(0, 0) = local.lambda;  // pressure component
    matR(1, 0) = std::sqrt(local.lambda / local.rho);  // normal velocity component
    
    if (faceType != FaceType::FreeSurface) {
      // For internal faces: Column 1 is incoming wave from neighbor
      matR(0, 1) = neighbor.lambda;  // pressure component
      matR(1, 1) = -std::sqrt(neighbor.lambda / neighbor.rho);  // normal velocity (opposite sign)
    } else {
      // For free surface: Column 1 should also use local material (no neighbor)
      matR(0, 1) = local.lambda;  // pressure component  
      matR(1, 1) = -std::sqrt(local.lambda / local.rho);  // normal velocity (opposite sign)
    }
    
    // Columns 2-3: transverse modes (eigenvalue 0, non-propagating)
    // These modes are needed for matrix invertibility but don't propagate
    matR(2, 2) = local.lambda;
    matR(3, 3) = local.lambda;

    if (faceType == FaceType::FreeSurface) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          qGodNeighbor(i, j) = std::numeric_limits<double>::signaling_NaN();
        }
      }
      qGodLocal.setZero();
      
      // Free surface boundary condition for acoustic materials
      // Following the full elastic free surface pattern (Common.h setBlocks function):
      // - Couple normal traction (pressure) to normal velocity: qGodLocal(0, 1) = -matR(1,0)/matR(0,0)
      // - Set ALL velocity diagonals to identity
      qGodLocal(0, 1) = -matR(1, 0) / matR(0, 0);
      qGodLocal(1, 1) = 1.0;  // normal velocity
      qGodLocal(2, 2) = 1.0;  // tangential velocity 1
      qGodLocal(3, 3) = 1.0;  // tangential velocity 2
    } else {
      // Godunov flux: select outgoing characteristics (positive eigenvalue)
      Matrix44 chi = Matrix44::Zero();
      chi(0, 0) = 1.0;  // Select column 0 (outgoing wave)

      const auto godunov = ((matR * chi) * matR.inverse()).eval();

      // Godunov matrices: qGodLocal = I - godunov^T, qGodNeighbor = godunov^T
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

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    matT.setZero();
    matTinv.setZero();

    // Pressure (row 0) is a scalar, doesn't rotate - set to identity
    matT(0, 0) = 1.0;
    matTinv(0, 0) = 1.0;
    
    // Velocity (rows 1-3) is a vector, rotate it
    seissol::transformations::tensor1RotationMatrix(normal, tangent1, tangent2, matT, 1, 1);
    seissol::transformations::inverseTensor1RotationMatrix(
        normal, tangent1, tangent2, matTinv, 1, 1);
  }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ACOUSTICSETUP_H_
