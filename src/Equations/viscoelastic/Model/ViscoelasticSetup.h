// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_

#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include "generated_code/init.h"
#include <yateto.h>

namespace seissol::model {
template <std::size_t N>
struct MaterialSetup<ViscoElasticMaterialParametrized<N>> {
  using MaterialT = ViscoElasticMaterialParametrized<N>;

  template <typename T>
  static void
      getTransposedViscoelasticCoefficientMatrix(real omega, unsigned dim, unsigned mech, T& M) {
    unsigned col = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
    switch (dim) {
    case 0:
      M(6, col) = -omega;
      M(7, col + 3) = -0.5 * omega;
      M(8, col + 5) = -0.5 * omega;
      break;

    case 1:
      M(7, col + 1) = -omega;
      M(6, col + 3) = -0.5 * omega;
      M(8, col + 4) = -0.5 * omega;
      break;

    case 2:
      M(8, col + 2) = -omega;
      M(7, col + 4) = -0.5 * omega;
      M(6, col + 5) = -0.5 * omega;
      break;
    }
  }

  template <typename T>
  static void getTransposedSourceCoefficientTensor(const MaterialT& material, T& sourceMatrix) {
    sourceMatrix.setZero();

    //       | E_1^T |
    // E^T = |  ...  |
    //       | E_L^T |
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      unsigned offset = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
      const double* theta = material.theta[mech];
      sourceMatrix(offset, 0) = theta[0];
      sourceMatrix(offset + 1, 0) = theta[1];
      sourceMatrix(offset + 2, 0) = theta[1];
      sourceMatrix(offset, 1) = theta[1];
      sourceMatrix(offset + 1, 1) = theta[0];
      sourceMatrix(offset + 2, 1) = theta[1];
      sourceMatrix(offset, 2) = theta[1];
      sourceMatrix(offset + 1, 2) = theta[1];
      sourceMatrix(offset + 2, 2) = theta[0];
      sourceMatrix(offset + 3, 3) = theta[2];
      sourceMatrix(offset + 4, 4) = theta[2];
      sourceMatrix(offset + 5, 5) = theta[2];
    }

    // E' = diag(-omega_1 I, ..., -omega_L I)
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      for (unsigned i = 0; i < MaterialT::NumberPerMechanism; ++i) {
        unsigned idx = MaterialT::NumElasticQuantities + MaterialT::NumberPerMechanism * mech + i;
        sourceMatrix(idx, idx) = -material.omega[mech];
      }
    }
  }

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, unsigned dim, T& AT) {
    getTransposedCoefficientMatrix(dynamic_cast<const ElasticMaterial&>(material), dim, AT);

    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      getTransposedViscoelasticCoefficientMatrix(material.omega[mech], dim, mech, AT);
    }
  }

  static void getTransposedGodunovState(const MaterialT& local,
                                        const MaterialT& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& QgodLocal,
                                        init::QgodNeighbor::view::type& QgodNeighbor) {
    seissol::model::getTransposedGodunovState(dynamic_cast<const ElasticMaterial&>(local),
                                              dynamic_cast<const ElasticMaterial&>(neighbor),
                                              faceType,
                                              QgodLocal,
                                              QgodNeighbor);
  }

  static void initializeSpecificLocalData(const MaterialT& material,
                                          real timeStepWidth,
                                          ViscoElasticLocalData* localData) {
    auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
    getTransposedSourceCoefficientTensor(material, sourceMatrix);
  }

  static void initializeSpecificNeighborData(const MaterialT& material,
                                             ViscoElasticLocalData* localData) {}

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    seissol::model::getFaceRotationMatrix<ElasticMaterial>(
        normal, tangent1, tangent2, matT, matTinv);

    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      const unsigned origin =
          MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
      seissol::transformations::symmetricTensor2RotationMatrix(
          normal, tangent1, tangent2, matT, origin, origin);
      seissol::transformations::inverseSymmetricTensor2RotationMatrix(
          normal, tangent1, tangent2, matTinv, origin, origin);
    }
  }

  static MaterialT getRotatedMaterialCoefficients(real rotationParameters[36],
                                                  MaterialT& material) {
    return material;
  }

  static void getPlaneWaveOperator(
      const MaterialT& material,
      const double n[3],
      std::complex<double> mdata[MaterialT::NumQuantities * MaterialT::NumQuantities]) {
    getElasticPlaneWaveOperator(material, n, mdata);
  }
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_
