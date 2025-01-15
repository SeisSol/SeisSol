// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
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
template <typename T>
inline void
    getTransposedViscoelasticCoefficientMatrix(real omega, unsigned dim, unsigned mech, T& M) {
  unsigned col = 9 + mech * 6;
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
inline void getTransposedSourceCoefficientTensor(const ViscoElasticMaterial& material,
                                                 T& sourceMatrix) {
  sourceMatrix.setZero();

  //       | E_1^T |
  // E^T = |  ...  |
  //       | E_L^T |
  for (unsigned mech = 0; mech < seissol::model::MaterialT::Mechanisms; ++mech) {
    unsigned offset = seissol::model::MaterialT::NumElasticQuantities +
                      mech * seissol::model::MaterialT::NumberPerMechanism;
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
  for (unsigned mech = 0; mech < seissol::model::MaterialT::Mechanisms; ++mech) {
    for (unsigned i = 0; i < 6; ++i) {
      unsigned idx = seissol::model::MaterialT::NumElasticQuantities +
                     seissol::model::MaterialT::NumberPerMechanism * mech + i;
      sourceMatrix(idx, idx) = -material.omega[mech];
    }
  }
}

template <typename T>
inline void
    getTransposedCoefficientMatrix(const ViscoElasticMaterial& material, unsigned dim, T& AT) {
  getTransposedCoefficientMatrix(dynamic_cast<const ElasticMaterial&>(material), dim, AT);

  for (unsigned mech = 0; mech < seissol::model::MaterialT::Mechanisms; ++mech) {
    getTransposedViscoelasticCoefficientMatrix(material.omega[mech], dim, mech, AT);
  }
}

template <>
inline void getTransposedGodunovState(const ViscoElasticMaterial& local,
                                      const ViscoElasticMaterial& neighbor,
                                      FaceType faceType,
                                      init::QgodLocal::view::type& QgodLocal,
                                      init::QgodNeighbor::view::type& QgodNeighbor) {
  getTransposedGodunovState(dynamic_cast<const ElasticMaterial&>(local),
                            dynamic_cast<const ElasticMaterial&>(neighbor),
                            faceType,
                            QgodLocal,
                            QgodNeighbor);
}

template <>
inline void initializeSpecificLocalData(const ViscoElasticMaterial& material,
                                        real timeStepWidth,
                                        ViscoElasticLocalData* localData) {
  auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
  getTransposedSourceCoefficientTensor(material, sourceMatrix);
}
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_
