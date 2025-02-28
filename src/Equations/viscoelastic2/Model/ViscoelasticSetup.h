// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_VISCOELASTICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_VISCOELASTICSETUP_H_

#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include "generated_code/init.h"
#include <Equations/viscoelastic2/Model/Datastructures.h>
#include <Equations/viscoelastic2/Model/IntegrationData.h>
#include <yateto.h>

namespace seissol::model {
template <std::size_t N>
struct MaterialSetup<ViscoElasticMaterialParametrized<N>> {
  using MaterialT = ViscoElasticMaterialParametrized<N>;

  using Matrix99 = Eigen::Matrix<double, 9, 9>;

  template <typename T>
  static void
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

  /*
   * The new implemenation of attenuation (viscoelastic2) is considered standard. This part will be
   * used unless the old attenuation (viscoelastic) implementation is chosen.
   */
  template <typename T>
  static void getTransposedSourceCoefficientTensor(const MaterialT& material, T& E) {
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      const double* theta = material.theta[mech];
      E(0, mech, 0) = theta[0];
      E(1, mech, 0) = theta[1];
      E(2, mech, 0) = theta[1];
      E(0, mech, 1) = theta[1];
      E(1, mech, 1) = theta[0];
      E(2, mech, 1) = theta[1];
      E(0, mech, 2) = theta[1];
      E(1, mech, 2) = theta[1];
      E(2, mech, 2) = theta[0];
      E(3, mech, 3) = theta[2];
      E(4, mech, 4) = theta[2];
      E(5, mech, 5) = theta[2];
    }
  }

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, unsigned dim, T& AT) {
    ::seissol::model::getTransposedCoefficientMatrix(
        dynamic_cast<const ElasticMaterial&>(material), dim, AT);

    getTransposedViscoelasticCoefficientMatrix(1.0, dim, 0, AT);
  }

  static void getTransposedGodunovState(const MaterialT& local,
                                        const MaterialT& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& QgodLocal,
                                        init::QgodNeighbor::view::type& QgodNeighbor) {
    ::seissol::model::getTransposedGodunovState<ElasticMaterial>(
        dynamic_cast<const ElasticMaterial&>(local),
        dynamic_cast<const ElasticMaterial&>(neighbor),
        faceType,
        QgodLocal,
        QgodNeighbor);
  }

  static void getPlaneWaveOperator(
      const MaterialT& material,
      const double n[3],
      std::complex<double> Mdata[MaterialT::NumQuantities * MaterialT::NumQuantities]) {
    yateto::DenseTensorView<2, std::complex<double>> M(
        Mdata, {MaterialT::NumQuantities, MaterialT::NumQuantities});
    M.setZero();

    double data[MaterialT::NumQuantities * MaterialT::NumQuantities];
    yateto::DenseTensorView<2, double> Coeff(data,
                                             {MaterialT::NumQuantities, MaterialT::NumQuantities});

    for (unsigned d = 0; d < 3; ++d) {
      Coeff.setZero();
      getTransposedCoefficientMatrix(material, d, Coeff);
      for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        getTransposedViscoelasticCoefficientMatrix(material.omega[mech], d, mech, Coeff);
      }

      for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
        for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
          M(i, j) += n[d] * Coeff(j, i);
        }
      }
    }
    double Edata[MaterialT::NumQuantities * MaterialT::NumQuantities];
    yateto::DenseTensorView<3, double> E(Edata, tensor::E::Shape);
    E.setZero();
    getTransposedSourceCoefficientTensor(material, E);
    Coeff.setZero();
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      unsigned offset = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
      for (unsigned i = 0; i < tensor::E::Shape[0]; ++i) {
        for (unsigned j = 0; j < tensor::E::Shape[2]; ++j) {
          Coeff(offset + i, j) = E(i, mech, j);
        }
      }
    }

    // E' = diag(-omega_1 I, ..., -omega_L I)
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      unsigned offset = MaterialT::NumElasticQuantities + MaterialT::NumberPerMechanism * mech;
      yateto::DenseTensorView<2, double> ETblock(data + offset + offset * MaterialT::NumQuantities,
                                                 {MaterialT::NumQuantities, 6});
      for (unsigned i = 0; i < MaterialT::NumberPerMechanism; ++i) {
        ETblock(i, i) = -material.omega[mech];
      }
    }

    for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
      for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
        M(i, j) -= std::complex<double>(0.0, Coeff(j, i));
      }
    }
  }

  static void initializeSpecificLocalData(const MaterialT& material,
                                          real timeStepWidth,
                                          ViscoElasticLocalData* localData) {
    auto E = init::E::view::create(localData->E);
    E.setZero();
    getTransposedSourceCoefficientTensor(material, E);

    auto w = init::w::view::create(localData->w);
    auto W = init::W::view::create(localData->W);
    W.setZero();
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      w(mech) = material.omega[mech];
      W(mech, mech) = -material.omega[mech];
    }
  }

  static void initializeSpecificNeighborData(const MaterialT& localMaterial,
                                             ViscoElasticNeighborData* neighborData) {
    // We only need the local omegas
    auto w = init::w::view::create(neighborData->w);
    for (unsigned mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      w(mech) = localMaterial.omega[mech];
    }
  }

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    ::seissol::model::getFaceRotationMatrix<ElasticMaterial>(
        normal, tangent1, tangent2, matT, matTinv);

    seissol::transformations::symmetricTensor2RotationMatrix(normal,
                                                             tangent1,
                                                             tangent2,
                                                             matT,
                                                             MaterialT::NumElasticQuantities,
                                                             MaterialT::NumElasticQuantities);
  }

  static MaterialT getRotatedMaterialCoefficients(real rotationParameters[36],
                                                  MaterialT& material) {
    return material;
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_VISCOELASTICSETUP_H_
