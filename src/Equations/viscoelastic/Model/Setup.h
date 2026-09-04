// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_SETUP_H_

#include "Equations/viscoelastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include <yateto.h>

namespace seissol::model {

template <std::size_t N>
struct MaterialSetup<ViscoElasticMaterial<N>,
                     std::enable_if_t<ViscoElasticMaterial<N>::ViscoMode ==
                                      ViscoImplementation::QuantityExtension>> {
  using MaterialT = ViscoElasticMaterial<N>;

  template <typename T>
  static void getTransposedViscoelasticCoefficientMatrix(double omega,
                                                         std::size_t dim,
                                                         std::size_t mech,
                                                         T& M) {
    const auto col = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
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
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      const std::size_t offset =
          MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
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
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      for (std::size_t i = 0; i < MaterialT::NumberPerMechanism; ++i) {
        const std::size_t idx =
            MaterialT::NumElasticQuantities + MaterialT::NumberPerMechanism * mech + i;
        sourceMatrix(idx, idx) = -material.omega[mech];
      }
    }
  }

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& AT) {
    getTransposedCoefficientMatrix(dynamic_cast<const ElasticMaterial&>(material), dim, AT);

    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      getTransposedViscoelasticCoefficientMatrix(material.omega[mech], dim, mech, AT);
    }
  }

  static void getTransposedGodunovState(const MaterialT& local,
                                        const MaterialT& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& qGodLocal,
                                        init::QgodNeighbor::view::type& qGodNeighbor) {
    seissol::model::getTransposedGodunovState(dynamic_cast<const ElasticMaterial&>(local),
                                              dynamic_cast<const ElasticMaterial&>(neighbor),
                                              faceType,
                                              qGodLocal,
                                              qGodNeighbor);
  }

  static void initializeSpecificLocalData(const MaterialT& material,
                                          double /*timeStepWidth*/,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
    getTransposedSourceCoefficientTensor(material, sourceMatrix);
  }

  static void initializeSpecificNeighborData(const MaterialT& material,
                                             typename MaterialT::Solver::NeighborData* localData) {}


  static MaterialT
      getRotatedMaterialCoefficients(const std::array<double, 36>& /*rotationParameters*/,
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

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC

template <std::size_t N>
struct MaterialSetup<
    ViscoElasticMaterial<N>,
    std::enable_if_t<ViscoElasticMaterial<N>::ViscoMode == ViscoImplementation::AnelasticTensor>> {
  using MaterialT = ViscoElasticMaterial<N>;

  template <typename T>
  static void getTransposedViscoelasticCoefficientMatrix(double omega,
                                                         std::size_t dim,
                                                         std::size_t mech,
                                                         T& M) {
    const std::size_t col = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
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
  static void getTransposedSourceCoefficientTensor(const MaterialT& material, T& E) {
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
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
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& AT) {
    ::seissol::model::getTransposedCoefficientMatrix(
        dynamic_cast<const ElasticMaterial&>(material), dim, AT);

    getTransposedViscoelasticCoefficientMatrix(1.0, dim, 0, AT);
  }

  static void getTransposedGodunovState(const MaterialT& local,
                                        const MaterialT& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& qGodLocal,
                                        init::QgodNeighbor::view::type& qGodNeighbor) {
    ::seissol::model::getTransposedGodunovState<ElasticMaterial>(
        dynamic_cast<const ElasticMaterial&>(local),
        dynamic_cast<const ElasticMaterial&>(neighbor),
        faceType,
        qGodLocal,
        qGodNeighbor);
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

    for (std::size_t d = 0; d < 3; ++d) {
      Coeff.setZero();
      getTransposedCoefficientMatrix(material, d, Coeff);
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        getTransposedViscoelasticCoefficientMatrix(material.omega[mech], d, mech, Coeff);
      }

      for (std::size_t i = 0; i < MaterialT::NumQuantities; ++i) {
        for (std::size_t j = 0; j < MaterialT::NumQuantities; ++j) {
          M(i, j) += n[d] * Coeff(j, i);
        }
      }
    }
    double Edata[MaterialT::NumQuantities * MaterialT::NumQuantities];
    yateto::DenseTensorView<3, double> E(Edata, tensor::E::Shape);
    E.setZero();
    getTransposedSourceCoefficientTensor(material, E);
    Coeff.setZero();
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      std::size_t offset = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
      for (std::size_t i = 0; i < tensor::E::Shape[0]; ++i) {
        for (std::size_t j = 0; j < tensor::E::Shape[2]; ++j) {
          Coeff(offset + i, j) = E(i, mech, j);
        }
      }
    }

    // E' = diag(-omega_1 I, ..., -omega_L I)
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      std::size_t offset = MaterialT::NumElasticQuantities + MaterialT::NumberPerMechanism * mech;
      yateto::DenseTensorView<2, double> ETblock(
          data + offset + offset * MaterialT::NumQuantities,
          {MaterialT::NumQuantities, MaterialT::NumberPerMechanism});
      for (std::size_t i = 0; i < MaterialT::NumberPerMechanism; ++i) {
        ETblock(i, i) = -material.omega[mech];
      }
    }

    for (std::size_t i = 0; i < MaterialT::NumQuantities; ++i) {
      for (std::size_t j = 0; j < MaterialT::NumQuantities; ++j) {
        M(i, j) -= std::complex<double>(0.0, Coeff(j, i));
      }
    }
  }

  static void initializeSpecificLocalData(const MaterialT& material,
                                          double timeStepWidth,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto E = init::E::view::create(localData->E);
    E.setZero();
    getTransposedSourceCoefficientTensor(material, E);

    auto w = init::w::view::create(localData->w);
    auto W = init::W::view::create(localData->W);
    W.setZero();
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      w(mech) = material.omega[mech];
      W(mech, mech) = -material.omega[mech];
    }
  }

  static void
      initializeSpecificNeighborData(const MaterialT& localMaterial,
                                     typename MaterialT::Solver::NeighborData* neighborData) {
    // We only need the local omegas
    auto w = init::w::view::create(neighborData->w);
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      w(mech) = localMaterial.omega[mech];
    }
  }


  static MaterialT
      getRotatedMaterialCoefficients(const std::array<double, 36>& /*rotationParameters*/,
                                     MaterialT& material) {
    return material;
  }
};

#endif

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_SETUP_H_
