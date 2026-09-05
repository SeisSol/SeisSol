// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_SETUP_H_

#include "Equations/viscoacoustic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include <yateto.h>

namespace seissol::model {

/**
 * The parts of the setup that do not depend on how the solver lays out the
 * memory variables: the Riemann problem is the base material's, and the
 * anelastic coupling block is the same matrix either way. Only its
 * coefficient and how often it is written differ, and that is the solver's
 * business.
 */
template <std::size_t N>
struct ViscoAcousticSetupCommon : public MaterialSetupDefaults<ViscoAcousticMaterial<N>> {
  using MaterialT = ViscoAcousticMaterial<N>;

  static void getTransposedGodunovState(const MaterialT& local,
                                        const MaterialT& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& qGodLocal,
                                        init::QgodNeighbor::view::type& qGodNeighbor) {
    seissol::model::getTransposedGodunovState(dynamic_cast<const AcousticMaterial&>(local),
                                              dynamic_cast<const AcousticMaterial&>(neighbor),
                                              faceType,
                                              qGodLocal,
                                              qGodNeighbor);
  }
  template <typename T>
  static void getTransposedViscoacousticCoefficientMatrix(double omega,
                                                          std::size_t dim,
                                                          std::size_t mech,
                                                          T& M) {
    const auto col = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
    switch (dim) {
    case 0:
      M(1, col) = -omega;
      break;

    case 1:
      M(2, col) = -omega;
      break;

    case 2:
      M(3, col) = -omega;
      break;
    }
  }
};

template <std::size_t N>
struct MaterialSetup<
    ViscoAcousticMaterial<N>,
    std::enable_if_t<ViscoAcousticMaterial<N>::ViscoMode == ViscoImplementation::QuantityExtension>>
    : public ViscoAcousticSetupCommon<N> {
  using MaterialT = ViscoAcousticMaterial<N>;
  using ViscoAcousticSetupCommon<N>::getTransposedViscoacousticCoefficientMatrix;

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
    ::seissol::model::getTransposedCoefficientMatrix(
        dynamic_cast<const AcousticMaterial&>(material), dim, AT);

    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      getTransposedViscoacousticCoefficientMatrix(material.omega[mech], dim, mech, AT);
    }
  }

  static void initializeSpecificLocalData(const MaterialT& material,
                                          double /*timeStepWidth*/,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
    getTransposedSourceCoefficientTensor(material, sourceMatrix);
  }
};

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC

template <std::size_t N>
struct MaterialSetup<
    ViscoAcousticMaterial<N>,
    std::enable_if_t<ViscoAcousticMaterial<N>::ViscoMode == ViscoImplementation::AnelasticTensor>>
    : public ViscoAcousticSetupCommon<N> {
  using MaterialT = ViscoAcousticMaterial<N>;
  using ViscoAcousticSetupCommon<N>::getTransposedViscoacousticCoefficientMatrix;

  template <typename T>
  static void getTransposedSourceCoefficientTensor(const MaterialT& material, T& E) {
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      const double* theta = material.theta[mech];
      E(0, mech, 0) = theta[0];
    }
  }

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& AT) {
    ::seissol::model::getTransposedCoefficientMatrix(
        dynamic_cast<const AcousticMaterial&>(material), dim, AT);

    getTransposedViscoacousticCoefficientMatrix(1.0, dim, 0, AT);
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
        getTransposedViscoacousticCoefficientMatrix(material.omega[mech], d, mech, Coeff);
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
};

#endif

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_SETUP_H_
