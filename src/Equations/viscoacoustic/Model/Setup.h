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
  static void getTransposedAnelasticCoefficientMatrix(double omega,
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
  using ViscoAcousticSetupCommon<N>::getTransposedAnelasticCoefficientMatrix;

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
      getTransposedAnelasticCoefficientMatrix(material.omega[mech], dim, mech, AT);
    }
  }
};

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC

template <std::size_t N>
struct MaterialSetup<
    ViscoAcousticMaterial<N>,
    std::enable_if_t<ViscoAcousticMaterial<N>::ViscoMode == ViscoImplementation::AnelasticTensor>>
    : public ViscoAcousticSetupCommon<N> {
  using MaterialT = ViscoAcousticMaterial<N>;
  using ViscoAcousticSetupCommon<N>::getTransposedAnelasticCoefficientMatrix;

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

    getTransposedAnelasticCoefficientMatrix(1.0, dim, 0, AT);
  }
};

#endif

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_SETUP_H_
