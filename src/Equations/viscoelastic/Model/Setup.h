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

/**
 * The parts of the setup that do not depend on how the solver lays out the
 * memory variables: the Riemann problem is the base material's, and the
 * anelastic coupling block is the same matrix either way. Only its
 * coefficient and how often it is written differ, and that is the solver's
 * business.
 */
template <std::size_t N>
struct ViscoElasticSetupCommon : public MaterialSetupDefaults<ViscoElasticMaterial<N>> {
  using MaterialT = ViscoElasticMaterial<N>;

  /**
   * The source entries one relaxation mechanism contributes, as (row, column,
   * value) within its own block. Where the block ends up -- appended to the
   * quantity axis or in a tensor dimension of its own -- is the solver's
   * business, so it is handed a writer rather than a matrix.
   */
  template <typename F>
  static void forEachSourceEntry(const MaterialT& material, std::size_t mech, F&& write) {
    const double* theta = material.theta[mech];
    write(0, 0, theta[0]);
    write(1, 0, theta[1]);
    write(2, 0, theta[1]);
    write(0, 1, theta[1]);
    write(1, 1, theta[0]);
    write(2, 1, theta[1]);
    write(0, 2, theta[1]);
    write(1, 2, theta[1]);
    write(2, 2, theta[0]);
    write(3, 3, theta[2]);
    write(4, 4, theta[2]);
    write(5, 5, theta[2]);
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
  template <typename T>
  static void getTransposedAnelasticCoefficientMatrix(double omega,
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
};

template <std::size_t N>
struct MaterialSetup<
    ViscoElasticMaterial<N>,
    std::enable_if_t<ViscoElasticMaterial<N>::ViscoMode == ViscoImplementation::QuantityExtension>>
    : public ViscoElasticSetupCommon<N> {
  using MaterialT = ViscoElasticMaterial<N>;
  using ViscoElasticSetupCommon<N>::getTransposedAnelasticCoefficientMatrix;

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& AT) {
    ::seissol::model::getTransposedCoefficientMatrix(
        dynamic_cast<const ElasticMaterial&>(material), dim, AT);

    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      getTransposedAnelasticCoefficientMatrix(material.omega[mech], dim, mech, AT);
    }
  }
};

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC

template <std::size_t N>
struct MaterialSetup<
    ViscoElasticMaterial<N>,
    std::enable_if_t<ViscoElasticMaterial<N>::ViscoMode == ViscoImplementation::AnelasticTensor>>
    : public ViscoElasticSetupCommon<N> {
  using MaterialT = ViscoElasticMaterial<N>;
  using ViscoElasticSetupCommon<N>::getTransposedAnelasticCoefficientMatrix;

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& AT) {
    ::seissol::model::getTransposedCoefficientMatrix(
        dynamic_cast<const ElasticMaterial&>(material), dim, AT);

    getTransposedAnelasticCoefficientMatrix(1.0, dim, 0, AT);
  }
};

#endif

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_SETUP_H_
