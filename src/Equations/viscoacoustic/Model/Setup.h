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

#include <array>
#include <cstddef>

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
struct MaterialSetup<ViscoAcousticMaterial<N>>
    : public MaterialSetupDefaults<ViscoAcousticMaterial<N>> {
  using MaterialT = ViscoAcousticMaterial<N>;

  /// The flux of the base material alone. How the anelastic blocks are added
  /// on top -- once per mechanism weighted by its relaxation frequency, or
  /// once with the frequency held elsewhere -- is the solver's decision.
  /// The flux is the base material's, and so is its decomposition.
  static constexpr std::size_t NumCoefficients =
      MaterialSetup<AcousticMaterial>::NumCoefficients;
  static constexpr auto CoefficientEntries = MaterialSetup<AcousticMaterial>::CoefficientEntries;

  static std::array<double, NumCoefficients> getCoefficients(const MaterialT& material) {
    return MaterialSetup<AcousticMaterial>::getCoefficients(
        dynamic_cast<const AcousticMaterial&>(material));
  }

  /// The coupling block one mechanism contributes, with its column relative
  /// to that mechanism's block. The weight is the solver's to supply.
  static constexpr std::array<AnelasticCoefficientEntry, 3> AnelasticEntries{{
      {0, 1, 0, -1.0},
      {1, 2, 0, -1.0},
      {2, 3, 0, -1.0},
  }};

  /// The single theta value of one relaxation mechanism.
  static constexpr std::size_t NumSourceCoefficients = 1;

  static std::array<double, NumSourceCoefficients> getSourceCoefficients(const MaterialT& material,
                                                                         std::size_t mech) {
    return {material.theta[mech][0]};
  }

  static constexpr std::array<SourceCoefficientEntry, 1> SourceEntries{{
      {0, 0, 0, 1.0},
  }};

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& matM) {
    MaterialSetup<AcousticMaterial>::getTransposedCoefficientMatrix(
        dynamic_cast<const AcousticMaterial&>(material), dim, matM);
  }

  /**
   * The source entries one relaxation mechanism contributes, as (row, column,
   * value) within its own block. Where the block ends up -- appended to the
   * quantity axis or in a tensor dimension of its own -- is the solver's
   * business, so it is handed a writer rather than a matrix.
   */
  template <typename F>
  static void forEachSourceEntry(const MaterialT& material, std::size_t mech, const F& write) {
    write(0, 0, material.theta[mech][0]);
  }

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
                                                      T& matM) {
    const auto col = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
    switch (dim) {
    case 0:
      matM(1, col) = -omega;
      break;

    case 1:
      matM(2, col) = -omega;
      break;

    case 2:
      matM(3, col) = -omega;
      break;

    default:
      break;
    }
  }
};

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC

#endif

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_SETUP_H_
