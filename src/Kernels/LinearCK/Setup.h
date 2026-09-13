// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_SETUP_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_SETUP_H_

#include "GeneratedCode/init.h"
#include "Kernels/LinearCK/Solver.h"
#include "Model/Common.h"
#include "Model/GodunovFaceFlux.h"

#include <cstddef>

namespace seissol::model {

/**
 * The memory variables, where there are any, sit in Q alongside the elastic
 * quantities. The whole system is therefore one operator, and the source term
 * is one matrix over all of it -- which is why the plane wave operator needs
 * nothing beyond the defaults.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::linearck::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::linearck::Solver, MaterialT> {
  static void assembleGodunovFaceFlux(bool rusanov,
                                      double fluxScale,
                                      FaceType faceType,
                                      real* aPlusT,
                                      real* aMinusT,
                                      const real* centralFlux,
                                      const real* rusanovPlus,
                                      const real* rusanovMinus,
                                      const real* qGodLocal,
                                      const real* qGodNeighbor,
                                      const real* rusanovPlusNull,
                                      const real* rusanovMinusNull,
                                      const real* matT,
                                      const real* matTinv,
                                      const real* matATtilde) {
    detail::assembleGodunovFaceFlux(rusanov,
                                    fluxScale,
                                    faceType,
                                    aPlusT,
                                    aMinusT,
                                    centralFlux,
                                    rusanovPlus,
                                    rusanovMinus,
                                    qGodLocal,
                                    qGodNeighbor,
                                    rusanovPlusNull,
                                    rusanovMinusNull,
                                    matT,
                                    matTinv,
                                    matATtilde);
  }

  /// One anelastic block per mechanism, each weighted by its own relaxation
  /// frequency, because the memory variables share the quantity axis.
  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, std::size_t dim, T& matM) {
    MaterialSetup<MaterialT>::getTransposedCoefficientMatrix(material, dim, matM);
    if constexpr (MaterialT::Mechanisms > 0) {
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        MaterialSetup<MaterialT>::getTransposedAnelasticCoefficientMatrix(
            material.omega[mech], dim, mech, matM);
      }
    }
  }

  /// E^T = [E_1^T ... E_L^T] stacked below the elastic quantities, with the
  /// relaxation on the diagonal.
  template <typename T>
  static void getTransposedSourceCoefficientTensor(const MaterialT& material, T& sourceMatrix) {
    sourceMatrix.setZero();
    if constexpr (MaterialT::Mechanisms == 0) {
      MaterialSetup<MaterialT>::getTransposedSourceCoefficientTensor(material, sourceMatrix);
      return;
    } else {
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        const std::size_t offset =
            MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
        MaterialSetup<MaterialT>::forEachSourceEntry(
            material, mech, [&](std::size_t i, std::size_t j, double value) {
              sourceMatrix(offset + i, j) = value;
            });
        for (std::size_t i = 0; i < MaterialT::NumberPerMechanism; ++i) {
          sourceMatrix(offset + i, offset + i) = -material.omega[mech];
        }
      }
    }
  }

  static void initializeSpecificLocalData(const MaterialT& material,
                                          double /*timeStepWidth*/,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
    sourceMatrix.setZero();
    getTransposedSourceCoefficientTensor(material, sourceMatrix);
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_LINEARCK_SETUP_H_
