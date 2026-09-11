// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SETUP_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SETUP_H_

#include "GeneratedCode/init.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Model/Common.h"

#include <complex>
#include <cstddef>
#include <yateto.h>

namespace seissol::model {

/**
 * The memory variables live in a tensor dimension of their own, so the star
 * matrix carries a single anelastic block and the relaxation frequencies are
 * held separately. Nothing of that is visible in the material's description of
 * itself: it supplies one coupling block and one source prototype, and what
 * follows is this solver's arithmetic.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::linearckanelastic::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::linearckanelastic::Solver, MaterialT> {
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
      MaterialSetup<MaterialT>::getTransposedCoefficientMatrix(material, d, Coeff);
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        MaterialSetup<MaterialT>::getTransposedAnelasticCoefficientMatrix(
            material.omega[mech], d, mech, Coeff);
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
    MaterialSetup<MaterialT>::getTransposedSourceCoefficientTensor(material, E);
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
    MaterialSetup<MaterialT>::getTransposedSourceCoefficientTensor(material, E);

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

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SETUP_H_
