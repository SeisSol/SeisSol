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

#include <array>
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
  /// The material's coefficients, then one relaxation frequency per
  /// mechanism, because each block carries its own weight here.
  static constexpr std::size_t NumCoefficients =
      MaterialSetup<MaterialT>::NumCoefficients + MaterialT::Mechanisms;

  /// The material's are fields; the relaxation frequencies are not. No
  /// material parameter moves them -- they follow the frequency band alone --
  /// so a cell does not have to carry them.
  static constexpr std::array<CoefficientOrigin, NumCoefficients> CoefficientOrigins = [] {
    std::array<CoefficientOrigin, NumCoefficients> origins{};
    const auto base = materialCoefficientOrigins<MaterialT>();
    for (std::size_t i = 0; i < base.size(); ++i) {
      origins[i] = base[i];
    }
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      origins[base.size() + mech] = CoefficientOrigin::Global;
    }
    return origins;
  }();

  static std::array<double, NumCoefficients> getCoefficients(const MaterialT& material) {
    std::array<double, NumCoefficients> coefficients{};
    const auto base = MaterialSetup<MaterialT>::getCoefficients(material);
    for (std::size_t i = 0; i < base.size(); ++i) {
      coefficients[i] = base[i];
    }
    if constexpr (MaterialT::Mechanisms > 0) {
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        coefficients[base.size() + mech] = material.omega[mech];
      }
    }
    return coefficients;
  }

  template <typename F>
  static void forEachCoefficientEntry(const F& write) {
    for (const auto& entry : MaterialSetup<MaterialT>::CoefficientEntries) {
      write(entry.coefficient, entry.dim, entry.row, entry.column, entry.factor);
    }
    if constexpr (MaterialT::Mechanisms > 0) {
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        const auto col = MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
        for (const auto& entry : MaterialSetup<MaterialT>::AnelasticEntries) {
          write(MaterialSetup<MaterialT>::NumCoefficients + mech,
                entry.dim,
                entry.row,
                col + entry.columnOffset,
                entry.factor);
        }
      }
    }
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
  /// The source of every mechanism, each with its own relaxation frequency on
  /// the diagonal: the material's scalars per block, then that frequency.
  static constexpr std::size_t SourcePerMechanism =
      MaterialSetup<MaterialT>::NumSourceCoefficients + 1;
  static constexpr std::size_t NumSourceCoefficients =
      MaterialT::Mechanisms > 0 ? SourcePerMechanism* MaterialT::Mechanisms
                                : MaterialSetup<MaterialT>::NumSourceCoefficients;

  static std::array<double, NumSourceCoefficients>
      getSourceCoefficients(const MaterialT& material) {
    std::array<double, NumSourceCoefficients> coefficients{};
    if constexpr (MaterialT::Mechanisms == 0) {
      if constexpr (NumSourceCoefficients > 0) {
        coefficients = MaterialSetup<MaterialT>::getSourceCoefficients(material, 0);
      }
    } else {
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        const auto block = MaterialSetup<MaterialT>::getSourceCoefficients(material, mech);
        for (std::size_t i = 0; i < block.size(); ++i) {
          coefficients[mech * SourcePerMechanism + i] = block[i];
        }
        coefficients[mech * SourcePerMechanism + block.size()] = material.omega[mech];
      }
    }
    return coefficients;
  }

  template <typename F>
  static void forEachSourceCoefficientEntry(const F& write) {
    if constexpr (MaterialT::Mechanisms == 0) {
      SolverSetupDefaults<kernels::solver::linearck::Solver,
                          MaterialT>::forEachSourceCoefficientEntry(write);
    } else {
      constexpr std::size_t PerBlock = MaterialSetup<MaterialT>::NumSourceCoefficients;
      for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
        const std::size_t offset =
            MaterialT::NumElasticQuantities + mech * MaterialT::NumberPerMechanism;
        for (const auto& entry : MaterialSetup<MaterialT>::SourceEntries) {
          write(mech * SourcePerMechanism + entry.coefficient,
                offset + entry.row,
                entry.column,
                entry.factor);
        }
        for (std::size_t i = 0; i < MaterialT::NumberPerMechanism; ++i) {
          write(mech * SourcePerMechanism + PerBlock, offset + i, offset + i, -1.0);
        }
      }
    }
  }

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
