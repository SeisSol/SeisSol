// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SETUP_DIRICHLETCONDITION_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SETUP_DIRICHLETCONDITION_H_

#include "Core/Common/Real.h"
#include "Core/Data/DataReader.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::initializer {

/**
 * The frame the affine boundary condition is stated in. Global is the default; face-aligned
 * lets a condition be stated in terms of the face normal, which a condition that mirrors or
 * fixes a direction needs on a boundary that is not axis-aligned.
 */
enum class BoundaryFrame { Global, FaceAligned };

class DirichletCondition {
  public:
  /// The condition that `model` states.
  explicit DirichletCondition(std::unique_ptr<data::DataReader> model);

  DirichletCondition() = default;

  /// Samples the condition at the barycenter of a face of a cell of the configuration `Cfg`, for
  /// the fused simulations `simulations` of the cell, which take it: their constant terms are those
  /// of the condition, while the map is one for all fused simulations of the cell.
  template <typename Cfg>
  [[nodiscard]] BoundaryFrame query(const double* barycenter,
                                    Real<Cfg>* mapTermsData,
                                    Real<Cfg>* constantTermsData,
                                    const std::vector<std::size_t>& simulations) const;

  private:
  std::unique_ptr<data::DataReader> model_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_WAVEPROPAGATION_SETUP_DIRICHLETCONDITION_H_
