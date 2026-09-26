// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_STAROPERANDS_H_
#define SEISSOL_SRC_KERNELS_STAROPERANDS_H_

#include "Common/Constants.h"
#include "Common/Offset.h"
#include "Equations/Setup.h"
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

#include <cstddef>
#include <yateto.h>

namespace seissol::kernels {

static_assert(!FactoredStar ||
                  StarCoefficientCount == model::SolverSetup<typename model::MaterialT::Solver,
                                                             model::MaterialT>::NumCoefficients,
              "the generated coefficient count and the solver's declaration disagree");

namespace internal {
/// The generator works out the origins from the same composition rule the
/// solver's declaration follows, so the two have to agree entry for entry.
constexpr bool originsAgree() {
  // a build that does not factor the star carries no coefficients, and the
  // generator emits an empty table for it
  if (generated::SolverCoefficientOrigins.empty()) {
    return true;
  }
  constexpr auto declared =
      model::SolverSetup<typename model::MaterialT::Solver, model::MaterialT>::CoefficientOrigins;
  if (declared.size() != generated::SolverCoefficientOrigins.size()) {
    return false;
  }
  for (std::size_t i = 0; i < declared.size(); ++i) {
    if (declared[i] != generated::SolverCoefficientOrigins[i]) {
      return false;
    }
  }
  return true;
}
} // namespace internal

static_assert(internal::originsAgree(),
              "the generated coefficient origins and the solver's declaration disagree");

/// Hands a kernel the operator of a cell, in whichever of the two shapes the
/// cell carries it.
template <typename KernelT, typename LocalIntegrationT>
void bindStarOperands(KernelT& krnl, const LocalIntegrationT& localIntegration) {
  if constexpr (NodalMaterial) {
    for (std::size_t a = 0; a < StarCoefficientCount; ++a) {
      krnl.nodalCoefficients(a) = localIntegration.materialCoefficients[a];
    }
    for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
      krnl.referenceGradients(dim) = localIntegration.referenceGradients[dim];
    }
  } else if constexpr (FactoredStar) {
    krnl.materialCoefficients = localIntegration.materialCoefficients[0];
    for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
      krnl.referenceGradients(dim) = localIntegration.referenceGradients[dim];
    }
  } else {
    for (std::size_t dim = 0; dim < yateto::numFamilyMembers<tensor::star>(); ++dim) {
      krnl.star(dim) = localIntegration.starMatrices[dim];
    }
  }
}

/// The same for a batch, where the operands are offsets into the cells rather
/// than pointers.
template <typename KernelT>
void bindStarOperandsBatched(KernelT& krnl, const real** localIntegrationPtrs) {
  if constexpr (NodalMaterial) {
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, materialCoefficients);
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, referenceGradients);
    for (std::size_t a = 0; a < StarCoefficientCount; ++a) {
      krnl.nodalCoefficients(a) = localIntegrationPtrs;
      krnl.extraOffset_nodalCoefficients(a) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, materialCoefficients, a);
    }
    for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
      krnl.referenceGradients(dim) = localIntegrationPtrs;
      krnl.extraOffset_referenceGradients(dim) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, referenceGradients, dim);
    }
  } else if constexpr (FactoredStar) {
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, materialCoefficients);
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, referenceGradients);
    krnl.materialCoefficients = localIntegrationPtrs;
    krnl.extraOffset_materialCoefficients =
        SEISSOL_ARRAY_OFFSET(LocalIntegrationData, materialCoefficients, 0);
    for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
      krnl.referenceGradients(dim) = localIntegrationPtrs;
      krnl.extraOffset_referenceGradients(dim) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, referenceGradients, dim);
    }
  } else {
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, starMatrices);
    for (std::size_t dim = 0; dim < yateto::numFamilyMembers<tensor::star>(); ++dim) {
      krnl.star(dim) = localIntegrationPtrs;
      krnl.extraOffset_star(dim) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, dim);
    }
  }
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_STAROPERANDS_H_
