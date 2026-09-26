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
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Equations/Setup.h"
#include "Kernels/Precision.h"

#include <cstddef>
#include <yateto.h>

namespace seissol::kernels {

static_assert(!FactoredStar || StarCoefficientCount ==
                                   model::SolverSetup<typename model::MaterialT::Solver,
                                                      model::MaterialT>::NumCoefficients,
              "the generated coefficient count and the solver's declaration disagree");

/// Hands a kernel the operator of a cell, in whichever of the two shapes the
/// cell carries it.
template <typename KernelT, typename LocalIntegrationT>
void bindStarOperands(KernelT& krnl, const LocalIntegrationT& localIntegration) {
  if constexpr (FactoredStar) {
    krnl.materialCoefficients = localIntegration.materialCoefficients;
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
  if constexpr (FactoredStar) {
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
      krnl.extraOffset_star(dim) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, dim);
    }
  }
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_STAROPERANDS_H_
