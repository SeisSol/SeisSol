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
#include "DynamicRupture/Typedefs.h"
#include "Equations/Setup.h"
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

#include <cstddef>
#include <yateto.h>

namespace seissol::kernel {
/// The kernel that applies the inverse weight of the medium to the flux buffer of a cell, which
/// the generator writes where a cell can be weight-adjusted (WeightAdjustable), and its batched
/// form on the device; named here so that code for every build can name them.
struct weightAdjust;
struct gpu_weightAdjust;
} // namespace seissol::kernel

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

/// Hands the kernel that applies the inverse weight of the medium (kernel::weightAdjust) what the
/// material of the cell says at its sample points: the coefficients of its operator, of which the
/// weight reads the material's own.
template <typename KernelT, typename LocalIntegrationT>
void bindWeightOperands(KernelT& krnl, const LocalIntegrationT& localIntegration) {
  for (std::size_t a = 0; a < StarCoefficientCount; ++a) {
    krnl.nodalCoefficients(a) = localIntegration.materialCoefficients[a];
  }
}

/// Hands a kernel the source term of a cell, where it is formed from what the
/// material says at the sample points rather than from one matrix.
template <typename KernelT, typename LocalIntegrationT>
void bindSourceOperands(KernelT& krnl, const LocalIntegrationT& localIntegration) {
  if constexpr (NodalSource) {
    for (std::size_t a = 0; a < SourceCoefficientCount; ++a) {
      krnl.sourceCoefficients(a) = localIntegration.sourceCoefficients[a];
    }
  }
}

/// The same for what a sample point deviates from the source term a cell
/// carries for itself.
template <typename KernelT, typename LocalIntegrationT>
void bindSourceDeviationOperands(KernelT& krnl, const LocalIntegrationT& localIntegration) {
  if constexpr (NodalSourceDeviation) {
    for (std::size_t a = 0; a < SourceDeviationCount; ++a) {
      krnl.sourceDeviation(a) = localIntegration.sourceDeviation[a];
    }
  }
}

/// Hands a flux kernel the rotation of a face: one for the face, or one per point of it where a
/// face may be curved.
template <typename KernelT>
void bindFaceRotation(KernelT& krnl, const real* rotation) {
  if constexpr (Curvilinear) {
    krnl.TNodes = rotation;
  } else {
    krnl.T = rotation;
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

/// The same as bindWeightOperands for a batch, where the operands are offsets into the cells
/// rather than pointers.
template <typename KernelT>
void bindWeightOperandsBatched(KernelT& krnl, const real** localIntegrationPtrs) {
  SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, materialCoefficients);
  for (std::size_t a = 0; a < StarCoefficientCount; ++a) {
    krnl.nodalCoefficients(a) = localIntegrationPtrs;
    krnl.extraOffset_nodalCoefficients(a) =
        SEISSOL_ARRAY_OFFSET(LocalIntegrationData, materialCoefficients, a);
  }
}

/// The source term of a batch, where it is formed at the sample points.
template <typename KernelT>
void bindSourceOperandsBatched(KernelT& krnl, const real** localIntegrationPtrs) {
  if constexpr (NodalSource) {
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, sourceCoefficients);
    for (std::size_t a = 0; a < SourceCoefficientCount; ++a) {
      krnl.sourceCoefficients(a) = localIntegrationPtrs;
      krnl.extraOffset_sourceCoefficients(a) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, sourceCoefficients, a);
    }
  }
}

/// The same for what a sample point deviates from the source term of its cell.
template <typename KernelT>
void bindSourceDeviationOperandsBatched(KernelT& krnl, const real** localIntegrationPtrs) {
  if constexpr (NodalSourceDeviation) {
    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, sourceDeviation);
    for (std::size_t a = 0; a < SourceDeviationCount; ++a) {
      krnl.sourceDeviation(a) = localIntegrationPtrs;
      krnl.extraOffset_sourceDeviation(a) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, sourceDeviation, a);
    }
  }
}

namespace internal {
/// Where the scalars of one flux coefficient of one face start in the data of
/// a cell, in reals.
template <typename DataT>
constexpr std::size_t fluxCoefficientOffset(std::size_t face, std::size_t coefficient) {
  constexpr std::size_t PerCoefficient = sizeof(DataT::fluxCoefficients[0][0]) / sizeof(real);
  return SEISSOL_ARRAY_OFFSET(DataT, fluxCoefficients, face) + coefficient * PerCoefficient;
}
} // namespace internal

/// What the correctors of the four faces of a cell apply to the jumps of their traces, for a batch
/// of the kernel that takes the four at once (see bindCorrectorOperands): each face the rotation of
/// the face and its scalars, or the matrix of the cell's own side, from the cell's own data.
template <typename KernelT>
void bindFaceCorrectorOperandsBatched(KernelT& krnl, const real** localIntegrationPtrs) {
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if constexpr (NodalFlux) {
      SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, faceRotation);
      SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, fluxCoefficients);
      if constexpr (Curvilinear) {
        krnl.faceTNodes(face) = localIntegrationPtrs;
        krnl.extraOffset_faceTNodes(face) =
            SEISSOL_ARRAY_OFFSET(LocalIntegrationData, faceRotation, face);
      } else {
        krnl.faceT(face) = localIntegrationPtrs;
        krnl.extraOffset_faceT(face) =
            SEISSOL_ARRAY_OFFSET(LocalIntegrationData, faceRotation, face);
      }
      for (std::size_t coefficient = 0; coefficient < FluxCoefficientCount; ++coefficient) {
        krnl.faceFluxCoefficients(face, coefficient) = localIntegrationPtrs;
        krnl.extraOffset_faceFluxCoefficients(face, coefficient) =
            internal::fluxCoefficientOffset<LocalIntegrationData>(face, coefficient);
      }
    } else {
      SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, nApNm1);
      krnl.faceAplusT(face) = localIntegrationPtrs;
      krnl.extraOffset_faceAplusT(face) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, nApNm1, face);
    }
  }
}

/// Hands the lift of a fault face the operator of one side, in whichever of the two shapes the
/// face stores it (see dr::FaultFluxLayout): the rotation of the face together with the scalars
/// at its points, or the matrix the two fold into.
template <typename KernelT>
void bindFaultFluxOperands(KernelT& krnl, const real* fluxSolver) {
  if constexpr (NodalFaultFlux) {
    if constexpr (Curvilinear) {
      krnl.TPoints = fluxSolver + dr::FaultFluxLayout::RotationOffset;
    } else {
      krnl.T = fluxSolver + dr::FaultFluxLayout::RotationOffset;
    }
    for (std::size_t coefficient = 0; coefficient < FaultFluxCoefficientCount; ++coefficient) {
      krnl.faultFluxCoefficients(coefficient) =
          fluxSolver + dr::FaultFluxLayout::coefficientOffset(coefficient);
    }
  } else {
    krnl.fluxSolver = fluxSolver;
  }
}

/// Hands a kernel that takes the values at the points of a fault face into the coordinates of the
/// fault the rotation of the face (DRGodunovData::dataTinvT): one for the face, or one per point
/// where the face may be curved.
template <typename KernelT>
void bindFaultRotation(KernelT& krnl, const real* rotation) {
  if constexpr (Curvilinear) {
    krnl.TinvTPoints = rotation;
  } else {
    krnl.TinvT = rotation;
  }
}

/// The same for a batch, where the operands are offsets into what each face stores.
template <typename KernelT>
void bindFaultFluxOperandsBatched(KernelT& krnl, const real** fluxSolvers) {
  if constexpr (NodalFaultFlux) {
    krnl.T = fluxSolvers;
    krnl.extraOffset_T = dr::FaultFluxLayout::RotationOffset;
    for (std::size_t coefficient = 0; coefficient < FaultFluxCoefficientCount; ++coefficient) {
      krnl.faultFluxCoefficients(coefficient) = fluxSolvers;
      krnl.extraOffset_faultFluxCoefficients(coefficient) =
          dr::FaultFluxLayout::coefficientOffset(coefficient);
    }
  } else {
    krnl.fluxSolver = fluxSolvers;
  }
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_STAROPERANDS_H_
