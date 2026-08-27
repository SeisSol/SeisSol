// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_

#include "Common/Literals.h"
#include "Kernels/Precision.h"

#include <cstdlib>

namespace seissol::recording::inner_keys {

/**
 * The structure contains encoded variables names
 * of the Wave Propagation (Wp) solver.
 */
struct Wp {
  using DataType = real*;
  enum struct Id : size_t {
    Dofs = 0,
    Idofs,
    LocalIntegrationData,
    NeighborIntegrationData,
    Buffers,
    Derivatives,
    Godunov,
    FluxSolver,
    Ivelocities, // 6th, 7the and 8th columns of Idofs
    FaceDisplacement,
    NodalStressTensor,
    Pstrains,
    InitialLoad,
    NodalAvgDisplacements,
    T,
    Tinv,
    EasiBoundaryMap,
    EasiBoundaryConstant,
    ZinvExtra,
    IdofsAne,
    DofsAne,
    DofsExt,
    DerivativesAne,
    DerivativesExt,
    Analytical,
    RotateDisplacementToFaceNormal,
    RotateDisplacementToGlobal,
    RotatedFaceDisplacement,
    DofsFaceNodal,
    PrevCoefficients,
    DofsFaceBoundaryNodal,
    Integrals,
    Count
  };
};

/**
 * The structure contains encoded variables names
 * of the Dynamic Rupture (Dr) solver.
 */
struct Dr {
  using DataType = real*;
  enum struct Id : size_t {
    DerivativesPlus = 0,
    DerivativesMinus,
    IdofsPlus,
    IdofsMinus,
    QInterpolatedPlus,
    QInterpolatedMinus,
    TinvT,
    Count
  };
};

struct Material {
  using DataType = double;
  enum struct Id : size_t { Rho = 0, Lambda, InvImpedances, Count };
};

struct Indices {
  using DataType = unsigned;
  enum struct Id : size_t { Cells = 0, Count };
};
} // namespace seissol::recording::inner_keys

namespace seissol::recording {
constexpr size_t AllBits = ~0_UZ;
constexpr size_t encodeAny(unsigned count) { return ~(AllBits << count); }

enum struct KernelNames : size_t {
  Time = 1_UZ << 0,
  Volume = 1_UZ << 1,
  LocalFlux = 1_UZ << 2,
  NeighborFlux = 1_UZ << 3,
  FaceDisplacements = 1_UZ << 4,
  Plasticity = 1_UZ << 5,
  DrSpaceMap = 1_UZ << 6,
  BoundaryConditions = 1_UZ << 7,
  Count = 9_UZ,
  Any = encodeAny(Count),
};

enum struct ComputationKind : size_t {
  WithoutDerivatives = 1_UZ << 0,
  WithDerivatives = 1_UZ << 1,
  WithLtsDerivatives = 1_UZ << 2,
  WithGtsDerivatives = 1_UZ << 3,
  WithGtsBuffers = 1_UZ << 4,
  WithLtsBuffers = 1_UZ << 5,
  FreeSurfaceGravity = 1_UZ << 6,
  Dirichlet = 1_UZ << 7,
  Analytical = 1_UZ << 8,
  Count = 9_UZ,
  None = encodeAny(Count),
};

enum struct FaceKinds : size_t {
  Regular = 1_UZ << 0,
  FreeSurface = 1_UZ << 1,
  Outflow = 1_UZ << 2,
  DynamicRupture = 1_UZ << 3,
  Count = 4_UZ,
  None = encodeAny(Count),
};

enum struct FaceId : size_t {
  Count = 4,
  Any = AllBits,
};
enum struct FaceRelations : size_t {
  Count = 48,
  Any = AllBits,
};
enum struct DrFaceRelations : size_t {
  Count = 16,
  Any = AllBits,
};

enum struct ExchangeInfo : size_t {
  Buffers = 1_UZ << 0,
  Derivatives = 1_UZ << 1,
  Count = 2_UZ,
  Any = encodeAny(Count),
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
