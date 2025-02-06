// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_

#include "Kernels/Precision.h"
#include <cstdlib>

namespace seissol::initializer::recording::inner_keys {

/**
 * The structure contains encoded variables names
 * of the Wave Propagation (Wp) solver.
 */
struct Wp {
  using DataType = real*;
  enum struct Id : size_t {
    Dofs = 0,
    Idofs,
    Star,
    Buffers,
    Derivatives,
    AplusT,
    AminusT,
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
    Analytical,
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
  enum struct Id : size_t { Rho = 0, Lambda, Count };
};

struct Indices {
  using DataType = unsigned;
  enum struct Id : size_t { Cells = 0, Count };
};
} // namespace seissol::initializer::recording::inner_keys

namespace seissol::initializer::recording {
constexpr size_t ALL_BITS = ~static_cast<size_t>(0);
constexpr size_t encodeAny(unsigned count) { return ~(ALL_BITS << count); }

enum struct KernelNames : size_t {
  Time = 1 << 0,
  Volume = 1 << 1,
  LocalFlux = 1 << 2,
  NeighborFlux = 1 << 3,
  FaceDisplacements = 1 << 4,
  Plasticity = 1 << 5,
  DrTime = 1 << 6,
  DrSpaceMap = 1 << 7,
  BoundaryConditions = 1 << 8,
  Count = 9,
  Any = encodeAny(Count)
};

enum struct ComputationKind : size_t {
  WithoutDerivatives = 1 << 0,
  WithDerivatives = 1 << 1,
  WithLtsDerivatives = 1 << 2,
  WithGtsDerivatives = 1 << 3,
  WithGtsBuffers = 1 << 4,
  WithLtsBuffers = 1 << 5,
  FreeSurfaceGravity = 1 << 6,
  Dirichlet = 1 << 7,
  Analytical = 1 << 8,
  Count = 9,
  None = encodeAny(Count)
};

enum struct FaceKinds : size_t {
  Regular = 1 << 0,
  FreeSurface = 1 << 1,
  Outflow = 1 << 2,
  DynamicRupture = 1 << 3,
  Periodic = 1 << 4,
  Count = 5,
  None = encodeAny(Count)
};

enum struct FaceId : size_t { Count = 4, Any = ALL_BITS };
enum struct FaceRelations : size_t { Count = 48, Any = ALL_BITS };
enum struct DrFaceRelations : size_t { Count = 16, Any = ALL_BITS };

enum struct ExchangeInfo : size_t {
  Buffers = 1 << 0,
  Derivatives = 1 << 1,
  Count = 2,
  Any = encodeAny(Count)
};

} // namespace seissol::initializer::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
