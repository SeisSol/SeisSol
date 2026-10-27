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
    Derivatives,
    Godunov,
    FluxSolver,
    Ivelocities, // Idofs at the first velocity column (model::MaterialT::VelocityOffset)
    FaceDisplacement,
    NodalStressTensor,
    Pstrains,
    InitialLoad,
    NodalAvgDisplacements,
    T,
    Tinv,
    FSGData,
    Rhos,
    DirichletOffset,
    ZinvExtra,
    IdofsAne,
    DofsAne,
    DofsExt,
    DerivativesAne,
    DerivativesExt,
    Analytical,
    QEtaNodal,
    Integrals,
    // what the faces of a cell lift into: its field, or the extended field of a solver that keeps
    // the memory variables apart; its flux buffer, where the cell is weight-adjusted
    FluxTarget,
    // the flux buffer of a weight-adjusted cell (see WeightAdjustable)
    FluxBuffer,
    // the trace a cell writes on one of its faces (LTS::FaceTraces)
    FaceTrace,
    // what the corrector of a face reads as the trace of the neighbor across it
    FaceTraceNeighbor,
    // what a cell writes on a face for the neighbor across it (LTS::FaceBuffers), or reads there
    // from the neighbor
    FaceBuffer,
    // the operator of the side of a face, for a kernel that takes the side as data
    // (LTS::TraceSideOperators)
    TraceSide,
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
    // the traces of the time derivatives a side writes on the fault face
    FaceDerivatives = 0,
    // where the trace of a side at one point in time is formed, for a kernel that needs it apart
    Trace,
    QInterpolated,
    TinvT,
    Count
  };
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
  Time = 1 << 0,
  Volume = 1 << 1,
  // what a cell writes on its faces: the traces of its time integral, of its time derivatives, and
  // their sum over the steps of a neighbor that takes fewer
  FaceTrace = 1 << 2,
  NeighborFlux = 1 << 3,
  FaceDisplacements = 1 << 4,
  Plasticity = 1 << 5,
  DrSpaceMap = 1 << 6,
  BoundaryConditions = 1 << 7,
  WeightAdjust = 1 << 8,
  // the corrector of a face, from the two traces of the face
  FaceCorrector = 1 << 9,
  Count = 10,
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
  Count = 4,
  None = encodeAny(Count)
};

enum struct FaceId : size_t { Count = 4, Any = AllBits };
enum struct DrFaceRelations : size_t { Count = 8, Any = AllBits };

enum struct ExchangeInfo : size_t {
  Buffers = 1 << 0,
  Derivatives = 1 << 1,
  Count = 2,
  Any = encodeAny(Count)
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
