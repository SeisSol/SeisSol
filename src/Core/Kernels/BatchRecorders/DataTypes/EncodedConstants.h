// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
#define SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_

#include "Core/Common/Literals.h"
#include "Core/Common/TypeId.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

namespace seissol::recording {

/**
 * A variable of the recorded batches (see BatchTable), named by the type `Tag`. A batch takes any
 * variable, so that every part of the program can name its own, e.g.
 *
 *   struct MyVariableTag {};
 *   inline constexpr BatchVariable<MyVariableTag> MyVariable{};
 */
template <typename Tag>
struct BatchVariable {
  /// The index of the variable in the tables of the batches (see TypeId).
  [[nodiscard]] static std::uint32_t id() { return TypeId<TypeDomain::Batch, Tag>; }
};

/**
 * A kernel that batches are recorded for, named by the type `Tag`: a bit of the keys of its batches
 * (see ConditionalKey); kernels that run together combine their bits (see operator||). As with the
 * variables, every part of the program can name its own kernels, up to 64 in all.
 */
template <typename Tag>
struct KernelName {
  [[nodiscard]] static std::size_t bit() {
    const auto id = TypeId<TypeDomain::Kernel, Tag>;
    assert(id < 64 && "too many kernel names for the bits of a key");
    return std::size_t{1} << id;
  }

  // (hidden friends, so that they do not hide the operators of the encoded constants)

  /// The bit of the kernel, for the key of a batch.
  friend std::size_t operator*(KernelName /*kernel*/) { return bit(); }

  /// The bits of two kernels that run together, for the key of their batch.
  template <typename OtherTag>
  friend std::size_t operator||(KernelName /*first*/, KernelName<OtherTag> /*second*/) {
    return bit() | KernelName<OtherTag>::bit();
  }
};

} // namespace seissol::recording

namespace seissol::recording::inner_keys {

/**
 * The variables of the batches of the Wave Propagation (Wp) solver.
 */
struct Wp {
  struct Id {
    struct DofsTag {};
    static constexpr BatchVariable<DofsTag> Dofs{};
    struct IdofsTag {};
    static constexpr BatchVariable<IdofsTag> Idofs{};
    struct LocalIntegrationDataTag {};
    static constexpr BatchVariable<LocalIntegrationDataTag> LocalIntegrationData{};
    struct NeighborIntegrationDataTag {};
    static constexpr BatchVariable<NeighborIntegrationDataTag> NeighborIntegrationData{};
    struct BuffersTag {};
    static constexpr BatchVariable<BuffersTag> Buffers{};
    struct DerivativesTag {};
    static constexpr BatchVariable<DerivativesTag> Derivatives{};
    struct GodunovTag {};
    static constexpr BatchVariable<GodunovTag> Godunov{};
    struct FluxSolverTag {};
    static constexpr BatchVariable<FluxSolverTag> FluxSolver{};
    // Idofs at the first velocity column (the VelocityOffset of the material)
    struct IvelocitiesTag {};
    static constexpr BatchVariable<IvelocitiesTag> Ivelocities{};
    struct FaceDisplacementTag {};
    static constexpr BatchVariable<FaceDisplacementTag> FaceDisplacement{};
    struct NodalStressTensorTag {};
    static constexpr BatchVariable<NodalStressTensorTag> NodalStressTensor{};
    struct PstrainsTag {};
    static constexpr BatchVariable<PstrainsTag> Pstrains{};
    struct InitialLoadTag {};
    static constexpr BatchVariable<InitialLoadTag> InitialLoad{};
    struct NodalAvgDisplacementsTag {};
    static constexpr BatchVariable<NodalAvgDisplacementsTag> NodalAvgDisplacements{};
    struct TTag {};
    static constexpr BatchVariable<TTag> T{};
    struct TinvTag {};
    static constexpr BatchVariable<TinvTag> Tinv{};
    struct FSGDataTag {};
    static constexpr BatchVariable<FSGDataTag> FSGData{};
    struct RhosTag {};
    static constexpr BatchVariable<RhosTag> Rhos{};
    struct DirichletOffsetTag {};
    static constexpr BatchVariable<DirichletOffsetTag> DirichletOffset{};
    struct ZinvExtraTag {};
    static constexpr BatchVariable<ZinvExtraTag> ZinvExtra{};
    struct IdofsAneTag {};
    static constexpr BatchVariable<IdofsAneTag> IdofsAne{};
    struct DofsAneTag {};
    static constexpr BatchVariable<DofsAneTag> DofsAne{};
    struct DofsExtTag {};
    static constexpr BatchVariable<DofsExtTag> DofsExt{};
    struct DerivativesAneTag {};
    static constexpr BatchVariable<DerivativesAneTag> DerivativesAne{};
    struct DerivativesExtTag {};
    static constexpr BatchVariable<DerivativesExtTag> DerivativesExt{};
    struct AnalyticalTag {};
    static constexpr BatchVariable<AnalyticalTag> Analytical{};
    struct QEtaNodalTag {};
    static constexpr BatchVariable<QEtaNodalTag> QEtaNodal{};
    struct IntegralsTag {};
    static constexpr BatchVariable<IntegralsTag> Integrals{};
    // a time integral of a neighbor of another configuration in the canonical form of the family
    struct CanonicalIdofsTag {};
    static constexpr BatchVariable<CanonicalIdofsTag> CanonicalIdofs{};
    // the coupling matrix of a face towards a neighbor of another family
    struct CouplingTag {};
    static constexpr BatchVariable<CouplingTag> Coupling{};
  };
};

/**
 * The variables of the batches of the Dynamic Rupture (Dr) solver.
 */
struct Dr {
  struct Id {
    struct DerivativesPlusTag {};
    static constexpr BatchVariable<DerivativesPlusTag> DerivativesPlus{};
    struct DerivativesMinusTag {};
    static constexpr BatchVariable<DerivativesMinusTag> DerivativesMinus{};
    struct IdofsPlusTag {};
    static constexpr BatchVariable<IdofsPlusTag> IdofsPlus{};
    struct IdofsMinusTag {};
    static constexpr BatchVariable<IdofsMinusTag> IdofsMinus{};
    struct QInterpolatedPlusTag {};
    static constexpr BatchVariable<QInterpolatedPlusTag> QInterpolatedPlus{};
    struct QInterpolatedMinusTag {};
    static constexpr BatchVariable<QInterpolatedMinusTag> QInterpolatedMinus{};
    struct TinvTTag {};
    static constexpr BatchVariable<TinvTTag> TinvT{};
  };
};

/**
 * The cells a batch works on, by their index in the layer.
 */
struct Indices {
  struct Id {
    struct CellsTag {};
    static constexpr BatchVariable<CellsTag> Cells{};
  };
};
} // namespace seissol::recording::inner_keys

namespace seissol::recording {
constexpr size_t AllBits = ~0_UZ;
constexpr size_t encodeAny(unsigned count) { return ~(AllBits << count); }

/// The kernels of the solvers that batches are recorded for.
struct KernelNames {
  struct TimeTag {};
  static constexpr KernelName<TimeTag> Time{};
  struct VolumeTag {};
  static constexpr KernelName<VolumeTag> Volume{};
  struct LocalFluxTag {};
  static constexpr KernelName<LocalFluxTag> LocalFlux{};
  struct NeighborFluxTag {};
  static constexpr KernelName<NeighborFluxTag> NeighborFlux{};
  struct FaceDisplacementsTag {};
  static constexpr KernelName<FaceDisplacementsTag> FaceDisplacements{};
  struct PlasticityTag {};
  static constexpr KernelName<PlasticityTag> Plasticity{};
  struct DrSpaceMapTag {};
  static constexpr KernelName<DrSpaceMapTag> DrSpaceMap{};
  struct BoundaryConditionsTag {};
  static constexpr KernelName<BoundaryConditionsTag> BoundaryConditions{};
  // the time integrals of the neighbors of another configuration
  struct ConfigBoundaryTag {};
  static constexpr KernelName<ConfigBoundaryTag> ConfigBoundary{};
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
enum struct FaceRelations : size_t { PerFace = 4, Count = 16, Any = AllBits };
enum struct DrFaceRelations : size_t { PerFace = 2, Count = 8, Any = AllBits };

enum struct ExchangeInfo : size_t {
  Buffers = 1 << 0,
  Derivatives = 1 << 1,
  Count = 2,
  Any = encodeAny(Count)
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_DATATYPES_ENCODEDCONSTANTS_H_
