#ifndef SEISSOL_ENCODINGCONSTANTS_HPP
#define SEISSOL_ENCODINGCONSTANTS_HPP

#include <cstdlib>

namespace seissol {
namespace initializers {
namespace recording {

enum struct EntityId : size_t {
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
  Displacements,
  NodalStressTensor,
  Pstrains,
  ElementsIds,
  InitialLoad,
  DrDerivativesPlus,
  DrDerivativesMinus,
  DrIdofsPlus,
  DrIdofsMinus,
  DrQInterpolatedPlus,
  DrQInterpolatedMinus,
  DrTinvT,
  Count
};

constexpr size_t ALL_BITS = ~static_cast<size_t>(0);
constexpr size_t encodeAny(unsigned count) {
  return ~(ALL_BITS << count);
}

enum struct KernelNames : size_t {
  Time = 1 << 0,
  Volume = 1 << 1,
  LocalFlux = 1 << 2,
  NeighborFlux = 1 << 3,
  Displacements = 1 << 4,
  Plasticity = 1 << 5,
  DrTime = 1 << 6,
  DrSpaceMap = 1 << 7,
  Count = 8,
  Any = encodeAny(Count)
};

enum struct ComputationKind : size_t {
  WithoutDerivatives = 1 << 0,
  WithDerivatives = 1 << 1,
  WithLtsDerivatives = 1 << 2,
  WithGtsDerivatives = 1 << 3,
  WithGtsBuffers = 1 << 4,
  WithLtsBuffers = 1 << 5,
  Count = 6,
  Any = encodeAny(Count)
};

enum struct FaceKinds : size_t {
  Regular = 1 << 0,
  FreeSurface = 1 << 1,
  Outflow = 1 << 2,
  DynamicRupture = 1 << 3,
  Periodic = 1 << 4,
  Count = 5,
  Any = encodeAny(Count)
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

} // namespace recording
} // namespace initializers
} // namespace seissol

#endif // SEISSOL_ENCODINGCONSTANTS_HPP
