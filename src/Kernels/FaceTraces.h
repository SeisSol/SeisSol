// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_FACETRACES_H_
#define SEISSOL_SRC_KERNELS_FACETRACES_H_

#include "Common/Constants.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/Condition.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "GeneratedCode/init.h"
#endif

namespace seissol::tensor {
// generated where a target takes the side of a face as data, on a device
struct traceSide;
} // namespace seissol::tensor

/// The corrector of a cell works from face traces: after the predictor, a cell writes the trace of
/// its time integral on each of its faces (LTS::FaceTraces), in the face basis; the corrector of a
/// face applies the flux to that trace and the one the neighbor wrote on the same face, and lifts
/// the result into the cell. See _configureTraces in the generator.
namespace seissol::kernels {

/// Writes the trace of the time integral of a cell on each of its faces (LTS::FaceTraces), each in
/// the parametrisation of the face its traces take (CellLocalInformation::traceSides).
template <typename KernelT>
void writeFaceTraces(const KernelT& prototype, const real* timeIntegrated, LTS::Ref& data) {
  KernelT krnl = prototype;
  krnl.I = timeIntegrated;
  const auto& cellInformation = data.get<LTS::CellInformation>();
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    krnl.faceTrace = data.get<LTS::FaceTraces>()[face];
    krnl.execute(face, cellInformation.traceSides[face]);
  }
}

/// Hands a kernel the traces of the time derivatives of a cell on a face, at `faceDerivatives`:
/// one trace per derivative, each with the modes it reaches, one after the other
/// (SpacetimeKernel::computeFaceDerivatives); for a predictor without derivatives, its one trace.
template <typename KernelT, typename PtrT>
void bindFaceDerivatives(KernelT& krnl, PtrT* faceDerivatives) {
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::faceDerivatives>(); ++i) {
    krnl.faceDerivatives(i) =
        faceDerivatives + yateto::computeFamilySize<tensor::faceDerivatives>(1, i);
  }
}

/// The same for a batch, where the traces of a member are an offset into what each pointer holds.
template <typename KernelT, typename PtrT>
void bindFaceDerivativesBatched(KernelT& krnl, PtrT** faceDerivatives) {
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::faceDerivatives>(); ++i) {
    krnl.faceDerivatives(i) = faceDerivatives;
    krnl.extraOffset_faceDerivatives(i) = yateto::computeFamilySize<tensor::faceDerivatives>(1, i);
  }
}

/// Hands the corrector of a face what it applies to the jump of the two traces of the face, in
/// whichever of the two shapes the cell carries it: the rotation of the face together with the one
/// set of scalars at its points, or the matrix of the cell's own side, with the rotation folded in
/// (see initializer::toCorrectorForm).
template <typename KernelT, typename LocalIntegrationT>
void bindCorrectorOperands(KernelT& krnl,
                           const LocalIntegrationT& localIntegration,
                           std::size_t face) {
  if constexpr (NodalFlux) {
    bindFaceRotation(krnl, localIntegration.faceRotation[face]);
    for (std::size_t coefficient = 0; coefficient < FluxCoefficientCount; ++coefficient) {
      krnl.fluxCoefficients(coefficient) = localIntegration.fluxCoefficients[face][coefficient];
    }
  } else {
    krnl.AplusT = localIntegration.nApNm1[face];
  }
}

/// Applies the corrector of every face of a cell that takes one (appliesFaceCorrector) into what
/// the kernel `prototype` has bound as its target: the cell's own trace of the face against
/// `neighborTraces`, the trace the neighbor wrote there over the cell's interval. A face without a
/// neighbor reads a trace of zeros there (LTS::ZeroFaceTrace): a boundary face, whose own flux the
/// face carries in its operator, and a fault face, where the corrector subtracts the normal flux of
/// the cell's own trace alone.
template <typename KernelT>
void applyFaceCorrectors(const KernelT& prototype,
                         const std::array<const real*, Cell::NumFaces>& neighborTraces,
                         LTS::Ref& data) {
  KernelT krnl = prototype;
  const auto& cellInformation = data.get<LTS::CellInformation>();
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (!appliesFaceCorrector(cellInformation.faceTypes[face])) {
      continue;
    }
    krnl.faceTrace = data.get<LTS::FaceTraces>()[face];
    krnl.faceTraceNeighbor = neighborTraces[face];
    bindCorrectorOperands(krnl, data.get<LTS::LocalIntegration>(), face);
    krnl.execute(face, cellInformation.traceSides[face]);
  }
}

/// The reals one operator of the side of a face takes (LTS::TraceSideOperators), as the device
/// kernels that take the side as data read it (tensor::traceSide); none where no kernel does.
constexpr std::size_t TraceSideOperatorSize = kernels::size<tensor::traceSide>(0);

/// Writes the operators of the two sides of a face, the identity and the reflection, one after the
/// other, where the device kernels that take the side as data read them (LTS::TraceSideOperators).
inline void writeTraceSideOperators([[maybe_unused]] real* operators) {
#ifdef ACL_DEVICE
  for (std::size_t side = 0; side < 2; ++side) {
    std::copy_n(init::traceSideOperator::Values[side],
                TraceSideOperatorSize,
                operators + side * TraceSideOperatorSize);
  }
#endif
}

namespace internal {
/// The pointers of one face in an array the recorders fill face after face for the cells of a
/// batch (see recording::LocalIntegrationRecorder): those of face `face` start `face * cells` in.
template <typename T>
T** ofFace(T** pointers, std::size_t face, std::size_t cells) {
  return pointers + face * cells;
}

/// Runs the member `face` of the kernel family of `prototype` on every batch the recorders noted
/// under `kernelName` and `kind` for that face, once `bind(krnl, entry, face)` has handed it the
/// batch, with the operator of the side each cell of the batch is on
/// (inner_keys::Wp::Id::TraceSide); `counted` is a pointer array of the batch. The batches run one
/// after the other on the stream of `runtime`.
template <typename KernelT, typename BindT>
void runFaceBatches(const KernelT& prototype,
                    std::size_t kernelName,
                    std::size_t kind,
                    recording::inner_keys::Wp::Id counted,
                    recording::ConditionalPointersToRealsTable& table,
                    parallel::runtime::StreamRuntime& runtime,
                    const BindT& bind) {
  using namespace seissol::recording;
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const auto found = table.find(ConditionalKey(kernelName, kind, face));
    if (found == table.end()) {
      continue;
    }
    auto& entry = found->second;
    KernelT krnl = prototype;
    krnl.numElements = entry.get(counted)->getSize();
    auto tmpMem = runtime.memoryHandle<real>(
        (KernelT::TmpMaxMemRequiredInBytes * krnl.numElements) / sizeof(real));
    krnl.linearAllocator.initialize(tmpMem.get());
    krnl.traceSide(face) =
        const_cast<const real**>(entry.get(inner_keys::Wp::Id::TraceSide)->getDeviceDataPtr());
    bind(krnl, entry, face);
    krnl.streamPtr = runtime.stream();
    krnl.execute(face);
  }
}
} // namespace internal

/// writeFaceTraces on the device: one kernel writes the traces of the four faces of every cell of
/// the batch the recorder noted, each face with the operator of the side its cell is on there
/// (see recording::LocalIntegrationRecorder).
template <typename KernelT>
void writeFaceTracesBatched(const KernelT& prototype,
                            recording::ConditionalPointersToRealsTable& table,
                            parallel::runtime::StreamRuntime& runtime) {
  using namespace seissol::recording;
  const auto found = table.find(ConditionalKey(*KernelNames::FaceTrace, *ComputationKind::None));
  if (found == table.end()) {
    return;
  }
  auto& entry = found->second;
  KernelT krnl = prototype;
  const auto cells = entry.get(inner_keys::Wp::Id::Idofs)->getSize();
  krnl.numElements = cells;
  krnl.I = const_cast<const real**>(entry.get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr());
  auto** traces = entry.get(inner_keys::Wp::Id::FaceTrace)->getDeviceDataPtr();
  auto** sides =
      const_cast<const real**>(entry.get(inner_keys::Wp::Id::TraceSide)->getDeviceDataPtr());
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    krnl.faceTraces(face) = internal::ofFace(traces, face, cells);
    krnl.traceSide(face) = internal::ofFace(sides, face, cells);
  }
  auto tmpMem =
      runtime.memoryHandle<real>((KernelT::TmpMaxMemRequiredInBytes * cells) / sizeof(real));
  krnl.linearAllocator.initialize(tmpMem.get());
  krnl.streamPtr = runtime.stream();
  krnl.execute();
}

/// applyFaceCorrectors on the device: one kernel applies the correctors of the four faces of every
/// cell of the batch the recorder noted, each face with the operator of the side its cell is on
/// there (see recording::NeighIntegrationRecorder); `bindTarget(krnl, entry)` hands the kernel what
/// the cells lift into. A face that takes no corrector (appliesFaceCorrector) reads a trace of
/// zeros on both sides, which its corrector turns into nothing.
template <typename KernelT, typename BindTargetT>
void applyFaceCorrectorsBatched(const KernelT& prototype,
                                recording::ConditionalPointersToRealsTable& table,
                                parallel::runtime::StreamRuntime& runtime,
                                const BindTargetT& bindTarget) {
  using namespace seissol::recording;
  const auto found =
      table.find(ConditionalKey(*KernelNames::FaceCorrector, *ComputationKind::None));
  if (found == table.end()) {
    return;
  }
  auto& entry = found->second;
  const auto pointers = [&](inner_keys::Wp::Id id) {
    return const_cast<const real**>(entry.get(id)->getDeviceDataPtr());
  };
  KernelT krnl = prototype;
  const auto cells = entry.get(inner_keys::Wp::Id::FluxTarget)->getSize();
  krnl.numElements = cells;
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    krnl.faceTraces(face) = internal::ofFace(pointers(inner_keys::Wp::Id::FaceTrace), face, cells);
    krnl.faceTracesNeighbor(face) =
        internal::ofFace(pointers(inner_keys::Wp::Id::FaceTraceNeighbor), face, cells);
    krnl.traceSide(face) = internal::ofFace(pointers(inner_keys::Wp::Id::TraceSide), face, cells);
  }
  bindFaceCorrectorOperandsBatched(krnl, pointers(inner_keys::Wp::Id::LocalIntegrationData));
  bindTarget(krnl, entry);
  auto tmpMem =
      runtime.memoryHandle<real>((KernelT::TmpMaxMemRequiredInBytes * cells) / sizeof(real));
  krnl.linearAllocator.initialize(tmpMem.get());
  krnl.streamPtr = runtime.stream();
  krnl.execute();
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_FACETRACES_H_
