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
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cstddef>
#include <yateto.h>

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

namespace internal {
/// Runs the member (face, side) of the kernel family of `prototype` on every batch the recorders
/// noted under `kernelName` and `kind` for that face and side, once `bind(krnl, entry, face)` has
/// handed it the batch; `counted` is a pointer array of the batch. The batches run one after the
/// other on the stream of `runtime`, which the correctors need: those of different faces add to
/// the same cells.
template <typename KernelT, typename BindT>
void runFaceSideBatches(const KernelT& prototype,
                        std::size_t kernelName,
                        std::size_t kind,
                        recording::inner_keys::Wp::Id counted,
                        recording::ConditionalPointersToRealsTable& table,
                        parallel::runtime::StreamRuntime& runtime,
                        const BindT& bind) {
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (std::size_t side = 0; side < 2; ++side) {
      const auto found = table.find(recording::ConditionalKey(kernelName, kind, face, side));
      if (found == table.end()) {
        continue;
      }
      auto& entry = found->second;
      KernelT krnl = prototype;
      krnl.numElements = entry.get(counted)->getSize();
      auto tmpMem = runtime.memoryHandle<real>(
          (KernelT::TmpMaxMemRequiredInBytes * krnl.numElements) / sizeof(real));
      krnl.linearAllocator.initialize(tmpMem.get());
      bind(krnl, entry, face);
      krnl.streamPtr = runtime.stream();
      krnl.execute(face, side);
    }
  }
}
} // namespace internal

/// writeFaceTraces on the device, for the batches the recorder noted by face and by side (see
/// recording::LocalIntegrationRecorder).
template <typename KernelT>
void writeFaceTracesBatched(const KernelT& prototype,
                            recording::ConditionalPointersToRealsTable& table,
                            parallel::runtime::StreamRuntime& runtime) {
  using namespace seissol::recording;
  internal::runFaceSideBatches(
      prototype,
      *KernelNames::FaceTrace,
      *ComputationKind::None,
      inner_keys::Wp::Id::FaceTrace,
      table,
      runtime,
      [](KernelT& krnl, auto& entry, std::size_t /*face*/) {
        krnl.I = const_cast<const real**>(entry.get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr());
        krnl.faceTrace = entry.get(inner_keys::Wp::Id::FaceTrace)->getDeviceDataPtr();
      });
}

/// applyFaceCorrectors on the device, for the batches the recorder noted by face and by side (see
/// recording::NeighIntegrationRecorder); `bindTarget(krnl, entry)` hands the kernel what the cells
/// of a batch lift into.
template <typename KernelT, typename BindTargetT>
void applyFaceCorrectorsBatched(const KernelT& prototype,
                                recording::ConditionalPointersToRealsTable& table,
                                parallel::runtime::StreamRuntime& runtime,
                                const BindTargetT& bindTarget) {
  using namespace seissol::recording;
  internal::runFaceSideBatches(
      prototype,
      *KernelNames::FaceCorrector,
      *ComputationKind::None,
      inner_keys::Wp::Id::FaceTrace,
      table,
      runtime,
      [&](KernelT& krnl, auto& entry, std::size_t face) {
        const auto pointers = [&](inner_keys::Wp::Id id) {
          return const_cast<const real**>(entry.get(id)->getDeviceDataPtr());
        };
        krnl.faceTrace = pointers(inner_keys::Wp::Id::FaceTrace);
        krnl.faceTraceNeighbor = pointers(inner_keys::Wp::Id::FaceTraceNeighbor);
        bindCorrectorOperandsBatched(
            krnl, pointers(inner_keys::Wp::Id::LocalIntegrationData), face);
        bindTarget(krnl, entry);
      });
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_FACETRACES_H_
