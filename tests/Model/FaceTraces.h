// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_FACETRACES_H_
#define SEISSOL_TESTS_MODEL_FACETRACES_H_

#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/FaceTraces.h"
#include "Kernels/Precision.h"

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace seissol::unit_test {

/// What a cell writes on one of its faces (LTS::FaceTraces).
using FaceTraceT = std::array<real, tensor::faceTrace::size()>;

namespace facecorrector {
/// Whether a kernel lifts into the extended field of a solver that keeps the memory variables
/// apart, rather than into the field.
template <typename KernelT, typename = void>
struct LiftsIntoExtended : std::false_type {};

template <typename KernelT>
struct LiftsIntoExtended<KernelT, std::void_t<decltype(std::declval<KernelT>().Qext)>>
    : std::true_type {};
} // namespace facecorrector

/// Hands a kernel what its faces lift into: the field, or the extended field.
template <typename KernelT>
void bindLiftTarget(KernelT& krnl, real* target) {
  if constexpr (facecorrector::LiftsIntoExtended<KernelT>::value) {
    krnl.Qext = target;
  } else {
    krnl.Q = target;
  }
}

/// The trace of `timeIntegrated` on `face`, in the parametrisation of side `side` of the face, as
/// a cell writes it (kernels::writeFaceTraces).
inline void
    writeFaceTrace(const real* timeIntegrated, std::size_t face, std::size_t side, real* trace) {
  kernel::faceTrace krnl{};
  krnl.bindGlobals(seissol::Pool::host());
  krnl.I = timeIntegrated;
  krnl.faceTrace = trace;
  krnl.execute(face, side);
}

/// Adds the corrector of `face` to `target`, as a cell applies it (kernels::applyFaceCorrectors):
/// from `ownTrace`, the trace the cell wrote there, and `neighborTrace`, the one the neighbor wrote
/// there, both in the parametrisation of side `side` of the face.
template <typename LocalIntegrationT>
void applyFaceCorrector(const real* ownTrace,
                        const real* neighborTrace,
                        const LocalIntegrationT& local,
                        std::size_t face,
                        std::size_t side,
                        real* target) {
  kernel::faceCorrector krnl{};
  krnl.bindGlobals(seissol::Pool::host());
  bindLiftTarget(krnl, target);
  krnl.faceTrace = ownTrace;
  krnl.faceTraceNeighbor = neighborTrace;
  kernels::bindCorrectorOperands(krnl, local, face);
  krnl.execute(face, side);
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_FACETRACES_H_
