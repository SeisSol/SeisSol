// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_RTCGPUDRIVER_H_
#define SEISSOL_SRC_EXPR_RTCGPUDRIVER_H_

// Runtime compilation and launch for the device backends.
//
// Split from RtcGpu.h on purpose. That header produces a string and answers
// whether a program can go to a device; it pulls in no vendor headers, so it
// builds and is TESTED everywhere -- the emitted source is compiled for the
// host and checked against the interpreter. This file is the part that cannot
// be: it needs nvrtc/hiprtc and the driver API, and it can only be exercised on
// hardware.
//
// Keeping the boundary here means a bug in the arithmetic is caught by a test
// that runs in ordinary CI, and only the mechanical part -- module loading,
// launch geometry, stream association -- waits for a GPU.
//
// When the toolkit is absent this compiles to a stub that returns nullptr, and
// makeKernel falls back to the interpreter with a warning naming the reason.
// The header itself is always safe to include.

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcGpu.h"

#include <memory>
#include <string>

namespace seissol::expr {

/// Whether this build has a runtime compiler for `target` at all.
[[nodiscard]] bool gpuRuntimeAvailable(GpuTarget target);

/// Compile and return a device kernel, or nullptr.
///
/// Never throws and never aborts on a compilation failure: a model that cannot
/// be compiled still has to run, so every failure path warns and returns null
/// for makeKernel to fall back on. The one thing it does NOT do is fall back
/// silently -- an unexplained drop to the interpreter is a log line nobody
/// acts on.
[[nodiscard]] std::unique_ptr<Kernel> makeRtcGpuKernel(const Program& program,
                                                       Binding& binding,
                                                       LoweredProgram lowered,
                                                       const BackendOptions& options,
                                                       GpuTarget target);

/// Number of distinct device kernels currently held.
[[nodiscard]] std::size_t rtcGpuCacheSize();

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_RTCGPUDRIVER_H_
