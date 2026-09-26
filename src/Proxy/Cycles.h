// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_CYCLES_H_
#define SEISSOL_SRC_PROXY_CYCLES_H_

#include <cstdint>
#include <string_view>

#if defined(__x86_64__) || defined(__i386__)
#define SEISSOL_PROXY_CYCLES_TSC
#elif defined(__aarch64__)
#define SEISSOL_PROXY_CYCLES_CNTVCT
#endif

namespace seissol::proxy {

/// Where the tick count comes from, so that a consumer of the numbers can tell
/// what they mean instead of having to guess from the platform.
enum class CycleSource : std::uint8_t {
  /// No counter on this platform; tick counts are reported as zero.
  None,
  /// x86 time-stamp counter, read with the RDTSC instruction.
  Tsc,
  /// AArch64 virtual counter CNTVCT_EL0, readable from EL0 on Linux.
  Cntvct,
};

/// Both counters available here run at a fixed reference frequency rather than
/// at the core clock, so a tick is not a core cycle: the counts do not follow
/// frequency scaling or turbo. They are a stable, cheap time base, and ratios
/// against them are comparable across runs on the same machine, which is what
/// the proxy reports them for.
constexpr auto cycleSource() -> CycleSource {
#if defined(SEISSOL_PROXY_CYCLES_TSC)
  return CycleSource::Tsc;
#elif defined(SEISSOL_PROXY_CYCLES_CNTVCT)
  return CycleSource::Cntvct;
#else
  return CycleSource::None;
#endif
}

constexpr auto cycleSourceName() -> std::string_view {
  switch (cycleSource()) {
  case CycleSource::Tsc:
    return "tsc";
  case CycleSource::Cntvct:
    return "cntvct";
  default:
    return "none";
  }
}

constexpr auto cyclesAvailable() -> bool { return cycleSource() != CycleSource::None; }

/// Read the tick counter. Returns zero where no counter is available, which
/// makes an elapsed count of zero the signal that no measurement happened.
inline auto readCycles() -> std::uint64_t {
#if defined(SEISSOL_PROXY_CYCLES_TSC)
  // RDTSC directly instead of the __rdtsc() intrinsic: NVHPC does not declare
  // the intrinsic in its x86intrin.h, and __builtin_ia32_rdtsc() compiles there
  // but does not link
  std::uint32_t low = 0;
  std::uint32_t high = 0;
  asm volatile("rdtsc" : "=a"(low), "=d"(high));
  return (static_cast<std::uint64_t>(high) << 32U) | low;
#elif defined(SEISSOL_PROXY_CYCLES_CNTVCT)
  std::uint64_t value = 0;
  asm volatile("mrs %0, cntvct_el0" : "=r"(value));
  return value;
#else
  return 0;
#endif
}

/// Frequency of the tick counter in Hz, or zero when it cannot be queried.
/// AArch64 exposes it directly; on x86 the nominal TSC rate is not available
/// without reading model-specific registers, which needs privileges the proxy
/// does not have.
inline auto cycleFrequency() -> std::uint64_t {
#if defined(SEISSOL_PROXY_CYCLES_CNTVCT)
  std::uint64_t value = 0;
  asm volatile("mrs %0, cntfrq_el0" : "=r"(value));
  return value;
#else
  return 0;
#endif
}

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_CYCLES_H_
