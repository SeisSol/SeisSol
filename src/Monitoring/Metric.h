// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MONITORING_METRIC_H_
#define SEISSOL_SRC_MONITORING_METRIC_H_

#include <cstdint>

namespace seissol {

/// Whose hardware FLOP and kernel bytes the estimates of the kernels count. On the host, those of
/// the host kernels, as yateto counts what their GEMM generators issue. On the device, those of
/// the device kernels where the device code generator reports the arithmetic it emits, which
/// TensorForge does (SEISSOL_DEVICE_CODEGEN_REPORTS_FLOP); gemmforge and chainforge report none,
/// and what yateto counts of the GEMMs it hands to them is not what they carry out, so a build with
/// them keeps counting the host kernels.
enum class HardwareFlopSource { HostKernels, DeviceKernels };

#if defined(ACL_DEVICE) && defined(SEISSOL_DEVICE_CODEGEN_REPORTS_FLOP)
constexpr auto KernelHardwareFlopSource = HardwareFlopSource::DeviceKernels;
#else
constexpr auto KernelHardwareFlopSource = HardwareFlopSource::HostKernels;
#endif

/// The name of a HardwareFlopSource in the logs and in the run summary.
constexpr auto hardwareFlopSourceName(HardwareFlopSource source) -> const char* {
  return source == HardwareFlopSource::DeviceKernels ? "device-kernels" : "host-kernels";
}

/**
  Collects metrics for measuring performance of the simulation.
  */
struct PerformanceEstimate {
  // flop count by instructions executed in hardware
  std::uint64_t hardwareFlop{0};

  // (logical) flop count for fully sparse, unpadded operations (i.e. only the operations that we
  // really "need")
  std::uint64_t nonzeroFlop{0};

  // legacy value; bytes transferred in/out by CPU kernels
  std::uint64_t bytes{0};

  // non-const bytes transferred in/out by generated (Yateto) kernels. Currently mostly relevant for
  // the GPU (due to kernel == loop correspondence).
  std::uint64_t kernelBytes{0};

  constexpr auto operator+(const PerformanceEstimate& other) const -> PerformanceEstimate {
    return PerformanceEstimate{.hardwareFlop = hardwareFlop + other.hardwareFlop,
                               .nonzeroFlop = nonzeroFlop + other.nonzeroFlop,
                               .bytes = bytes + other.bytes,
                               .kernelBytes = kernelBytes + other.kernelBytes};
  }

  constexpr auto operator+=(const PerformanceEstimate& other) -> PerformanceEstimate& {
    *this = *this + other;
    return *this;
  }

  constexpr auto operator*(std::uint64_t other) const -> PerformanceEstimate {
    return PerformanceEstimate{.hardwareFlop = hardwareFlop * other,
                               .nonzeroFlop = nonzeroFlop * other,
                               .bytes = bytes * other,
                               .kernelBytes = kernelBytes * other};
  }

  constexpr auto operator*=(std::uint64_t other) -> PerformanceEstimate& {
    *this = *this * other;
    return *this;
  }

  template <typename T, typename... Args>
  constexpr static auto fromKernel(Args... kargs) -> PerformanceEstimate {
    return PerformanceEstimate{T::hardwareFlops(kargs...),
                               T::nonZeroFlops(kargs...),
                               0,
                               T::outboundBytes(kargs...) + T::inboundBytes(kargs...)};
  }

  /// The estimate of a kernel that a device build runs as its twin `DeviceT`: the hardware FLOP
  /// and the kernel bytes of the kernel that runs where KernelHardwareFlopSource says so. The
  /// nonzero FLOP stay those of the host kernel: the operation they count is the same in every
  /// build, and the device kernels count it in the layout their code generator chose, so builds
  /// with different code generators would differ in how much useful work they did. Called through
  /// SEISSOL_KERNEL_ESTIMATE.
  template <typename HostT, typename DeviceT, typename... Args>
  constexpr static auto fromKernels(Args... kargs) -> PerformanceEstimate {
    auto estimate = fromKernel<HostT>(kargs...);
    if constexpr (KernelHardwareFlopSource == HardwareFlopSource::DeviceKernels) {
      estimate.hardwareFlop = DeviceT::hardwareFlops(kargs...);
      estimate.kernelBytes = DeviceT::outboundBytes(kargs...) + DeviceT::inboundBytes(kargs...);
    }
    return estimate;
  }
};

} // namespace seissol

// The estimate of a host kernel and of its twin on the device (PerformanceEstimate::fromKernels),
// to be called with the family index of the kernels; a macro, as a build for the host has no
// device kernels to name.
#ifdef ACL_DEVICE
#define SEISSOL_KERNEL_ESTIMATE(HostKernel, DeviceKernel)                                          \
  ::seissol::PerformanceEstimate::fromKernels<HostKernel, DeviceKernel>
#else
#define SEISSOL_KERNEL_ESTIMATE(HostKernel, DeviceKernel)                                          \
  ::seissol::PerformanceEstimate::fromKernel<HostKernel>
#endif
#endif // SEISSOL_SRC_MONITORING_METRIC_H_
