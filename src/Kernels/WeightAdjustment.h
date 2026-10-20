// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_WEIGHTADJUSTMENT_H_
#define SEISSOL_SRC_KERNELS_WEIGHTADJUSTMENT_H_

#include "Common/Marker.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/device.h>
#include <yateto.h>
#endif

namespace seissol::kernels {

/// Adds the inverse weight of the medium applied to the flux buffer of a weight-adjusted cell to
/// the field its faces lift into otherwise (see WeightAdjustable). A class of its own, so that a
/// build without the kernels never instantiates what names their operands.
///
/// On the device, the weight-adjusted cells of a layer form a batch of their own
/// (recording::KernelNames::WeightAdjust): the local integral clears their buffers before their
/// faces lift into them, and the neighbour integral applies the weight once all of them are in.
template <bool Enabled, typename KernelT, typename DeviceKernelT>
class WeightAdjustment {
  public:
  void setGlobalData(const CompoundGlobalData& global) {
    prototype_.bindGlobals(*global.onHost);
#ifdef ACL_DEVICE
    devicePrototype_.bindGlobals(*global.onDevice);
#endif
  }

  /// `bindTarget` hands the kernel that field: the degrees of freedom, or the extended field of a
  /// solver that keeps the memory variables apart.
  template <typename BindT>
  void apply(LTS::Ref& data, const BindT& bindTarget) const {
    KernelT krnl = prototype_;
    bindTarget(krnl);
    krnl.fluxBuffer = data.get<LTS::FluxBuffer>();
    kernels::bindWeightOperands(krnl, data.get<LTS::LocalIntegration>());
    krnl.execute();
  }

  /// Clears the flux buffers of the weight-adjusted cells of a batch.
  static void clearBatched(SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
                           SEISSOL_GPU_PARAM parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
    const recording::ConditionalKey key(*recording::KernelNames::WeightAdjust);
    if (table.find(key) != table.end()) {
      auto* buffers = table[key].get(recording::inner_keys::Wp::Id::FluxBuffer);
      device::DeviceInstance::getInstance().algorithms.setToValue(
          buffers->getDeviceDataPtr(),
          static_cast<real>(0),
          kernels::size<tensor::fluxBuffer>(),
          buffers->getSize(),
          runtime.stream());
    }
#endif
  }

  /// The same as apply for the weight-adjusted cells of a batch; `bindTarget` hands the kernel the
  /// field from the entry of the batch.
  template <typename BindT>
  void applyBatched(SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
                    SEISSOL_GPU_PARAM parallel::runtime::StreamRuntime& runtime,
                    SEISSOL_GPU_PARAM const BindT& bindTarget) const {
#ifdef ACL_DEVICE
    const recording::ConditionalKey key(*recording::KernelNames::WeightAdjust);
    if (table.find(key) != table.end()) {
      auto& entry = table[key];
      const auto numElements = entry.get(recording::inner_keys::Wp::Id::FluxBuffer)->getSize();
      DeviceKernelT krnl = devicePrototype_;
      krnl.numElements = numElements;
      bindTarget(krnl, entry);
      krnl.fluxBuffer = const_cast<const real**>(
          entry.get(recording::inner_keys::Wp::Id::FluxBuffer)->getDeviceDataPtr());
      kernels::bindWeightOperandsBatched(
          krnl,
          const_cast<const real**>(
              entry.get(recording::inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
      auto tmpMem = runtime.memoryHandle<real>((yateto::getMaxTmpMemRequired(krnl) * numElements) /
                                               sizeof(real));
      krnl.linearAllocator.initialize(tmpMem.get());
      krnl.streamPtr = runtime.stream();
      krnl.execute();
    }
#endif
  }

  [[nodiscard]] static PerformanceEstimate metrics() {
    return PerformanceEstimate::fromKernel<KernelT>();
  }

  private:
  KernelT prototype_;
#ifdef ACL_DEVICE
  DeviceKernelT devicePrototype_;
#endif
};

template <typename KernelT, typename DeviceKernelT>
class WeightAdjustment<false, KernelT, DeviceKernelT> {
  public:
  void setGlobalData(const CompoundGlobalData& /*global*/) {}
  template <typename BindT>
  void apply(LTS::Ref& /*data*/, const BindT& /*bindTarget*/) const {}
  static void clearBatched(recording::ConditionalPointersToRealsTable& /*table*/,
                           parallel::runtime::StreamRuntime& /*runtime*/) {}
  template <typename BindT>
  void applyBatched(recording::ConditionalPointersToRealsTable& /*table*/,
                    parallel::runtime::StreamRuntime& /*runtime*/,
                    const BindT& /*bindTarget*/) const {}
  [[nodiscard]] static PerformanceEstimate metrics() { return {}; }
};

/// The one of this build.
using WeightAdjustmentT =
    WeightAdjustment<WeightAdjustable, kernel::weightAdjust, kernel::gpu_weightAdjust>;

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_WEIGHTADJUSTMENT_H_
