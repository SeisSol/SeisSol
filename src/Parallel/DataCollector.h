// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_
#define SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_

#include <cstddef>
#include <vector>

#include "Kernels/Precision.h"

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::parallel {

// wrapper class for sparse device-to-host transfers (e.g. needed for receivers)
class DataCollector {
  public:
  DataCollector(const std::vector<real*>& indexDataHost,
                size_t elemSize,
                bool hostAccessible = false)
      : indexCount(indexDataHost.size()), indexDataHost(indexDataHost), elemSize(elemSize),
        hostAccessible(hostAccessible) {
#ifndef ACL_DEVICE
    // in case we want to use this class in a host-only scenario
    this->hostAccessible = true;
#endif
    if (!this->hostAccessible && indexCount > 0) {
#ifdef ACL_DEVICE
      indexDataDevice = reinterpret_cast<real**>(
          device::DeviceInstance::getInstance().api->allocGlobMem(sizeof(real*) * indexCount));
      device::DeviceInstance::getInstance().api->copyTo(
          indexDataDevice, indexDataHost.data(), sizeof(real*) * indexCount);
      copiedData =
          reinterpret_cast<real*>(device::DeviceInstance::getInstance().api->allocPinnedMem(
              sizeof(real) * elemSize * indexCount));
      copiedDataDevice = reinterpret_cast<real*>(
          device::DeviceInstance::getInstance().api->devicePointer(copiedData));
#endif
    }
  }

  DataCollector(const DataCollector&) = delete;

  ~DataCollector() {
    if (!hostAccessible && indexCount > 0) {
#ifdef ACL_DEVICE
      device::DeviceInstance::getInstance().api->freeMem(indexDataDevice);
      device::DeviceInstance::getInstance().api->freePinnedMem(copiedData);
#endif
    }
  }

  // gathers and sends data to the device to the host
  void gatherToHost(void* stream) {
    if (!hostAccessible && indexCount > 0) {
#ifdef ACL_DEVICE
      device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
          indexDataDevice, copiedDataDevice, elemSize, elemSize, indexCount, stream);
#endif
    }
  }

  // sends and scatters data from the host to the device
  void scatterFromHost(void* stream) {
    if (!hostAccessible && indexCount > 0) {
#ifdef ACL_DEVICE
      device::DeviceInstance::getInstance().algorithms.copyUniformToScatter(
          copiedDataDevice, indexDataDevice, elemSize, elemSize, indexCount, stream);
#endif
    }
  }

  real* get(size_t index) {
    if (hostAccessible) {
      return indexDataHost[index];
    } else {
      return copiedData + index * elemSize;
    }
  }

  const real* get(size_t index) const {
    if (hostAccessible) {
      return indexDataHost[index];
    } else {
      return copiedData + index * elemSize;
    }
  }

  private:
  bool hostAccessible;
  real** indexDataDevice;
  std::vector<real*> indexDataHost;
  size_t indexCount;
  size_t elemSize;
  real* copiedData;
  real* copiedDataDevice;
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_
