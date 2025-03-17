// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_HOSTEVENT_H_
#define SEISSOL_SRC_PARALLEL_HOSTEVENT_H_

#include <cstdint>
#include <cstring>
#include <device.h>
#include <memory>
#include <optional>
namespace seissol::parallel {

class HostEvent {
  public:
  HostEvent() : managed(true), value(nullptr) {
    value = reinterpret_cast<uint32_t*>(device::DeviceInstance::getInstance().api->allocPinnedMem(
        sizeof(uint32_t), device::Destination::CurrentDevice));
#pragma omp atomic write
    *value = 0;
  }

  HostEvent(const HostEvent&) = delete;
  HostEvent& operator=(const HostEvent&) = delete;

  HostEvent(HostEvent&&) = default;
  HostEvent& operator=(HostEvent&&) = default;

  ~HostEvent() {
    if (managed) {
      device::DeviceInstance::getInstance().api->freePinnedMem(value);
      value = nullptr;
    }
  }

  HostEvent(uint32_t* location) : value(location) { *location = 0; }

  void complete() {
#pragma omp atomic write
    *value = 1;
  }

  void streamWait(void* stream) {
    device::DeviceInstance::getInstance().api->streamWaitMemory(stream, value, 1);
  }

  void hostWait() {
    while (!completed()) {
    }
  }

  bool completed() {
    // TODO: check again
    return *reinterpret_cast<volatile uint32_t*>(value) > 0;
  }

  [[nodiscard]] uint32_t* location() const { return value; }

  private:
  uint32_t* value;
  uint32_t target{1};
  bool managed{false};
};

class HostEventPool {
  public:
  HostEventPool() { buffers.emplace_back(size); }

  std::shared_ptr<HostEvent> get() {
    for (std::size_t i = 0; i < buffers.size(); ++i) {
      const auto tryGet = buffers[i].get();
      if (tryGet.has_value()) {
        return std::make_shared<HostEvent>(tryGet.value());
      }
    }
    buffers.emplace_back(size);
    return std::make_shared<HostEvent>(buffers.back().get().value());
  }

  void put(HostEvent& event) {
    for (std::size_t i = 0; i < buffers.size(); ++i) {
      if (buffers[i].put(event.location())) {
        return;
      }
    }
  }

  HostEventPool(const HostEventPool&) = delete;
  HostEventPool& operator=(const HostEventPool&) = delete;

  HostEventPool(HostEventPool&&) = default;
  HostEventPool& operator=(HostEventPool&&) = default;

  ~HostEventPool() = default;

  private:
  struct HostEventBuffer {
    HostEventBuffer(const HostEventBuffer&) = delete;
    HostEventBuffer& operator=(const HostEventBuffer&) = delete;

    HostEventBuffer(HostEventBuffer&&) = default;
    HostEventBuffer& operator=(HostEventBuffer&&) = default;

    uint32_t* data{nullptr};
    std::vector<bool> taken;

    std::optional<uint32_t*> get() {
      for (std::size_t i = 0; i < taken.size(); ++i) {
        if (!taken[i]) {
          taken[i] = true;
          return data + i;
        }
      }
      return {};
    }

    bool put(const uint32_t* location) {
      const auto ptrLocation = reinterpret_cast<intptr_t>(location);
      const auto ptrData = reinterpret_cast<intptr_t>(data);
      const auto ptrDiff = ptrLocation - ptrData;
      if (ptrDiff >= 0 && ptrDiff < taken.size()) {
        taken[ptrDiff] = false;
        return true;
      }
      return false;
    }

    HostEventBuffer(std::size_t size) : taken(size) {
      data = reinterpret_cast<uint32_t*>(
          device::DeviceInstance::getInstance().api->allocPinnedMem(size * sizeof(uint32_t)));
      std::memset(data, 0, size * sizeof(uint32_t));
    }

    ~HostEventBuffer() { device::DeviceInstance::getInstance().api->freePinnedMem(data); }
  };

  std::size_t size;
  std::vector<HostEventBuffer> buffers;
};

} // namespace seissol::parallel
#endif // SEISSOL_SRC_PARALLEL_HOSTEVENT_H_
