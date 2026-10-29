// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_HOSTBUFFERPOOL_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_HOSTBUFFERPOOL_H_

#include <cstddef>
#include <mutex>
#include <vector>

namespace seissol::parallel::runtime {

/**
 * Buffers of pinned host memory of one size, which the device writes for work on the host (see
 * HostActor): one for each piece of work that has not run yet. A buffer gets taken before the
 * device work that writes it is enqueued, and given back by the host work once it has read it. As
 * the device never waits for the host, there are as many buffers as pieces of work are in flight
 * at once; the pool grows to that, and keeps its buffers until it goes.
 *
 * Buffers may be given back on any thread (e.g. by host work that runs on the thread of the
 * outputs). Without a device, the buffers are plain host memory.
 */
class HostBufferPool {
  public:
  struct Buffer {
    /// where the host reads it
    void* host{nullptr};
    /// where the device writes it
    void* device{nullptr};
  };

  explicit HostBufferPool(std::size_t bytes);
  ~HostBufferPool();

  HostBufferPool(const HostBufferPool&) = delete;
  HostBufferPool& operator=(const HostBufferPool&) = delete;
  HostBufferPool(HostBufferPool&&) = delete;
  HostBufferPool& operator=(HostBufferPool&&) = delete;

  /// a buffer that nobody uses; a new one if there is none
  Buffer take();

  /// gives a buffer back once nobody reads it any more
  void give(const Buffer& buffer);

  [[nodiscard]] std::size_t bytes() const { return bytes_; }

  /// the buffers of the pool, used or not; for the thread that takes them
  [[nodiscard]] std::size_t size() const { return all_.size(); }

  private:
  std::size_t bytes_;
  std::vector<Buffer> all_;
  std::vector<Buffer> free_;
  std::mutex mutex_;
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_HOSTBUFFERPOOL_H_
