// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MEMORY_HOSTSNAPSHOT_H_
#define SEISSOL_SRC_MEMORY_HOSTSNAPSHOT_H_

#include "Memory/HostData.h"

#include <cstddef>
#include <vector>

namespace seissol::initializer {

/**
 * Copies data from the device to the host while the device goes on and changes it: the data gets
 * copied within the device first, into a buffer of the snapshot, which is fast; and from there to
 * the host later, e.g. on the thread of the outputs.
 *
 * The buffer grows to the largest snapshot, as far as the memory of the device allows (see take()).
 * The host memory has to be pinned: a copy into other memory blocks the thread until the device
 * gets to it, and the driver holds up the other threads meanwhile (e.g. the main thread when it
 * enqueues host functions).
 */
class HostSnapshot {
  public:
  HostSnapshot() = default;
#ifdef ACL_DEVICE
  ~HostSnapshot();
#else
  ~HostSnapshot() = default;
#endif

  HostSnapshot(const HostSnapshot&) = delete;
  HostSnapshot& operator=(const HostSnapshot&) = delete;
  HostSnapshot(HostSnapshot&&) = delete;
  HostSnapshot& operator=(HostSnapshot&&) = delete;

  /**
   * Copies the device parts of `copies` into the buffer, on `stream`, and waits for that; the last
   * snapshot has to be copied to the host before (see copyToHost()).
   *
   * False, with nothing done, if the buffer would not fit into the memory of the device (it leaves
   * a tenth of that, and at least 2 GiB, to what SeisSol does not allocate itself, e.g. the context
   * and the communication libraries), or if host memory of a copy is not pinned.
   */
  bool take(const std::vector<HostCopy>& copies, void* stream);

  /// Copies the last snapshot to the host parts of its copies, and waits for that; on any thread.
  void copyToHost();

#ifdef ACL_DEVICE
  private:
  // the copies of the last snapshot, from the buffer
  std::vector<HostCopy> pending_;
  void* buffer_{nullptr};
  std::size_t capacity_{0};
  void* stream_{nullptr};
  bool warnedSize_{false};
  bool warnedPinned_{false};
#endif
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_HOSTSNAPSHOT_H_
