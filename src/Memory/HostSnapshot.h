// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MEMORY_HOSTSNAPSHOT_H_
#define SEISSOL_SRC_MEMORY_HOSTSNAPSHOT_H_

#include <vector>

namespace seissol::initializer {

struct DualMemoryContainer;

/**
 * The copies to the host of a synchronization point, through the snapshots of the data (see
 * DualMemoryContainer::snapshot), while the device goes on and changes the data: take() takes the
 * snapshots, which is fast, and copyToHost() copies those that are apart from the host parts
 * there, later and e.g. on the thread of the outputs. On unified memory, the snapshots are the
 * host parts, and copyToHost() has nothing to do.
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
   * Takes the snapshots of the memory of `containers` on `stream`, copies the memory without one to
   * the host right away, and waits for that; the last snapshots have to be copied to the host
   * before (see copyToHost()).
   */
  void take(const std::vector<const DualMemoryContainer*>& containers, void* stream);

  /// Copies the snapshots of the last take() to the host parts, and waits for that; on any thread.
  void copyToHost();

#ifdef ACL_DEVICE
  private:
  // the memory whose snapshot still has to be copied to the host
  std::vector<const DualMemoryContainer*> pending_;
  void* stream_{nullptr};
#endif
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_HOSTSNAPSHOT_H_
