// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MEMORY_HOSTDATA_H_
#define SEISSOL_SRC_MEMORY_HOSTDATA_H_

#include <cstddef>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace seissol::initializer {

/// A copy of device memory to its counterpart on the host.
struct HostCopy {
  void* host{nullptr};
  const void* device{nullptr};
  std::size_t size{0};
  //! whether the host memory is pinned: a copy into other memory blocks until it is done
  bool pinned{false};
};

/**
 * The variables of the storages that the host reads, e.g. the outputs due at a synchronization
 * point. Each gets copied to the host once, however many ask for it.
 *
 * Everything (all variables of all storages) is what a reader gets that does not say what it
 * reads.
 */
class HostData {
  public:
  /// all variables of all storages
  static HostData everything();

  /// the variable (or bucket) `StorageT` of `storage`
  template <typename StorageT, typename StorageType>
  HostData& add(StorageType& storage) {
    const auto variable = storage.template info<StorageT>().index;
    variables_.try_emplace(std::make_pair(static_cast<const void*>(&storage), variable),
                           [&storage](std::vector<HostCopy>& copies) {
                             storage.template varHostCopies<StorageT>(copies);
                           });
    return *this;
  }

  /// all variables of `storage` (not its buckets and scratchpads)
  template <typename StorageType>
  HostData& addVariables(StorageType& storage) {
    variables_.try_emplace(
        std::make_pair(static_cast<const void*>(&storage), AllVariables),
        [&storage](std::vector<HostCopy>& copies) { storage.variablesHostCopies(copies); });
    return *this;
  }

  /// adds what `other` holds
  HostData& merge(const HostData& other);

  [[nodiscard]] bool isEverything() const { return everything_; }

  [[nodiscard]] bool empty() const { return !everything_ && variables_.empty(); }

  /// the number of variables; zero for everything
  [[nodiscard]] std::size_t size() const { return variables_.size(); }

  /// the copies to the host the variables need (none for those the host and the device share);
  /// not for everything
  [[nodiscard]] std::vector<HostCopy> copies() const;

  /// enqueues the copies to the host on the stream; not for everything
  void copyToHost(void* stream) const;

  private:
  // the key of all variables of a storage
  static constexpr std::size_t AllVariables = static_cast<std::size_t>(-1);

  bool everything_{false};
  // by the storage and the index of the variable in it: what adds the copies of the variable
  std::map<std::pair<const void*, std::size_t>, std::function<void(std::vector<HostCopy>&)>>
      variables_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_HOSTDATA_H_
