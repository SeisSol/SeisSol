// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_COMMON_TYPEID_H_
#define SEISSOL_SRC_CORE_COMMON_TYPEID_H_

#include <cstdint>
#include <mutex>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

namespace seissol {

/**
 * A numbering of types: every type it is asked for gets the next id, beginning with 0, so that the
 * ids can index tables over the types.
 */
class TypeNumbering {
  public:
  /// The id of `type`; a type it has not seen yet gets the next one.
  std::uint32_t idOf(const std::type_info& type);

  /// The number of ids given out so far.
  [[nodiscard]] std::uint32_t size() const;

  private:
  mutable std::mutex mutex_;
  std::unordered_map<std::type_index, std::uint32_t> ids_;
};

/// The purposes types are numbered for, each on its own, so that each numbering stays dense.
enum class TypeDomain : std::uint8_t {
  /// the variables, buckets and scratchpads of the storages (see initializer::Storage::add)
  Storage,
  /// the variables of the recorded batches (see recording::BatchTable)
  Batch,
  /// the kernels the batches are recorded for (see recording::KernelName)
  Kernel,
};

/**
 * The numbering of the process for `domain`. It numbers types, not data: a type has the same id in
 * every library of the process and for every simulation in it, like its `typeid`.
 */
TypeNumbering& typeNumbering(TypeDomain domain);

/**
 * The id of the type `T` in `Domain`. Every id a program names is given out during the static
 * initialization, so that typeCount() covers them all from then on.
 */
template <TypeDomain Domain, typename T>
inline const std::uint32_t TypeId = typeNumbering(Domain).idOf(typeid(T));

/// The number of ids given out in `domain` so far.
inline std::uint32_t typeCount(TypeDomain domain) { return typeNumbering(domain).size(); }

} // namespace seissol

#endif // SEISSOL_SRC_CORE_COMMON_TYPEID_H_
