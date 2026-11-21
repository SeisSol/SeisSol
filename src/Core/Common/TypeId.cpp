// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TypeId.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <typeindex>
#include <typeinfo>

namespace seissol {

std::uint32_t TypeNumbering::idOf(const std::type_info& type) {
  const std::scoped_lock lock(mutex_);
  // (std::type_index compares by the name of the type, so the copies of a type_info in several
  // libraries count as one)
  const auto [entry, inserted] =
      ids_.try_emplace(std::type_index(type), static_cast<std::uint32_t>(ids_.size()));
  return entry->second;
}

std::uint32_t TypeNumbering::size() const {
  const std::scoped_lock lock(mutex_);
  return static_cast<std::uint32_t>(ids_.size());
}

TypeNumbering& typeNumbering(TypeDomain domain) {
  // (a function-local static, as the ids are given out during the static initialization)
  static std::array<TypeNumbering, static_cast<std::size_t>(TypeDomain::Kernel) + 1> numberings;
  return numberings[static_cast<std::size_t>(domain)];
}

} // namespace seissol
