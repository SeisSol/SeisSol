// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_COMMON_ALIGNEDALLOCATOR_H_
#define SEISSOL_SRC_CORE_COMMON_ALIGNEDALLOCATOR_H_

#include <cstddef>
#include <new>

namespace seissol {

/**
 * Allocates memory of the host aligned to `AlignmentV` bytes (e.g. to `Alignment`, for the vector
 * loads of the generated kernels), with the aligned new of C++17; for the containers of the
 * standard library, e.g. `std::vector<double, AlignedAllocator<double, Alignment>>`. Other kinds of
 * memory (high-bandwidth memory, the device) come from the memory of the core (see
 * memory::MemkindArray).
 */
template <typename T, std::size_t AlignmentV>
class AlignedAllocator {
  public:
  // the names the standard library expects of an allocator
  // NOLINTBEGIN(readability-identifier-naming)
  using value_type = T;

  // (the default of std::allocator_traits only rebinds type parameters)
  template <typename U>
  struct rebind {
    using other = AlignedAllocator<U, AlignmentV>;
  };
  // NOLINTEND(readability-identifier-naming)

  AlignedAllocator() noexcept = default;

  template <typename U>
  explicit AlignedAllocator(const AlignedAllocator<U, AlignmentV>& /*other*/) noexcept {}

  [[nodiscard]] T* allocate(std::size_t count) {
    return static_cast<T*>(::operator new(count * sizeof(T), std::align_val_t(AlignmentV)));
  }

  void deallocate(T* pointer, std::size_t /*count*/) noexcept {
    ::operator delete(pointer, std::align_val_t(AlignmentV));
  }

  friend bool operator==(const AlignedAllocator& /*lhs*/,
                         const AlignedAllocator& /*rhs*/) noexcept {
    return true;
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_CORE_COMMON_ALIGNEDALLOCATOR_H_
