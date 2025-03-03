// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LTSSYNC_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LTSSYNC_H_

#include <cstddef>
#include <cstring>
#include <type_traits>

#include "Initializer/MemoryManager.h"
#include "LTSTree.h"
#include "Lut.h"

namespace seissol {
class SeisSol;
namespace initializer {

/*
Assigns the given value to the target object, initializing the memory in the process.

NOTE: std::copy (or the likes) do not work here, since they do not initialize the _vptr for virtual
function calls (rather, they leave it undefined), since they do merely assign `value` to `target`.
*/

template <typename T>
void initAssign(T& target, const T& value) {
  if constexpr (std::is_trivially_copyable_v<T>) {
    // if the object is trivially copyable, we may just memcpy it (it's safe to do that in this
    // case).
    std::memcpy(&target, &value, sizeof(T));
  } else {
    // otherwise, call the class/struct initializer.
    // problem: we may have an array here... So we unwrap it.
    if constexpr (std::is_array_v<T>) {
      // unwrap array, dimension by dimension...
      // example: T[N][M] yields SubT=T[M]
      using SubT = std::remove_extent_t<T>;
      auto subExtent = std::extent_v<T>;

      // for now, init element-wise... (TODO(David): we could look for something faster here, in
      // case it should ever matter)
      for (size_t i = 0; i < subExtent; ++i) {
        initAssign<SubT>(target[i], value[i]);
      }
    } else {
      // now call new here.
      new (&target) T(value);
    }
  }
  // (these two methods cannot be combined, unless we have some way for C-style arrays, i.e. S[N]
  // for <typename S, size_t N>, to use a copy constructor as well)
}

template <typename T>
void synchronizeLTSTreeDuplicates(const seissol::initializer::Variable<T>& handle,
                                  seissol::initializer::MemoryManager& memoryManager) {
  const auto& meshToLts = memoryManager.getLtsLut()->getMeshToLtsLut(handle.mask);
  const unsigned* duplicatedMeshIds = memoryManager.getLtsLut()->getDuplicatedMeshIds(handle.mask);
  const unsigned numberOfDuplicatedMeshIds =
      memoryManager.getLtsLut()->getNumberOfDuplicatedMeshIds(handle.mask);
  T* var = memoryManager.getLtsTree()->var(handle);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned dupMeshId = 0; dupMeshId < numberOfDuplicatedMeshIds; ++dupMeshId) {
    const unsigned meshId = duplicatedMeshIds[dupMeshId];
    const T& ref = var[meshToLts[0][meshId]];
    for (unsigned dup = 1; dup < seissol::initializer::Lut::MaxDuplicates &&
                           meshToLts[dup][meshId] != std::numeric_limits<unsigned>::max();
         ++dup) {

      T& target = var[meshToLts[dup][meshId]];
      // copy data on a byte-wise level (we need to initialize memory here as well)
      initAssign(target, ref);
    }
  }
}

} // namespace initializer
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_TREE_LTSSYNC_H_
