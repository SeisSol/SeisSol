// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITIONALKEY_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITIONALKEY_H_

#include <functional>
#include <limits>
#include <utility>

namespace seissol::initializer::recording {
struct ConditionalKey {
  ConditionalKey(size_t kernel,
                 size_t type = std::numeric_limits<size_t>::max(),
                 size_t face = std::numeric_limits<size_t>::max(),
                 size_t relation = std::numeric_limits<size_t>::max())
      : kernelId(kernel), typeId(type), faceId(face), faceRelationId(relation) {};

  bool operator==(const ConditionalKey& other) const {
    return ((kernelId == other.kernelId) && (typeId == other.typeId) && (faceId == other.faceId) &&
            (faceRelationId == other.faceRelationId));
  }

  size_t kernelId;
  size_t typeId;
  size_t faceId;
  size_t faceRelationId;
};

template <class T>
inline void hashCombine(std::size_t& seed, const T& value) {
  std::hash<T> hasher;
  seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <class T>
struct ConditionalHash;

template <>
struct ConditionalHash<ConditionalKey> {
  std::size_t operator()(const ConditionalKey& key) const {
    std::size_t result = 0;
    hashCombine(result, key.kernelId);
    hashCombine(result, key.typeId);
    hashCombine(result, key.faceId);
    hashCombine(result, key.faceRelationId);
    return result;
  }
};
} // namespace seissol::initializer::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITIONALKEY_H_
