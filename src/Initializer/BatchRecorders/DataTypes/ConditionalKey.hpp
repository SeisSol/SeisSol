#ifndef SEISSOL_CONDITIONALKEY_HPP
#define SEISSOL_CONDITIONALKEY_HPP

#include <functional>
#include <limits>
#include <utility>

namespace seissol::initializers::recording {
struct ConditionalKey {
  ConditionalKey(size_t kernel,
                 size_t type = std::numeric_limits<size_t>::max(),
                 size_t face = std::numeric_limits<size_t>::max(),
                 size_t relation = std::numeric_limits<size_t>::max())
      : kernelId(kernel), typeId(type), faceId(face), faceRelationId(relation){};

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
  std::size_t operator()(ConditionalKey const& key) const {
    std::size_t result = 0;
    hashCombine(result, key.kernelId);
    hashCombine(result, key.typeId);
    hashCombine(result, key.faceId);
    hashCombine(result, key.faceRelationId);
    return result;
  }
};
} // namespace seissol::initializers::recording

#endif // SEISSOL_CONDITIONALKEY_HPP
