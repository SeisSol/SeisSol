#ifndef SEISSOL_CONDITIONALTABLE_HPP
#define SEISSOL_CONDITIONALTABLE_HPP

#include "Table.hpp"
#include "Condition.hpp"
#include "ConditionalKey.hpp"
#include "EncodedConstants.hpp"

namespace seissol::initializers::recording {
using ConditionalPointersToRealsTable =
    std::unordered_map<ConditionalKey, PointersToRealsTable, ConditionalHash<ConditionalKey>>;

using DrConditionalPointersToRealsTable =
    std::unordered_map<ConditionalKey, DrPointersToRealsTable, ConditionalHash<ConditionalKey>>;

using ConditionalIndicesTable =
    std::unordered_map<ConditionalKey, IndicesTable, ConditionalHash<ConditionalKey>>;
} // namespace seissol::initializers::recording

#endif // SEISSOL_CONDITIONALTABLE_HPP
