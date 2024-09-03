#ifndef SEISSOL_CONDITIONALTABLE_HPP
#define SEISSOL_CONDITIONALTABLE_HPP

#include "Condition.h"
#include "ConditionalKey.h"
#include "EncodedConstants.h"
#include "Table.h"

namespace seissol::initializer::recording {
using ConditionalPointersToRealsTable =
    std::unordered_map<ConditionalKey, PointersToRealsTable, ConditionalHash<ConditionalKey>>;

using DrConditionalPointersToRealsTable =
    std::unordered_map<ConditionalKey, DrPointersToRealsTable, ConditionalHash<ConditionalKey>>;

using ConditionalMaterialTable =
    std::unordered_map<ConditionalKey, MaterialTable, ConditionalHash<ConditionalKey>>;

using ConditionalIndicesTable =
    std::unordered_map<ConditionalKey, IndicesTable, ConditionalHash<ConditionalKey>>;
} // namespace seissol::initializer::recording

#endif // SEISSOL_CONDITIONALTABLE_HPP
