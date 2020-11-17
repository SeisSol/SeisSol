#ifndef SEISSOL_CONDITIONALTABLE_HPP
#define SEISSOL_CONDITIONALTABLE_HPP

#include "Initializer/recording/BatchTable.h"
#include "Initializer/recording/Condition.h"
#include "Initializer/recording/ConditionalKey.h"
#include "Initializer/recording/EncodedConstants.h"

using ConditionalTableT =
    std::unordered_map<ConditionalKey, BatchTable, ConditionalHash<ConditionalKey>>;

#endif // SEISSOL_CONDITIONALTABLE_HPP
