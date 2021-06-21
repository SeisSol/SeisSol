#ifndef SEISSOL_CONDITIONALTABLE_HPP
#define SEISSOL_CONDITIONALTABLE_HPP

#include "BatchTable.hpp"
#include "Condition.hpp"
#include "ConditionalKey.hpp"
#include "EncodedConstants.hpp"

namespace seissol {
namespace initializers {
namespace recording {
using ConditionalBatchTableT =
    std::unordered_map<ConditionalKey, BatchTable, ConditionalHash<ConditionalKey>>;
}
} // namespace initializers
} // namespace seissol

#endif // SEISSOL_CONDITIONALTABLE_HPP
