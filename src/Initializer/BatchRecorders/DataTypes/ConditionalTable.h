// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITIONALTABLE_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITIONALTABLE_H_

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

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITIONALTABLE_H_
