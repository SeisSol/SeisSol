// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_DATATYPES_CONDITIONALTABLE_H_
#define SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_DATATYPES_CONDITIONALTABLE_H_

#include "Condition.h"
#include "ConditionalKey.h"
#include "EncodedConstants.h"
#include "Table.h"

namespace seissol::recording {
/// The batches of a layer, by their keys.
using ConditionalTable =
    std::unordered_map<ConditionalKey, BatchTable, ConditionalHash<ConditionalKey>>;
} // namespace seissol::recording

#endif // SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_DATATYPES_CONDITIONALTABLE_H_
