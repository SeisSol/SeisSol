// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BUCKETS_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BUCKETS_H_

#include <Common/Real.h>
#include <Initializer/TimeStepping/Halo.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/LTSTree.h>
#include <Solver/TimeStepping/HaloCommunication.h>

namespace seissol::initializer::internal {

solver::HaloCommunication bucketsAndCommunication(LTS::Storage& storage, const MeshLayout& layout);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BUCKETS_H_
