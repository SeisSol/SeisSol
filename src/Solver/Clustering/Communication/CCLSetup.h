// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLSETUP_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLSETUP_H_
#include <vector>
namespace seissol::solver::clustering::communication {

std::vector<void*> createComms(std::size_t count);

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLSETUP_H_
