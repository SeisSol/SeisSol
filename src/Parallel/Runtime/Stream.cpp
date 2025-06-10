// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Stream.h"
namespace seissol::parallel::runtime {

std::mutex StreamRuntime::mutexCPU = std::mutex();

} // namespace seissol::parallel::runtime
