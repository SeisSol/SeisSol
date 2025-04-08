// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_EXECUTOR_H_
#define SEISSOL_SRC_COMMON_EXECUTOR_H_

namespace seissol {

enum class Executor { Host, Device };

constexpr auto executorEnabled(Executor executor) -> bool {
#ifdef ACL_DEVICE
  return true;
#else
  return executor == Executor::Host;
#endif
}

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_EXECUTOR_H_
