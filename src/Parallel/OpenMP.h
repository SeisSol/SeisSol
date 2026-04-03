// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_OPENMP_H_
#define SEISSOL_SRC_PARALLEL_OPENMP_H_

#include <cstddef>
namespace seissol {

class OpenMP {
  public:
  static bool enabled();
  static std::size_t threadId();
  static std::size_t threadCount();
};

} // namespace seissol
#endif // SEISSOL_SRC_PARALLEL_OPENMP_H_
