// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAARCH64_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAARCH64_H_

#include <Kernels/Precision.h>

#ifdef __aarch64__

// cf. https://stackoverflow.com/a/61248308

#define DMO_INCREMENT64 2
#define DMO_STREAM64(IN, OUT)                                                                      \
  uint64_t v1 = *(reinterpret_cast<const uint64_t*>(IN));                                          \
  uint64_t v2 = *(reinterpret_cast<const uint64_t*>(IN) + 1);                                      \
  asm volatile(                                                                                    \
      "stnp %[IN1],%[IN2],[%[OUTaddr]]" ::[IN1] "r"(v1), [IN2] "r"(v2), [OUTaddr] "r"(OUT)         \
      :);

#define DMO_INCREMENT32 4
#define DMO_STREAM32(IN, OUT)                                                                      \
  uint64_t v1 = *(reinterpret_cast<const uint64_t*>(IN));                                          \
  uint64_t v2 = *(reinterpret_cast<const uint64_t*>(IN) + 1);                                      \
  asm volatile(                                                                                    \
      "stnp %[IN1],%[IN2],[%[OUTaddr]]" ::[IN1] "r"(v1), [IN2] "r"(v2), [OUTaddr] "r"(OUT)         \
      :);

#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAARCH64_H_
