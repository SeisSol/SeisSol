// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
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

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 2
#define DMO_STREAM(IN, OUT)                                                                        \
  uint64_t v1 = *(reinterpret_cast<const uint64_t*>(IN));                                          \
  uint64_t v2 = *(reinterpret_cast<const uint64_t*>(IN) + 1);                                      \
  asm volatile(                                                                                    \
      "stnp %[IN1],%[IN2],[%[OUTaddr]]" ::[IN1] "r"(v1), [IN2] "r"(v2), [OUTaddr] "r"(OUT)         \
      :);

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 4
#define DMO_STREAM(IN, OUT)                                                                        \
  uint64_t v1 = *(reinterpret_cast<const uint64_t*>(IN));                                          \
  uint64_t v2 = *(reinterpret_cast<const uint64_t*>(IN) + 1);                                      \
  asm volatile(                                                                                    \
      "stnp %[IN1],%[IN2],[%[OUTaddr]]" ::[IN1] "r"(v1), [IN2] "r"(v2), [OUTaddr] "r"(OUT)         \
      :);

#else
#error no precision was defined
#endif

#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAARCH64_H_
