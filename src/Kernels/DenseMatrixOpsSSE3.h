// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSSE3_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSSE3_H_

#include <immintrin.h>

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 2
#define DMO_STREAM(IN, OUT) _mm_stream_pd(OUT, _mm_load_pd(IN));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 4
#define DMO_STREAM(IN, OUT) _mm_stream_ps(OUT, _mm_load_ps(IN));

#else
#error no precision was defined
#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSSE3_H_
