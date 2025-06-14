// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAVX_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAVX_H_

#include <immintrin.h>

#define DMO_INCREMENT64 4
#define DMO_STREAM64(IN, OUT) _mm256_stream_pd(OUT, _mm256_load_pd(IN));

#define DMO_INCREMENT32 8
#define DMO_STREAM32(IN, OUT) _mm256_stream_ps(OUT, _mm256_load_ps(IN));

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAVX_H_
