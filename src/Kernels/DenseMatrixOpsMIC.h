// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSMIC_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSMIC_H_

#include <immintrin.h>

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 8
#define DMO_STREAM(IN, OUT) _mm512_storenrngo_pd(OUT, _mm512_load_pd(IN));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 16
#define DMO_STREAM(IN, OUT) _mm512_storenrngo_ps(OUT, _mm512_load_ps(IN));

#else
#error no precision was defined
#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSMIC_H_

