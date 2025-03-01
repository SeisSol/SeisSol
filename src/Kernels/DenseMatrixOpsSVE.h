// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSVE_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSVE_H_

#include <Kernels/Precision.h>

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT svcntd()
#define DMO_STREAM(IN, OUT) svstnt1_f64(svptrue_b64(), OUT, svld1_f64(svptrue_b64(), IN));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT svcntw()
#define DMO_STREAM(IN, OUT) svstnt1_f32(svptrue_b32(), OUT, svld1_f32(svptrue_b32(), IN));

#else
#error no precision was defined
#endif

#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSVE_H_
