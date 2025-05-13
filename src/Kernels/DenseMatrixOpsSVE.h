// SPDX-FileCopyrightText: 2023 SeisSol Group
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

#define DMO_INCREMENT64 svcntd()
#define DMO_STREAM64(IN, OUT) svstnt1_f64(svptrue_b64(), OUT, svld1_f64(svptrue_b64(), IN));

#define DMO_INCREMENT32 svcntw()
#define DMO_STREAM32(IN, OUT) svstnt1_f32(svptrue_b32(), OUT, svld1_f32(svptrue_b32(), IN));

#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSSVE_H_
