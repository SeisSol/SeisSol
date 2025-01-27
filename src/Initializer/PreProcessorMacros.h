// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PREPROCESSORMACROS_H_
#define SEISSOL_SRC_INITIALIZER_PREPROCESSORMACROS_H_

#include "Monitoring/Instrumentation.h"
#include <cstddef>

#include "Common/Constants.h"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// for now, turn all #defines which are left into constexprs

// workaround for old NVHPC versions (the output would cause errors there)
#ifdef __NVCOMPILER
// we'll leave the comment in the next line in for now, until a NVHPC version is fixed
#if (__NVCOMPILER_MAJOR__ > 24) || (__NVCOMPILER_MAJOR__ == 24 && __NVCOMPILER_MINOR__ >= 7)
#define NVHPC_AVOID_OMP 0
#else
#define NVHPC_AVOID_OMP 1
#endif
#else
#define NVHPC_AVOID_OMP 0
#endif

#endif // SEISSOL_SRC_INITIALIZER_PREPROCESSORMACROS_H_
