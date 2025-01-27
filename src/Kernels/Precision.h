// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_PRECISION_H_
#define SEISSOL_SRC_KERNELS_PRECISION_H_

#ifdef USE_MPI
#include <mpi.h>
#endif

#if REAL_SIZE == 8
#define DOUBLE_PRECISION
#elif REAL_SIZE == 4
#define SINGLE_PRECISION
#else
#error REAL_SIZE not supported.
#endif

// (real should be lower-case)

#ifdef SINGLE_PRECISION
// NOLINTNEXTLINE
using real = float;
#endif
#ifdef DOUBLE_PRECISION
// NOLINTNEXTLINE
using real = double;
#endif

#ifdef USE_MPI
#ifdef SINGLE_PRECISION
#define MPI_C_REAL MPI_FLOAT
#endif
#ifdef DOUBLE_PRECISION
#define MPI_C_REAL MPI_DOUBLE
#endif
#endif

#ifdef USE_HDF
#ifdef SINGLE_PRECISION
#define HDF_C_REAL H5T_NATIVE_FLOAT
#endif
#ifdef DOUBLE_PRECISION
#define HDF_C_REAL H5T_NATIVE_DOUBLE
#endif
#endif

#endif // SEISSOL_SRC_KERNELS_PRECISION_H_
