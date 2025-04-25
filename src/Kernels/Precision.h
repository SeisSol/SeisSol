// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_PRECISION_H_
#define SEISSOL_SRC_KERNELS_PRECISION_H_

#include <Common/Real.h>
#include <Config.h>

// (real should be lower-case)

// NOLINTNEXTLINE
using real = seissol::RealT<seissol::Config::Precision>;

#endif // SEISSOL_SRC_KERNELS_PRECISION_H_
