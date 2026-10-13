// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_TOUCH_H_
#define SEISSOL_SRC_KERNELS_TOUCH_H_

#include "Common/Real.h"

namespace seissol::kernels {

/// Zeroes the time integrated buffers and the derivatives of cells of the configuration `Cfg`.
template <typename Cfg>
void touchBuffersDerivatives(Real<Cfg>** buffers, Real<Cfg>** derivatives, unsigned numberOfCells);

/// Fills `buffer` with arbitrary finite values that are not zero.
template <typename RealT>
void fillWithStuff(RealT* buffer, unsigned nValues, bool onDevice);

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TOUCH_H_
