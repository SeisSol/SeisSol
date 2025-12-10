// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_TOUCH_H_
#define SEISSOL_SRC_KERNELS_TOUCH_H_

#include "Kernels/Precision.h"

namespace seissol::kernels {

template <typename Cfg>
void touchBuffersDerivatives(Real<Cfg>** buffers, Real<Cfg>** derivatives, unsigned numberOfCells);

template <typename RealT>
void fillWithStuff(RealT* buffer, unsigned nValues, bool onDevice);

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TOUCH_H_
