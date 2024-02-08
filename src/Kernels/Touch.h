// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_TOUCH_H_
#define KERNELS_TOUCH_H_

#include <Kernels/precision.hpp>

namespace seissol::kernels {

void touchBuffersDerivatives(real** buffers, real** derivatives, unsigned numberOfCells);
void fillWithStuff(real* buffer, unsigned nValues, bool onDevice);

} // namespace seissol::kernels

#endif // KERNELS_TOUCH_H_
