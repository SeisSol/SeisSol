// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_CLOCKAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_CLOCKAUX_H_

namespace seissol::kernels::device::aux::clock {

/**
 * Sets a clock (see solver::ClusterClock) to `time`, in the order of the stream: the value on the
 * device, and the one on the host, which the device writes itself (`host` is its address on the
 * device), so that no copy engine is involved.
 */
void set(double* device, double* host, double time, void* stream);

/**
 * Advances a clock by `timeStepSize`, in the order of the stream: both values become the sum of the
 * one on the device and the step.
 */
void advance(double* device, double* host, double timeStepSize, void* stream);

} // namespace seissol::kernels::device::aux::clock

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_CLOCKAUX_H_
