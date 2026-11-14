// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_INTEGRALAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_INTEGRALAUX_H_

#include <cstddef>

namespace seissol::kernels::device::aux::integral {

/**
 * Adds block `i` of `source` to block `i` of `target`, for the `count` blocks of `size` values, in
 * binary64 and in the order of the stream: the time integrals of the predictor of a step to the
 * time integrals since the start (cf. LTS::Integrals). Both tables of blocks are in device memory.
 *
 * What the generated kernel gpu_integrateOutput does where the device code generator converts
 * between precisions (TensorForge): for the older ones, which do not (see
 * codegen/kernels/outputintegral.py).
 */
void accumulate(const float* const* source,
                double* const* target,
                std::size_t size,
                std::size_t count,
                void* stream);
void accumulate(const double* const* source,
                double* const* target,
                std::size_t size,
                std::size_t count,
                void* stream);

} // namespace seissol::kernels::device::aux::integral

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_INTEGRALAUX_H_
