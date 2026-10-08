// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_CONFIGBOUNDARYCHECK_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_CONFIGBOUNDARYCHECK_H_

#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer::internal {

/**
 * Checks the faces between cells of different configurations, and aborts with a combined
 * diagnostic if some of them cannot be computed.
 *
 * Across such a face, the time integral of the neighbor is converted into the configuration of
 * the cell, and the neighbor enters the Riemann problem as a material of the cell. That needs the
 * materials of both to pose the Riemann problem in the same material, and the same number of fused
 * simulations. Dynamic rupture faces between configurations are not supported, nor are such faces
 * on GPUs.
 */
void checkConfigBoundaries(LTS::Storage& storage);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_CONFIGBOUNDARYCHECK_H_
