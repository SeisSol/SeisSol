// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BOUNDARY_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BOUNDARY_H_

#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Solver/FreeSurfaceIntegrator.h"

namespace seissol::initializer::internal {

void initBoundaryStorage(Boundary::Storage& boundaryStorage, LTS::Storage& storage);
void initSurfaceStorage(SurfaceLTS::Storage& surfaceStorage,
                        LTS::Storage& storage,
                        solver::FreeSurfaceIntegrator& freeSurfaceIntegrator,
                        int refinement);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BOUNDARY_H_
