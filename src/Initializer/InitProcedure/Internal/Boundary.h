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

#include <cstddef>

namespace seissol::initializer::internal {

void initBoundaryStorage(Boundary::Storage& boundaryStorage, LTS::Storage& storage);
/// `derivedState`: values per face of the state of the derived outputs of the free surface.
void initSurfaceStorage(SurfaceLTS::Storage& surfaceStorage,
                        LTS::Storage& storage,
                        solver::FreeSurfaceIntegrator& freeSurfaceIntegrator,
                        std::size_t derivedState);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_BOUNDARY_H_
