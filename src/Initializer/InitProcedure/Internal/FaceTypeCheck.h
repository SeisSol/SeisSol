// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_FACETYPECHECK_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_FACETYPECHECK_H_

#include "Initializer/Parameters/InitializationParameters.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer::internal {

/**
 * Checks every face type occurring in the mesh against the material model and the solver this
 * binary was built for, and aborts with a combined diagnostic if any of them is unavailable.
 *
 * The analytical boundary additionally needs the configured scenario to supply a state at
 * arbitrary times, so it is checked against the scenario registry as well.
 *
 * Without this, an unavailable boundary condition degenerates into a silent no-op: the local
 * kernel falls through its switch and the face contributes nothing, which is indistinguishable
 * from a correct run except in the results.
 */
void checkFaceTypeSupport(LTS::Storage& storage, parameters::InitializationType scenarioType);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_FACETYPECHECK_H_
