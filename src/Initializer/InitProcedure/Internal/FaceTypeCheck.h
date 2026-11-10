// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_FACETYPECHECK_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_FACETYPECHECK_H_

#include "Initializer/BasicTypedefs.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Memory/Descriptor/LTS.h"

#include <cstdint>
#include <string>
#include <vector>

namespace seissol::initializer::internal {

//! @brief The bit of @p faceType in the masks of a FaceTypeCensus.
constexpr std::uint32_t faceTypeBit(FaceType faceType) {
  return std::uint32_t{1} << static_cast<std::uint8_t>(faceType);
}

//! @brief The face types the cells of a configuration have, as masks of their faceTypeBit.
struct FaceTypeCensus {
  //! the face types the cells have
  std::uint32_t present{0};
  //! those of them that a cell has which does not qualify for them (see faceTypeCellAdmissible)
  std::uint32_t cellRejected{0};
};

/**
 * @brief What of the face types that the cells of each configuration have (by the id of the
 * configuration) this build cannot handle: one text per face type and configuration, empty if
 * nothing.
 *
 * Needs neither MPI nor the storage, so that a check of the run can tell beforehand from the face
 * types the mesh file gives.
 */
std::vector<std::string> faceTypeProblems(const std::vector<FaceTypeCensus>& census,
                                          parameters::InitializationType scenarioType);

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
