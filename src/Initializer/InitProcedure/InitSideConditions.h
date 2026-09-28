// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITSIDECONDITIONS_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITSIDECONDITIONS_H_

#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Typedefs.h"

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initSideConditions(seissol::SeisSol& seissolInstance);

/**
 * The acoustic travelling wave with ITM as the parameter file sets it up. The wave passes the
 * mirror only if ITM is enabled. Exposed for testing.
 */
AcousticTravellingWaveParametersITM getAcousticTravellingWaveITMInformation(
    const parameters::InitializationParameters& initConditionParams,
    const parameters::ITMParameters& itmParams);
} // namespace seissol::initializer::initprocedure

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITSIDECONDITIONS_H_
