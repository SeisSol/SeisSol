// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FaceCluster.h"

#include "Core/TimeStepping/Actor/ActorState.h"

namespace seissol::solver {

DataReadiness FaceCluster::dataReadiness() const { return DataReadiness::AfterCorrection; }

void FaceCluster::correct() { interact(stepParams()); }

} // namespace seissol::solver
