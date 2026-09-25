// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FaceCluster.h"

#include "Common/Executor.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"

namespace seissol::time_stepping {

FaceCluster::FaceCluster(double maxTimeStepSize, long timeStepRate, Executor executor)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate, executor) {}

DataReadiness FaceCluster::dataReadiness() const { return DataReadiness::AfterCorrection; }

void FaceCluster::correct() { interact(stepParams()); }

} // namespace seissol::time_stepping
