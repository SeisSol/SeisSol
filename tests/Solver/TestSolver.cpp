// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>
#include <doctest/trompeloeil.hpp>

#include "Estimator.t.h"
#include "TimeStepping/Actor/AbstractTimeCluster.t.h"
#include "TimeStepping/Actor/ActorState.t.h"
#include "TimeStepping/Actor/ConcurrentClusters.t.h"
#include "TimeStepping/Actor/StepParams.t.h"
#include "TimeStepping/Compute/ClusterClock.t.h"
#include "TimeStepping/Halo/GhostCluster.t.h"
#include "TimeStepping/Halo/Stream/ExchangeScheduler.t.h"
#include "TimeStepping/Halo/Stream/StreamOrderedExchange.t.h"
#include "TimeStepping/Plan/TimeSteppingPlan.t.h"
