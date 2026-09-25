// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>
#include <doctest/trompeloeil.hpp>

#include "Estimator.t.h"
#include "TimeStepping/AbstractTimeCluster.t.h"
#include "TimeStepping/ActorState.t.h"
#include "TimeStepping/ClusterClock.t.h"
#include "TimeStepping/ConcurrentClusters.t.h"
#include "TimeStepping/ExchangeScheduler.t.h"
#include "TimeStepping/GhostCluster.t.h"
#include "TimeStepping/StepParams.t.h"
#include "TimeStepping/StreamOrderedExchange.t.h"
#include "TimeStepping/TimeSteppingPlan.t.h"
