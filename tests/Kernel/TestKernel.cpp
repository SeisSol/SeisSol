// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common.t.h"
#include "FlopCounting.t.h"
#include "Plasticity.t.h"
#include "PointSourceCluster.t.h"

#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "NodalBoundary.t.h"
#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#ifdef SEISSOL_KERNELS_STP
#include "STP.t.h"
#endif // SEISSOL_KERNELS_STP
