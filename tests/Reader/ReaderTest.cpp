// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "FSRMReader.t.h"

#ifdef USE_NETCDF
#include "NRFReader.t.h"
#endif

#ifdef USE_EASI
#include "EasiReader.t.h"
#include "EasiWalker.t.h"
#endif

#include "Datafield/AsagiLite.t.h"
#include "Datafield/GridStore.t.h"
#include "Datafield/GridUpdateModule.t.h"
#include "Datafield/Interpolation.t.h"
#include "Datafield/TimeWindow.t.h"
#include "LuaFrontend.t.h"
#include "LuaTracer.t.h"
