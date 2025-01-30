// SPDX-FileCopyrightText: 2014-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MONITORING_INSTRUMENTATION_H_
#define SEISSOL_SRC_MONITORING_INSTRUMENTATION_H_

// Manual instrumentation for Scalasca with epik.
// Override function calls if not compiled with EPIK.
#if defined(EPIK)
#include "epik_user.h"
#else
#define EPIK_FUNC_REG(str)
#define EPIK_FUNC_START()
#define EPIK_FUNC_END()
#define EPIK_USER_REG(id, str)
#define EPIK_USER_START(id)
#define EPIK_USER_END(id)
#define EPIK_TRACER(str)
#endif

// empty score-p definitions without usage
#if defined(SCOREP_USER_ENABLE)
#include <scorep/SCOREP_User.h>
#else
#define SCOREP_USER_REGION(name, type)
#define SCOREP_USER_REGION_DEFINE(handle)
#define SCOREP_USER_OA_PHASE_BEGIN(handle, name, type)
#define SCOREP_USER_OA_PHASE_END(handle)
#define SCOREP_USER_REGION_BEGIN(handle, name, type)
#define SCOREP_USER_REGION_INIT(handle, name, type)
#define SCOREP_USER_REGION_END(handle)
#define SCOREP_USER_REGION_ENTER(handle)
#define SCOREP_USER_FUNC_BEGIN()
#define SCOREP_USER_FUNC_END()
#define SCOREP_GLOBAL_REGION_DEFINE(handle)
#define SCOREP_GLOBAL_REGION_EXTERNAL(handle)
#define SCOREP_USER_PARAMETER_DEFINE(handle)
#define SCOREP_USER_PARAMETER_INT64(name, value)
#define SCOREP_USER_PARAMETER_UINT64(name, value)
#define SCOREP_USER_PARAMETER_STRING(name, value)
#define SCOREP_USER_METRIC_GLOBAL(metricHandle)
#define SCOREP_USER_METRIC_EXTERNAL(metricHandle)
#define SCOREP_USER_METRIC_LOCAL(metricHandle)
#define SCOREP_USER_METRIC_INIT(metricHandle, name, unit, type, context)
#define SCOREP_USER_METRIC_INT64(metricHandle, value)
#define SCOREP_USER_METRIC_UINT64(metricHandle, value)
#define SCOREP_USER_METRIC_DOUBLE(metricHandle, value)
#define SCOREP_RECORDING_ON()
#define SCOREP_RECORDING_OFF()
#define SCOREP_RECORDING_IS_ON() 0
#endif

// likwid
#if defined(LIKWID_PERFMON)
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

#endif // SEISSOL_SRC_MONITORING_INSTRUMENTATION_H_
