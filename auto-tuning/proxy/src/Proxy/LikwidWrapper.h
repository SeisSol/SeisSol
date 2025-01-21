// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_LIKWIDWRAPPER_H_
#define SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_LIKWIDWRAPPER_H_

#ifdef LIKWID_PERFMON
#include <likwid.h>
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

inline void registerMarkers() {
#pragma omp parallel
  {
    LIKWID_MARKER_REGISTER("ader");
    LIKWID_MARKER_REGISTER("localwoader");
    LIKWID_MARKER_REGISTER("local");
    LIKWID_MARKER_REGISTER("neighboring");
  }
}

#endif // SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_LIKWIDWRAPPER_H_
