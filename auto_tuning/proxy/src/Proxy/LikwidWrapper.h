#ifndef LIKWID_WRAPPER_20210526_H
#define LIKWID_WRAPPER_20210526_H

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

#endif // LIKWID_WRAPPER_20210526_H
