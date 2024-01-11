#include "Stopwatch.h"

#include <time.h>
#include "Parallel/MPI.h"
#include "utils/logger.h"

#include <string>
#include "Unit.hpp"

namespace seissol {

Stopwatch::Stopwatch() : time(0) {}

/**
 * Reset the stopwatch to zero
 */
void Stopwatch::reset() { time = 0; }

/**
 * starts the time measuring
 */
void Stopwatch::start() { clock_gettime(CLOCK_MONOTONIC, &startTime); }

/**
 * get time measuring
 *
 * @return measured time (until now) in seconds
 */
double Stopwatch::split() {
  struct timespec end;
  clock_gettime(CLOCK_MONOTONIC, &end);

  return seconds(difftime(startTime, end));
}

/**
 * pauses the measuring
 *
 * @return measured time (until now) in seconds
 */
double Stopwatch::pause() {
  struct timespec end;
  clock_gettime(CLOCK_MONOTONIC, &end);

  time += difftime(startTime, end);
  return seconds(time);
}

/**
 * stops time measuring
 *
 * @return measured time in seconds
 */
double Stopwatch::stop() {
  double time = pause();
  reset();
  return time;
}

/**
 * Collective operation, printing avg, min and max time
 */
void Stopwatch::printTime(const char* text, MPI_Comm comm) const {
  int rank = 0;
  double avg = seconds(time);

#ifdef USE_MPI
  double min = seconds(time);
  double max = seconds(time);

  if (comm == MPI_COMM_NULL)
    comm = seissol::MPI::mpi.comm();

  MPI_Comm_rank(comm, &rank);

  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    int size;
    MPI_Comm_size(comm, &size);
    avg /= size;
  } else {
    MPI_Reduce(&avg, 0L, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&min, 0L, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&max, 0L, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  }
#endif // USE_MPI

  logInfo(rank) << text << UnitTime.formatTime(avg).c_str()
#ifdef USE_MPI
                << "(min:" << utils::nospace << UnitTime.formatTime(min).c_str()
                << ", max: " << UnitTime.formatTime(max).c_str() << ')'
#endif // USE_MPI
      ;
}

} // namespace seissol
