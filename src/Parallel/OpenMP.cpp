#include "OpenMP.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seissol {

bool OpenMP::enabled() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

std::size_t OpenMP::threadId() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

std::size_t OpenMP::threadCount() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

} // namespace seissol
