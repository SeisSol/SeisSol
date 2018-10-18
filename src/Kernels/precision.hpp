#ifndef KERNELS_PRECISION_HPP_
#define KERNELS_PRECISION_HPP_

#ifdef USE_MPI
  #include <mpi.h>
#endif

#if REAL_SIZE == 8
#  define DOUBLE_PRECISION
#elif REAL_SIZE == 4
#  define SINGLE_PRECISION
#else
#  error REAL_SIZE not supported.
#endif

#ifdef SINGLE_PRECISION
typedef float real;
#endif
#ifdef DOUBLE_PRECISION
typedef double real;
#endif


#ifdef USE_MPI
#ifdef SINGLE_PRECISION
#define MPI_C_REAL MPI_FLOAT
#endif
#ifdef DOUBLE_PRECISION
#define MPI_C_REAL MPI_DOUBLE
#endif
#endif

#endif
