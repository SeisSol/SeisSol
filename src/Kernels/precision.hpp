#ifndef KERNELS_PRECISION_HPP_
#define KERNELS_PRECISION_HPP_

#ifdef USE_MPI
  #include <mpi.h>
#endif
#include <initialization/precision.h>

#ifdef SINGLE_PRECISION
#define REAL_BYTES sizeof(float)
typedef float real;
#endif
#ifdef DOUBLE_PRECISION
#define REAL_BYTES sizeof(double)
typedef double real;
#endif


#ifdef USE_MPI
#ifdef SINGLE_PRECISION
static MPI_Datatype real_mpi = MPI_FLOAT;
#endif
#ifdef DOUBLE_PRECISION
static MPI_Datatype real_mpi = MPI_DOUBLE;
#endif
#endif

#endif
