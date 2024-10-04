#define DOCTEST_CONFIG_IMPLEMENT
#include "Parallel/MPI.h"
#include "doctest.h"

int main(int argc, char** argv) {
  seissol::MPI::mpi.init(argc, argv);
  doctest::Context context;

  context.applyCommandLine(argc, argv);

  const int returnValue = context.run();

  seissol::MPI::finalize();

  return returnValue;
}