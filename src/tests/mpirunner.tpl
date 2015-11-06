#ifdef USE_MPI
#include <mpi.h>
#endif
#include <cxxtest/ErrorPrinter.h>

int main(int argc, char** argv)
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif
    
    int status = CxxTest::ErrorPrinter().run();

#ifdef USE_MPI
    MPI_Finalize();
#endif
    
    return status;
}

<CxxTest world>
