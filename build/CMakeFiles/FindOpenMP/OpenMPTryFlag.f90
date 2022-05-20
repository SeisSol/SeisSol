
      program test
      implicit none
      include 'omp_lib.h'
  !$  integer :: n
      n = omp_get_num_threads()
      end program test
  
