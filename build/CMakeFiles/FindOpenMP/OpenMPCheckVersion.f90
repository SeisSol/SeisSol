
      program omp_ver
      use omp_lib
      implicit none
      integer, parameter :: zero = ichar('0')
      integer, parameter :: ompv = openmp_version
      character, dimension(24), parameter :: ompver_str =&
      (/ 'I', 'N', 'F', 'O', ':', 'O', 'p', 'e', 'n', 'M', 'P', '-',&
         'd', 'a', 't', 'e', '[',&
         char(zero + mod(ompv/100000, 10)),&
         char(zero + mod(ompv/10000, 10)),&
         char(zero + mod(ompv/1000, 10)),&
         char(zero + mod(ompv/100, 10)),&
         char(zero + mod(ompv/10, 10)),&
         char(zero + mod(ompv/1, 10)), ']' /)
      print *, ompver_str
      end program omp_ver
