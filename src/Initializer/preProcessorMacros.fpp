#if 0
!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!! @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
!!
!! @section LICENSE
!! Copyright (c) 2012-2015, SeisSol Group
!! All rights reserved.
!! 
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!! 
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!! 
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!! 
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived from this
!!    software without specific prior written permission.
!! 
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!!
!! @section DESCRIPTION
!! Defines preprocessor variables and macros used in SeisSol.
#endif

#ifndef PRE_PROCESSOR_MACROS_FPP
#define PRE_PROCESSOR_MACROS_FPP

#include "Monitoring/instrumentation.fpp"

#if 0
!
! Fortran Logging
!
#endif
#if 0
! Unfortunately on BG is not that easy to know if we are compiling C or Fortran.
! In this case, we check for the flag (BG) that specified in the compilation string and
! if __bg__ is defined (the latter is only defined on XL C anc XL C++ compilers, not F90) 
#endif
#if !defined(__STDC__) || (defined(BG) && !defined(__bg__))

#ifndef LOGLEVEL0
#define LOGLEVEL0 LOGLEVEL
#endif

#define FORTRAN_STDERR 0
#define FORTRAN_STDOUT 6

#define FORTRAN_LINE_SIZE 1500

#define logDebug(format)   if (LOGLEVEL .ge. 3) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Debug   |'; if (LOGLEVEL .ge. 3) write(FORTRAN_STDERR, format)
#define logInfo(format)    if (LOGLEVEL .ge. 2) write(FORTRAN_STDOUT, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Info    |'; if (LOGLEVEL .ge. 2) write(FORTRAN_STDOUT, format)
#define logWarning(format) if (LOGLEVEL .ge. 1) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Warning |'; if (LOGLEVEL .ge. 1) write(FORTRAN_STDERR, format)
#define logError(format)   if (LOGLEVEL .ge. 0) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Error   |'; if (LOGLEVEL .ge. 0) write(FORTRAN_STDERR, format)

#define logDebug0(format)   if (LOGLEVEL0 .ge. 3 .and. myrank .eq. 0) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', 0, ' | Debug   |'; if (LOGLEVEL0 .ge. 3 .and. myrank .eq. 0) write(FORTRAN_STDERR, format)
#define logInfo0(format)    if (LOGLEVEL0 .ge. 2 .and. myrank .eq. 0) write(FORTRAN_STDOUT, '(A, I8, A)', advance='no') 'Rank: ', 0, ' | Info    |'; if (LOGLEVEL0 .ge. 2 .and. myrank .eq. 0) write(FORTRAN_STDOUT, format)
#define logWarning0(format) if (LOGLEVEL0 .ge. 1 .and. myrank .eq. 0) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', 0, ' | Warning |'; if (LOGLEVEL0 .ge. 1 .and. myrank .eq. 0) write(FORTRAN_STDERR, format)

#endif

#if 0
!
! Generated Kernels
!
#endif


#if 0
! preprocessor concatenation
#endif
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#if 0
! zero tolerance constant
#endif
#ifndef ZEROTOLERANCE
#define ZEROTOLERANCE 1e-10
#endif

#ifndef CONVERGENCE_ORDER
#error  Preprocessor flag CONVERGENCE_ORDER not set.
#endif

#if 0
// define number of basis functions relative to the convergence order
#endif
#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_BASIS_FUNCTIONS 4

#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_BASIS_FUNCTIONS 10

#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_BASIS_FUNCTIONS 20

#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_BASIS_FUNCTIONS 35

#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_BASIS_FUNCTIONS 56

#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_BASIS_FUNCTIONS 84

#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_BASIS_FUNCTIONS 120

#else
#error Preprocessor flag CONVERGENCE_ORDER is not in {2, 3, 4, 5, 6, 7, 8}.
#endif


#if 0
! page sizes
! heap: perfect use of transparent huge pages
#endif
#define PAGESIZE_HEAP 2097152
#define PAGESIZE_STACK 4096

#if 0
! fortran specific variables
#endif
#if !defined(__STDC__) || (defined(BG) && !defined(__bg__))
! position of the stderr stream
#ifndef STDERR
#define STDERR 0
#endif

#if (NUMBER_OF_BASIS_FUNCTIONS !=  1) && (NUMBER_OF_BASIS_FUNCTIONS !=   4)  &&\
    (NUMBER_OF_BASIS_FUNCTIONS != 10) && (NUMBER_OF_BASIS_FUNCTIONS !=  20)  &&\
    (NUMBER_OF_BASIS_FUNCTIONS != 35) && (NUMBER_OF_BASIS_FUNCTIONS !=  56)  &&\
    (NUMBER_OF_BASIS_FUNCTIONS != 84) && (NUMBER_OF_BASIS_FUNCTIONS != 120)
#error Preprocessor flag NUMBER_OF_BASIS_FUNCTIONS is not in {1, 4, 10, 20, 35, 56, 84, 120}.
#endif

#endif


#endif

#if defined(REAL_SIZE)
#define REAL_TYPE real(kind=REAL_SIZE)
#else
#error Unknown real size.
#endif

#define ALLOW_POSSILBE_ZERO_LENGTH_ARRAY(X) X == 0 ? 1 : X
