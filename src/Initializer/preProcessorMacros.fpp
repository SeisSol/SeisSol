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
#include "generated_code/initialization/precision.h"

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

#ifdef GENERATEDKERNELS

#if 0
! preprocessor concatenation
#endif
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define CONCAT_HELPER_4(a,b,c,d) a ## b ## c ## d
#define CONCAT_4(a,b,c,d) CONCAT_HELPER_4(a,b,c,d)
#define CONCAT_HELPER_6(a,b,c,d,e,f) a ## b ## c ## d ## e ## f
#define CONCAT_6(a,b,c,d,e,f) CONCAT_HELPER_6(a,b,c,d, e, f)

#if defined(__STDC__)
#include <initialization/bind.h>
#endif

#if 0
! zero tolerance constant
#endif
#ifndef ZEROTOLERANCE
#define ZEROTOLERANCE 1e-10
#endif

#if 0
! check if the necessary precompiler macros are defined
#endif
#ifndef NUMBER_OF_QUANTITIES
#error Preprocessor flag NUMBER_OF_QUANTITIES not set.
#endif

#ifndef CONVERGENCE_ORDER
#error  Preprocessor flag ONVERGENCE_ORDER not set.
#endif


#if 0
! check valid preprocessor flags.
! * The number of variables depends on the degree of attenuation.
!   Whereas 9 variables correspond to the elastic wave equations.
! * The number of basis functions depends on the order of the discontinuous Galerkin method.
!   Therfore valid values are given by the formula (O)*(O+1)*(O+2)/6, where O is the order of method.
#endif
#if NUMBER_OF_QUANTITIES != 9
#error Preprocessor flag NUMBER_OF_QUANTITIES is not equal 9 (elastic wave equations).
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
! memory subsystem specifications
#endif

#ifdef USE_MEMKIND
#define MEMKIND_GLOBAL 1
#define MEMKIND_CONSTANT 1
#define MEMKIND_DOFS 1
#define MEMKIND_TIMEDOFS 1
#else
#define MEMKIND_GLOBAL 0
#define MEMKIND_CONSTANT 0
#define MEMKIND_DOFS 0
#define MEMKIND_TIMEDOFS 0
#endif

#if 0
! aligned number of basis functions
#endif

#if ALIGNMENT == 16 && defined(DOUBLE_PRECISION)

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 4
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 6
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 10
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 16
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 20
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 36
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 36
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 72
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 128
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 84
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 212
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 332
#endif

#elif ( ALIGNMENT == 32 && defined(DOUBLE_PRECISION) ) || ( ALIGNMENT == 16 && defined(SINGLE_PRECISION) )

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 4
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 8
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 12
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 20
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 20
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 40
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 36
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 76
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 132
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 84
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 216
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 336
#endif

#elif ( ALIGNMENT == 64 && defined(DOUBLE_PRECISION) ) || ( ALIGNMENT == 32 && defined(SINGLE_PRECISION) )

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 8
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 16
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 32
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 24
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 56
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 40
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 96
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 152
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 88
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 240
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 360
#endif

#elif ALIGNMENT == 64 && defined(SINGLE_PRECISION)

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 32
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 48
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 32
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 80
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 48
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 128
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 64
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 192
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 96
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 288
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 128
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS 416
#endif

#else

#error alignment-precision combination not implemented.

#endif

#if 0
! number of (aligned) dofs.
#endif
#define NUMBER_OF_DOFS         (NUMBER_OF_BASIS_FUNCTIONS             * NUMBER_OF_QUANTITIES)
#define NUMBER_OF_ALIGNED_DOFS (NUMBER_OF_ALIGNED_BASIS_FUNCTIONS     * NUMBER_OF_QUANTITIES)
#define NUMBER_OF_ALIGNED_DERS (NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS * NUMBER_OF_QUANTITIES)

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

#endif
