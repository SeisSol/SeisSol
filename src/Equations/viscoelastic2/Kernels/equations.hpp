/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

/* If we enable C++11 at some time, the following defines may be replaced
 * by the following code:
constexpr unsigned numberOfBasisFunctions(unsigned O) {
  return O * (O + 1) * (O + 2) / 6;
}

constexpr unsigned numberOfAlignedBasisFunctions(unsigned O) {
  return (numberOfBasisFunctions(O) * REAL_BYTES + (ALIGNMENT - (numberOfBasisFunctions(O) * REAL_BYTES) % ALIGNMENT) % ALIGNMENT) / REAL_BYTES;
}

constexpr unsigned numberOfAlignedDerBasisFunctions(unsigned O) {
  return (O > 0) ? numberOfAlignedBasisFunctions(O) + numberOfAlignedDerBasisFunctions(O-1) : 0;
}

#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS numberOfAlignedBasisFunctions(CONVERGENCE_ORDER)
#define NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS numberOfAlignedDerBasisFunctions(CONVERGENCE_ORDER)
 
*/
#if 0
! aligned number of basis functions
#endif

#if ALIGNMENT == 16 && defined(DOUBLE_PRECISION)

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 4
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 10
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 20
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 36
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 84
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#endif

#elif ( ALIGNMENT == 32 && defined(DOUBLE_PRECISION) ) || ( ALIGNMENT == 16 && defined(SINGLE_PRECISION) )

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 4
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 12
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 20
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 36
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 84
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#endif

#elif ( ALIGNMENT == 64 && defined(DOUBLE_PRECISION) ) || ( ALIGNMENT == 32 && defined(SINGLE_PRECISION) )

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 8
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 24
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 40
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 88
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#endif

#elif ALIGNMENT == 64 && defined(SINGLE_PRECISION)

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 32
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 48
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 64
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 96
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 128
#endif

#else

#error alignment-precision combination not implemented.

#endif

#define NUMBER_OF_ALIGNED_STRESS_DOFS     (NUMBER_OF_ALIGNED_BASIS_FUNCTIONS     * 6)
