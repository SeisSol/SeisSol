/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Common kernel-level functions
 **/

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <algorithm>
#include <Initializer/typedefs.hpp>
#include <generated_code/init.h>
#include <generated_code/kernel.h>
#include <cassert>

namespace seissol {
  namespace kernels {
    /**
     * Gets the number of basis functions for the given convergence order.
     *
     * @param i_convergenceOrder convergence order.
     * @return number of basis funcitons.
     **/
    inline unsigned int getNumberOfBasisFunctions( unsigned int i_convergenceOrder = CONVERGENCE_ORDER ) {
      return i_convergenceOrder*(i_convergenceOrder+1)*(i_convergenceOrder+2)/6;
    }

    /**
     * Gets the number of aligned reals.
     *
     * @param i_alignment alignment in bytes.
     * @return aligned number of reals.
     **/
    inline unsigned int getNumberOfAlignedReals( unsigned int i_numberOfReals,
                                                 unsigned int i_alignment = ALIGNMENT ) {
      unsigned int const alignment = i_alignment / sizeof(real);
      return i_numberOfReals + (alignment-(i_numberOfReals % alignment))%alignment;
    }

    /**
     * Get the # of basis functions aligned to the given boundaries.
     *
     * @param i_convergenceOrder convergence order.
     * @param i_alignment alignment in bytes.
     * @return aligned number of basis functions.
     **/
    inline unsigned int getNumberOfAlignedBasisFunctions( unsigned int i_convergenceOrder = CONVERGENCE_ORDER,
                                                   unsigned int i_alignment        = ALIGNMENT ) {
      unsigned int l_numberOfBasisFunctions = getNumberOfBasisFunctions( i_convergenceOrder);
      return getNumberOfAlignedReals( l_numberOfBasisFunctions );
    }

    /**
     * Converts memory aligned degrees of freedom (with zero padding) to unaligned (compressed, without zero padding) storage.
     *
     * @param i_alignedDofs aligned degrees of freedom (zero padding in the basis functions / columns).
     * @param o_unalignedDofs unaligned degrees of freedom.
     **/
    template<typename real_from, typename real_to>
    void convertAlignedDofs( const real_from i_alignedDofs[tensor::Q::size()],
                                   real_to   o_unalignedDofs[tensor::QFortran::size()] ) {
      kernel::copyQToQFortran krnl;
      krnl.Q = i_alignedDofs;
#ifdef MULTIPLE_SIMULATIONS
      krnl.multSimToFirstSim = init::multSimToFirstSim::Values;
#endif

      if (std::is_same<real_from, real_to>::value) {
        krnl.QFortran = reinterpret_cast<real_from*>(o_unalignedDofs);
        krnl.execute();
      } else {
        real_from unalignedDofs[tensor::QFortran::size()];
        krnl.QFortran = unalignedDofs;
        krnl.execute();
        std::copy(unalignedDofs, unalignedDofs + tensor::QFortran::size(), o_unalignedDofs);
      }
    }
  }
}

#endif
