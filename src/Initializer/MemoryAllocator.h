/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * Aligned memory allocation.
 **/

#ifndef MEMORYALLOCATOR_H_
#define MEMORYALLOCATOR_H_

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef USE_MEMKIND
#include <hbwmalloc.h>
#endif

#include <utils/logger.h>

namespace seissol {
  class MemoryAllocator;
}

/**
 * Allocates aligned memory for dynamic chunk sizes.
 **/
class seissol::MemoryAllocator {
  private:
    typedef std::pair<unsigned, void*> Address; 
    typedef std::vector<Address>       AddressVector;
  
    //! holds all memory addresses, which point to data arrays and have been returned by mallocs calling functions of the memory allocator.
    AddressVector m_dataMemoryAddresses;

    /**
     * Prints the memory alignment of in terms of relative start and ends in bytes.
     *
     * @param i_memoryAlignment memory alignment.
     **/
    void printMemoryAlignment( std::vector< std::vector<unsigned long long> > i_memoryAlignment );

  public:
    enum Memkind {
      Standard = 0,
      HighBandwidth = 1
    };
  
    MemoryAllocator();
    ~MemoryAllocator();

    /**
     * Allocates a single chunk of memory with the given size and alignment.
     *
     * @param  i_size size of the chunk in byte.
     * @param  i_alignment alignment of the memory chunk in byte.
     * @return pointer, which points to the aligned memory of the given size.
     *
     * @todo Allow this function to be called directly (update m_dataMemoryAddresses and free correctly)
     **/
    inline void* allocateMemory( size_t i_size,
                                 size_t i_alignment,
                                 int    i_memkind = 0 )
    {
      void* l_ptrBuffer;
      bool error = false;

#ifdef USE_MEMKIND
      if( i_memkind == 0 ) {
#endif
        if (i_alignment % (sizeof(void*)) != 0) {
          l_ptrBuffer = malloc( i_size );
          error = (l_ptrBuffer == NULL);
        } else {
          error = (posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
        }
#ifdef USE_MEMKIND
      } else {
        if (i_alignment % (sizeof(void*)) != 0) {
          l_ptrBuffer = hbw_malloc( i_size );
          error = (l_ptrBuffer == NULL);
        } else {
          error = (hbw_posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
        } 
      }
#endif
      
      if (!error) {
        m_dataMemoryAddresses.push_back( Address(i_memkind, l_ptrBuffer) );
      } else {
        logError() << "The malloc failed (bytes: " << i_size << ", alignment: " << i_alignment << ", memkind: " << i_memkind << ").";
      }

      return l_ptrBuffer;
    }

    /**
     * Frees all memory, which was allocated by functions of the MemoryAllocator.
     **/
    void freeMemory();
};

#endif
