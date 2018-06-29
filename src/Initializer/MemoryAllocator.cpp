/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
#include "MemoryAllocator.h"

#include <utils/logger.h>

void* seissol::memory::allocate(size_t i_size, size_t i_alignment, enum Memkind i_memkind)
{
    void* l_ptrBuffer;
    bool error = false;

    /* handle zero allocation */
    if ( i_size == 0 ) {
      //logWarning() << "allocation of size 0 requested, returning NULL; (alignment: " << i_alignment << ", memkind: " << i_memkind << ").";
      l_ptrBuffer = NULL;
      return l_ptrBuffer;
    }

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
    
    if (error) {
      logError() << "The malloc failed (bytes: " << i_size << ", alignment: " << i_alignment << ", memkind: " << i_memkind << ").";
    }

    return l_ptrBuffer;
}

void seissol::memory::free(void* i_pointer, enum Memkind i_memkind)
{
#ifdef USE_MEMKIND
    if (i_memkind == Standard) {
#endif
      ::free( i_pointer );
#ifdef USE_MEMKIND
    } else {
      hbw_free( i_pointer );
    }
#endif
}

void seissol::memory::printMemoryAlignment( std::vector< std::vector<unsigned long long> > i_memoryAlignment ) {
  logDebug() << "printing memory alignment per struct";
  for( unsigned long long l_i = 0; l_i < i_memoryAlignment.size(); l_i++ ) {
    logDebug() << i_memoryAlignment[l_i][0] << ", " << i_memoryAlignment[l_i][1];
  }
}

seissol::memory::ManagedAllocator::~ManagedAllocator()
{
  for (AddressVector::const_iterator it = m_dataMemoryAddresses.begin(); it != m_dataMemoryAddresses.end(); ++it) {
    seissol::memory::free(it->second, it->first);
  }

  // reset memory vectors
  m_dataMemoryAddresses.clear();
}

void* seissol::memory::ManagedAllocator::allocateMemory( size_t i_size, size_t i_alignment, enum Memkind i_memkind )
{
  void* l_ptrBuffer = seissol::memory::allocate(i_size, i_alignment, i_memkind);
  m_dataMemoryAddresses.push_back( Address(i_memkind, l_ptrBuffer) );
  return l_ptrBuffer;
}
