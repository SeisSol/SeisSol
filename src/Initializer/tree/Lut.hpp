/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 * Tree for managing lts data.
 **/
 
#ifndef INITIALIZER_TREE_LUT_HPP_
#define INITIALIZER_TREE_LUT_HPP_

#include "foreach.hpp"
#include "LTSTree.hpp"

#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace initializers {
    template<typename Spec>
    class Lut {
    private:
      unsigned*                         m_ltsToMesh;
      unsigned*                         m_meshToCell[1 << Spec::NUM_LAYERS];
      LTSTree<Spec>*                    m_ltsTree;
      seissol::memory::ManagedAllocator m_allocator;
      
      typedef BoolToType<true> ContinueRecursion;
      typedef BoolToType<false>  TerminateRecursion;
      template<unsigned INDEX, unsigned LAYER>
      static unsigned calcLutNumber(ContinueRecursion) {
        return (Spec::Available[INDEX][LAYER] << LAYER) | calcLutNumber<INDEX, LAYER+1>(BoolToType<LAYER+1 < Spec::NUM_LAYERS>());
      };
      
      template<unsigned INDEX, unsigned LAYER>
      static unsigned calcLutNumber(TerminateRecursion) {
        return 0;
      };

      template<unsigned INDEX>
      static unsigned lutNumber() {
        return calcLutNumber<INDEX, 0>(BoolToType<0 < Spec::NUM_LAYERS>());
      }
      
      /// \todo Could be generalized to work with any tree
      template<unsigned INDEX>
      struct createLutForIndex {
        void operator()(  Lut<Spec>*            self,
                          unsigned              numberOfMeshIds  ) const {
          unsigned const lut = lutNumber<INDEX>();
          if (self->m_meshToCell[lut] == NULL) {
            unsigned startLtsId = 0;
            unsigned offset = 0;
            self->m_meshToCell[lut] = static_cast<unsigned*>(self->m_allocator.allocateMemory(numberOfMeshIds * sizeof(unsigned)));
            std::fill(self->m_meshToCell[lut], self->m_meshToCell[lut] + numberOfMeshIds, std::numeric_limits<unsigned>::max());
            
            // Every time cluster
            for (unsigned tc = 0; tc < self->m_ltsTree->numChildren(); ++tc) {
              TimeCluster<Spec> const& cluster = self->m_ltsTree->child(tc);
              // Every layer
              for (unsigned l = 0; l < cluster.numChildren(); ++l) {
                Layer<Spec> const& layer = cluster.child(l);
                // If data is saved on layer
                if (Spec::Available[INDEX][l]) {
                  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {      
                    unsigned meshId = self->m_ltsToMesh[startLtsId + cell];
                    if (meshId != std::numeric_limits<unsigned>::max()) {
                      self->m_meshToCell[lut][meshId] = offset + cell;
                    }
                  }
                  offset += layer.getNumberOfCells();
                }
                startLtsId += layer.getNumberOfCells();
              }
            }
          }
        }
      };

    public:
      Lut() : m_ltsToMesh(NULL), m_ltsTree(NULL) {
        for (unsigned i = 0; i < (1 << Spec::NUM_LAYERS); ++i) {
          m_meshToCell[i] = NULL;
        }
      }
      
      /// \todo Could be generalized to work with any tree
      void createLuts(  LTSTree<Spec>*  ltsTree,
                        unsigned*       ltsToMesh,
                        unsigned        numberOfCells,
                        unsigned        numberOfMeshIds ) {
        assert(numberOfCells == ltstree->getNumberOfCells());

        m_ltsToMesh = static_cast<unsigned*>(m_allocator.allocateMemory(numberOfCells * sizeof(unsigned)));
        std::copy(ltsToMesh, ltsToMesh + numberOfCells, m_ltsToMesh);
        m_ltsTree = ltsTree;
        ForEachClassMethod<createLutForIndex, Spec::NUM_VARIABLES>(this, numberOfMeshIds);
      }
      
      inline unsigned meshId(unsigned ltsId) const {
        return m_ltsToMesh[ltsId];
      }
      
      template<unsigned INDEX>
      typename get_type<typename Spec::Types, INDEX>::type& lookup(unsigned meshId) {
        unsigned offset = m_meshToCell[lutNumber<INDEX>()][meshId];
        return m_ltsTree->template var<INDEX>()[offset];
      }
    };
  }
}

#endif
