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
 
#ifndef INITIALIZER_TREE_LTSTREE_HPP_
#define INITIALIZER_TREE_LTSTREE_HPP_

#include "typelist.hpp"
#include "foreach.hpp"
#include "Node.hpp"

#include <Initializer/MemoryAllocator.h>

#include <cstring>

namespace seissol {
  namespace initializers {
    template<typename Spec>
    class Layer {
    private:
      unsigned m_layer;
      unsigned m_numberOfCells;
      void* m_vars[Spec::NUM_VARIABLES];
      void* m_buckets[Spec::NUM_BUCKETS];
      size_t m_bucketSizes[Spec::NUM_BUCKETS];

    public:
      Layer() : m_layer(-1), m_numberOfCells(0) {    
        for (unsigned bucket = 0; bucket < Spec::NUM_BUCKETS; ++bucket) {
          m_bucketSizes[bucket] = 0;
        }
      }

      template<unsigned INDEX>
      typename get_type<typename Spec::Types, INDEX>::type* var() {
        return static_cast<typename get_type<typename Spec::Types, INDEX>::type*>(m_vars[INDEX]);
      }

      inline unsigned getLayerType() const {
        return m_layer;
      }      

      inline void setLayerType(unsigned layerType) {
        m_layer = layerType;
      }
      
      inline unsigned getNumberOfCells() const {
        return m_numberOfCells;
      }
      
      inline void setNumberOfCells(unsigned numberOfCells) {
        m_numberOfCells = numberOfCells;
      }
      
      inline void* bucket(unsigned index) {
        return m_buckets[index];
      }
      
      inline void setBucketSize(unsigned index, size_t size) {
        m_bucketSizes[index] = size;
      }
      
      void addBytes(size_t varBytes[Spec::NUM_VARIABLES], size_t bucketBytes[Spec::NUM_BUCKETS]) {
        ForEachClassMethod<addBytesForIndex, Spec::NUM_VARIABLES>(this, varBytes);
        for (unsigned bucket = 0; bucket < Spec::NUM_BUCKETS; ++bucket) {
          bucketBytes[bucket] += m_bucketSizes[bucket];
        }
      }

      template<unsigned INDEX>
      struct addBytesForIndex {
        void operator()(Layer<Spec>* self, size_t varBytes[Spec::NUM_VARIABLES]) const {
          if (Spec::Available[INDEX][self->m_layer]) {
            varBytes[INDEX] += sizeof(typename get_type<typename Spec::Types, INDEX>::type) * self->m_numberOfCells;
          }
        }
      };
      
      void setMemoryRegions(void* vars[Spec::NUM_VARIABLES], void* buckets[Spec::NUM_BUCKETS]) {
        ForEachClassMethod<setMemoryRegionForIndex, Spec::NUM_VARIABLES>(this, vars);
        for (unsigned b = 0; b < Spec::NUM_BUCKETS; ++b) {
          m_buckets[b] = (m_bucketSizes[b] != 0) ? buckets[b] : NULL;
        }
      }
      
      template<unsigned INDEX>
      struct setMemoryRegionForIndex {
        void operator()(Layer<Spec>* self, void* vars[Spec::NUM_VARIABLES]) const {
          self->m_vars[INDEX] = (Spec::Available[INDEX][self->m_layer]) ? vars[INDEX] : NULL;
        }
      };
      
      void touch() {
        ForEachClassMethod<touchForIndex, Spec::NUM_VARIABLES>(this);
      }
      
      template<unsigned INDEX>
      struct touchForIndex {
        void operator()(Layer<Spec>* self) const {
          if (Spec::Available[INDEX][self->m_layer]) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
            for (unsigned cell = 0; cell < self->m_numberOfCells; ++cell) {
              memset(self->var<INDEX>() + cell, 0, sizeof(typename get_type<typename Spec::Types, INDEX>::type));
            }
          }
        }
      };
    };

    template<typename Spec>
    class TimeCluster : public Node<Layer<Spec>, Spec::NUM_LAYERS, Spec::NUM_VARIABLES, Spec::NUM_BUCKETS> {
    public: 
      TimeCluster() {
        this->m_children[0].setLayerType(Spec::Ghost);
        this->m_children[1].setLayerType(Spec::Copy);
        this->m_children[2].setLayerType(Spec::Interior);
      }
    };

    template<typename Spec>
    class LTSTree : public Node<TimeCluster<Spec>, 0, Spec::NUM_VARIABLES, Spec::NUM_BUCKETS> {
    private:
      void* m_vars[Spec::NUM_VARIABLES];
      void* m_buckets[Spec::NUM_BUCKETS];
      unsigned m_numberOfClusters;
      seissol::memory::ManagedAllocator m_allocator;

    public:
      LTSTree() {
        for (unsigned i = 0; i < Spec::NUM_VARIABLES; ++i) {
          m_vars[i] = NULL;
        }
        for (unsigned i = 0; i < Spec::NUM_BUCKETS; ++i) {
          m_buckets[i] = NULL;
        }
      }
      
      /// \todo remove and replace by a proper concept
      template<unsigned INDEX>
      typename get_type<typename Spec::Types, INDEX>::type* var() {
        return static_cast<typename get_type<typename Spec::Types, INDEX>::type*>(m_vars[INDEX]);
      }
      
      inline void setNumberOfTimeClusters(unsigned numberOfTimeClusters) {
        this->m_children.allocate(numberOfTimeClusters);
      }
      
      void printRequiredMemory() {
        size_t varBytes[Spec::NUM_VARIABLES] = {0};
        size_t bucketBytes[Spec::NUM_BUCKETS] = {0};
        addBytes(varBytes, bucketBytes);
        
        size_t total = 0;
        for (unsigned i = 0; i < Spec::NUM_VARIABLES; ++i) {
          std::cout << "Var " << i << ": " << varBytes[i] << std::endl;
          total += varBytes[i];
        }
        
        for (unsigned i = 0; i < Spec::NUM_BUCKETS; ++i) {
          std::cout << "Bucket " << i << ": " << bucketBytes[i] << std::endl;
          total += bucketBytes[i];
        }
        
        std::cout << total / 1073741824.0 << " GiB" << std::endl;
      }
      
      void allocateMemory() {
        size_t varBytes[Spec::NUM_VARIABLES] = {0};
        size_t bucketBytes[Spec::NUM_BUCKETS] = {0};
        addBytes(varBytes, bucketBytes);
        for (unsigned i = 0; i < Spec::NUM_VARIABLES; ++i) {
          m_vars[i] = m_allocator.allocateMemory(varBytes[i], Spec::VarAlignment[i], Spec::VarMemkind[i]);
        }
        for (unsigned i = 0; i < Spec::NUM_BUCKETS; ++i) {
          m_buckets[i] = m_allocator.allocateMemory(bucketBytes[i], Spec::BucketAlignment[i], Spec::BucketMemkind[i]);
        }
            
        setMemoryRegions(m_vars, m_buckets);
      }
    };
  }
}

#endif
