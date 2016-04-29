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
 * General purpose node for a data tree.
 **/

#ifndef INITIALIZER_TREE_NODE_HPP_
#define INITIALIZER_TREE_NODE_HPP_

#include <cstddef>

namespace seissol {
  namespace initializers {
    template<typename T, unsigned NUM_CHILDREN>
    class Children {
    private:
      T m_storage[NUM_CHILDREN];
    public:
      static const unsigned size = NUM_CHILDREN;
      
      T& operator[](unsigned index) { return m_storage[index]; }
      T const& operator[](unsigned index) const { return m_storage[index]; }
    };

    template<typename T>
    class Children<T, 0> {
    private:
      T* m_storage;
    public:
      unsigned size;

      T& operator[](unsigned index) { return m_storage[index]; }
      T const& operator[](unsigned index) const { return m_storage[index]; }
      
      Children() : m_storage(NULL), size(0) {}
      ~Children() { delete[] m_storage; }
      
      void allocate(unsigned numChildren) {
        delete[] m_storage;
        m_storage = new T[numChildren];
        size = numChildren;
      }
    };

    template<typename T, unsigned NUM_CHILDREN, unsigned NUM_VARIABLES, unsigned NUM_BUCKETS>
    class Node {
    protected:
      Children<T, NUM_CHILDREN> m_children;

    public:
      void addBytes(size_t varBytes[NUM_VARIABLES], size_t bucketBytes[NUM_BUCKETS]) {
        for (unsigned i = 0; i < m_children.size; ++i) {
          m_children[i].addBytes(varBytes, bucketBytes);
        }
      }
      
      void setMemoryRegions(void* vars[NUM_VARIABLES], void* buckets[NUM_BUCKETS]) {
        size_t varBytes[NUM_VARIABLES] = {0};
        size_t bucketBytes[NUM_BUCKETS] = {0};
        void* varPtrs[NUM_VARIABLES];
        void* bucketPtrs[NUM_BUCKETS];
        for (unsigned i = 0; i < m_children.size; ++i) {
          for (unsigned j = 0; j < NUM_VARIABLES; ++j) {
            varPtrs[j] = static_cast<char*>(vars[j]) + varBytes[j];
          }
          for (unsigned j = 0; j < NUM_BUCKETS; ++j) {
            bucketPtrs[j] = static_cast<char*>(buckets[j]) + bucketBytes[j];
          }
          m_children[i].setMemoryRegions(varPtrs, bucketPtrs);
          m_children[i].addBytes(varBytes, bucketBytes);
        }
      }
      
      void touch() {
        for (unsigned i = 0; i < m_children.size; ++i) {
          m_children[i].touch();
        }
      }

      unsigned getNumberOfCells() {
        unsigned numberOfCells = 0;
        for (unsigned i = 0; i < m_children.size; ++i) {
          numberOfCells += m_children[i].getNumberOfCells();
        }
        return numberOfCells;
      }
      
      inline T& child(unsigned index) {
        return m_children[index];
      }
      
      inline T const& child(unsigned index) const {
        return m_children[index];
      }
      
      inline unsigned numChildren() const {
        return m_children.size;
      }
    };
  }
}

#endif
