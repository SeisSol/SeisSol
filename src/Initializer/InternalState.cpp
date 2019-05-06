/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
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
 * Setup of SeisSol's internal state.
 **/

#include "InternalState.h"
#include <limits>
#include <cstddef>
#include <cassert>
#include <yateto.h>

void seissol::initializers::InternalState::deriveLayerLayout(       unsigned int                  i_numberOfClusters,
                                                                    unsigned int                 *i_numberOfRegions,
                                                                    unsigned int                **i_numberOfRegionCells,
                                                              const struct CellLocalInformation  *i_cellLocalInformation,
                                                                    unsigned int                **o_numberOfBuffers,
                                                                    unsigned int                **o_numberOfDerivatives ) {
  unsigned int l_firstRegionCell = 0;

  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    // iterate over all regions
    for( unsigned int l_region = 0; l_region < i_numberOfRegions[l_cluster]; l_region++ ) {
      // reset values
      o_numberOfBuffers[l_cluster][l_region] = 0; o_numberOfDerivatives[l_cluster][l_region] = 0;

      unsigned int l_firstNonRegionCell = l_firstRegionCell +
                                          i_numberOfRegionCells[l_cluster][l_region];

      for( unsigned int l_cell = l_firstRegionCell; l_cell < l_firstNonRegionCell; l_cell++ ) {
        if( (i_cellLocalInformation[l_cell].ltsSetup >> 8 ) % 2 == 1 ) o_numberOfBuffers[l_cluster][l_region]++;
        if( (i_cellLocalInformation[l_cell].ltsSetup >> 9 ) % 2 == 1 ) o_numberOfDerivatives[l_cluster][l_region]++;
      }

      l_firstRegionCell = l_firstNonRegionCell;
    }
  }
}

void seissol::initializers::InternalState::deriveInteriorLayout(       unsigned int                  i_numberOfClusters,
                                                                       unsigned int                 *i_numberOfInteriorCells,
                                                                 const struct CellLocalInformation  *i_cellLocalInformation,
                                                                       unsigned int                 *o_numberOfBuffers,
                                                                       unsigned int                 *o_numberOfDerivatives ) {
  // set up dummy regions of size 1 and dummy cells
  unsigned int  *l_numberOfRegions     = new unsigned int [i_numberOfClusters];
  unsigned int **l_numberOfRegionCells = new unsigned int*[i_numberOfClusters];
  unsigned int **l_numberOfBuffers     = new unsigned int*[i_numberOfClusters];
  unsigned int **l_numberOfDerivatives = new unsigned int*[i_numberOfClusters];

  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    l_numberOfRegions[l_cluster] = 1;

    l_numberOfRegionCells[l_cluster]    = new unsigned int[1];
    l_numberOfBuffers[l_cluster]        = &o_numberOfBuffers[l_cluster];
    l_numberOfDerivatives[l_cluster]    = &o_numberOfDerivatives[l_cluster];

    l_numberOfRegionCells[l_cluster][0] = i_numberOfInteriorCells[l_cluster];
  }

  // interior is a special case with a single layer
  deriveLayerLayout(  i_numberOfClusters,
                      l_numberOfRegions,
                      l_numberOfRegionCells,
                      i_cellLocalInformation,
                      l_numberOfBuffers,
                      l_numberOfDerivatives );

  // free memory
  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    delete[] l_numberOfRegionCells[l_cluster];
  }
  delete[] l_numberOfBuffers;
  delete[] l_numberOfDerivatives;
  delete[] l_numberOfRegions;
  delete[] l_numberOfRegionCells;
}

void seissol::initializers::InternalState::setUpLayerPointers(       unsigned int                 i_numberOfRegions,
                                                               const unsigned int                *i_numberOfRegionCells,
                                                               const struct CellLocalInformation *i_cellLocalInformation,
                                                               const unsigned int                *i_numberOfBuffers,
                                                               const unsigned int                *i_numberOfDerivatives,
                                                                     real                        *i_layerMemory,
                                                                     real                       **o_buffers,
                                                                     real                       **o_derivatives ) {
  // first cell of the current region
  unsigned int l_firstRegionCell = 0;

  // offset in the layer to the current region
  unsigned int l_offset = 0;

  // iterate over all regions
  for( unsigned int l_region = 0; l_region < i_numberOfRegions; l_region++ ) {
    // first cell belonging to the next region, if present
    unsigned int l_firstNonRegionCell = l_firstRegionCell +
                                        i_numberOfRegionCells[l_region];

    // counters to dertermine the current position in memory of the buffers/derivatives
    unsigned int l_bufferCounter = 0;
    unsigned int l_derivativeCounter = 0;

    // iterate over this particular region
    for( unsigned int l_cell = l_firstRegionCell; l_cell < l_firstNonRegionCell; l_cell++ ) {
      // set pointers and increase conunters
      if( (i_cellLocalInformation[l_cell].ltsSetup >> 8 ) % 2 ) {
        o_buffers[l_cell] = i_layerMemory + l_offset
                                          + l_bufferCounter * tensor::I::size();
        l_bufferCounter++;
      }
      else o_buffers[l_cell] = NULL;

      if( (i_cellLocalInformation[l_cell].ltsSetup >> 9 ) % 2 ) {
        o_derivatives[l_cell] = i_layerMemory + l_offset 
                                              + i_numberOfBuffers[l_region] * tensor::I::size()
                                              + l_derivativeCounter * yateto::computeFamilySize<tensor::dQ>();
        l_derivativeCounter++;
      }
      else o_derivatives[l_cell] = NULL;
    }

    // check that we have all buffers and derivatives
    assert( l_bufferCounter     == i_numberOfBuffers[l_region] );
    assert( l_derivativeCounter == i_numberOfDerivatives[l_region] );

    // update offsets
    l_firstRegionCell = l_firstNonRegionCell;
    l_offset += i_numberOfBuffers[l_region]     * tensor::I::size() +
                i_numberOfDerivatives[l_region] * yateto::computeFamilySize<tensor::dQ>();
  }
}

void seissol::initializers::InternalState::setUpInteriorPointers(       unsigned int                  i_numberOfInteriorCells,
                                                                  const struct CellLocalInformation  *i_cellLocalInformation,
                                                                        unsigned int                  i_numberOfBuffers,
                                                                        unsigned int                  i_numberOfDerivatives,
                                                                        real                         *i_interiorMemory,
                                                                        real                        **o_buffers,
                                                                        real                        **o_derivatives ) {
  // interior is a special layered case with a single region
  setUpLayerPointers(  1,
                      &i_numberOfInteriorCells,
                       i_cellLocalInformation,
                      &i_numberOfBuffers,
                      &i_numberOfDerivatives,
                       i_interiorMemory,
                       o_buffers,
                       o_derivatives );
}
