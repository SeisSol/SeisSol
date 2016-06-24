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

#ifndef INTERNALSTATE_H_
#define INTERNALSTATE_H_

#include <cstddef>
#include <Initializer/typedefs.hpp>

namespace seissol {
  namespace initializers {
    class InternalState;
  }
}

/**
 * Setup of SeisSols internal state.
 **/
class seissol::initializers::InternalState {
  //private:

  public:
    /**
     * Derives the layout of either the ghost or copy layer or interior.
     *
     * @param i_numberOfCluster number of time stepping clusters.
     * @param i_numberOfRegions number of communication regions per cluster.
     * @param i_numberOfRegionCells number of cells in the communication regions.
     * @param i_cellLocalInformation cell local information; (points to information of the first cell in the first cluster and communication region).
     * @param o_numberOfBuffers number of cells with time buffers per cluster and region.
     * @param o_numberOfDerivatives number of cells with derivatives per cluster and region.
     **/
    static void deriveLayerLayout(        unsigned int                  i_numberOfClusters,
                                          unsigned int                 *i_numberOfRegions,
                                          unsigned int                **i_numberOfRegionCells,
                                    const struct CellLocalInformation  *i_cellLocalInformation,
                                          unsigned int                **o_numberOfBuffers,
                                          unsigned int                **o_numberOfDerivatives );

    /**
     * Derives the size of the interior time buffers/derivatives.
     *
     * @param i_numberOfCluster number of time stepping clusters.
     * @param i_numberOfInteriorCells number of interior cells per cluster.
     * @param i_cellLocalInformation cell local information; (points to first interior cell).
     * @param o_numberOfBuffers set to: number of cells with time buffers in every cluster.
     * @param o_numberOfDerivatives set to: number of cells with derivatives in every cluster.
     **/
    static void deriveInteriorLayout(       unsigned int                 i_numberOfClusters,
                                            unsigned int                *i_numberOfInteriorCells,
                                      const struct CellLocalInformation *i_cellLocalInformation,
                                            unsigned int                *o_numberOfBuffers,
                                            unsigned int                *o_numberOfDerivatives );

    /**
     * Sets up the pointers to time buffers/derivatives of the ghost or copy layer or interior.
     *
     * @param i_numberOfRegions number of communication regions.
     * @param i_numberOfRegionCells number of cells in the regions.
     * @param i_cellLocalInformation cell local information (points to first cell in the layer).
     * @param i_numberOfBuffers number of cells with time buffers per region.
     * @param i_numberOfDerivatives number of cells with derivatives per region.
     * @param i_layerMemory layer in memory.
     * @param o_buffers pointers will be set to time buffers (first pointer belongs to first region cell); set to NULL if no buffer exists.
     * @param o_derivatives pointers will be set to time derivatives (first pointer belongs to first region cell); set to NULL if no derivative exists.
     **/
    static void setUpLayerPointers(       unsigned int                 i_numberOfRegions,
                                    const unsigned int                *i_numberOfRegionCells,
                                    const struct CellLocalInformation *i_cellLocalInformation,
                                    const unsigned int                *i_numberOfBuffers,
                                    const unsigned int                *i_numberOfDerivatives,
                                          real                        *i_layerMemory,
                                          real                       **o_buffers,
                                          real                       **o_derivatives );

    /**
     * Sets up the pointers to time buffers/derivatives in the interior of the computational domain.
     *
     * @param i_numberOfInteriorCells number of cells in the interior.
     * @param i_cellLocalInformation cell local information (points to first cell local information in the interior).
     * @param i_numberOfBuffers number of cells with buffers in the interior.
     * @param i_numberOfDerivatives number of cells with derivatives in the interioir.
     * @param i_interiorMemory chunk of memory for the internal state of the interior.
     * @param o_buffers will be set time buffers (first pointer belongs to first cell in the interior); set to NULL if no buffer for a cell exists.
     * @param o_derivatives will be set to time derivatives (first pointer belongs to first cell interior); set to NULL if no derivatives exist for a cell.
     **/
    static void setUpInteriorPointers(       unsigned int                  i_numberOfInteriorCells,
                                       const struct CellLocalInformation  *i_cellLocalInformation,
                                             unsigned int                  i_numberOfBuffers,
                                             unsigned int                  i_numberOfDerivatives,
                                             real                         *i_interiorMemory,
                                             real                        **o_buffers,
                                             real                        **o_derivatives );
};

#endif
