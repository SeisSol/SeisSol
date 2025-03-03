// SPDX-FileCopyrightText: 2014-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_INITIALIZER_INTERNALSTATE_H_
#define SEISSOL_SRC_INITIALIZER_INTERNALSTATE_H_

#include <cstddef>
#include "Initializer/Typedefs.h"

namespace seissol {
  namespace initializer {
    class InternalState;
  }
}

/**
 * Setup of SeisSols internal state.
 **/
class seissol::initializer::InternalState {
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


#endif // SEISSOL_SRC_INITIALIZER_INTERNALSTATE_H_

