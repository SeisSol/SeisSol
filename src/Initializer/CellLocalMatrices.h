// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_
#define SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Boundary.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/LTS.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/Typedefs.h"

namespace seissol {
  namespace initializer {
      class EasiBoundary;
      /**
      * Computes the star matrices A*, B*, and C*, and solves the Riemann problems at the interfaces.
      **/
     void initializeCellLocalMatrices( seissol::geometry::MeshReader const&      i_meshReader,
                                       LTSTree*               io_ltsTree,
                                       LTS*                   i_lts,
                                       Lut*                   i_ltsLut,
                                       TimeStepping const&    timeStepping,
                                       const parameters::ModelParameters& modelParameters );
                                       
     void initializeBoundaryMappings(seissol::geometry::MeshReader const& i_meshReader,
                                     const EasiBoundary* easiBoundary,
                                     LTSTree* io_ltsTree,
                                     LTS* i_lts,
                                     Lut* i_ltsLut);
 
     void initializeDynamicRuptureMatrices( seissol::geometry::MeshReader const&      i_meshReader,                                                    
                                            LTSTree*               io_ltsTree,
                                            LTS*                   i_lts,
                                            Lut*                   i_ltsLut,
                                            LTSTree*               dynRupTree,
                                            DynamicRupture*        dynRup,
                                            unsigned*              ltsFaceToMeshFace,
                                            GlobalData const&      global,
                                            double etaHack );

      void copyCellMatricesToDevice(LTSTree*          ltsTree,
                                    LTS*              lts,
                                    LTSTree*          dynRupTree,
                                    DynamicRupture*   dynRup,
                                    LTSTree*          boundaryTree,
                                    Boundary*         boundary);
  }
}


#endif // SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_

