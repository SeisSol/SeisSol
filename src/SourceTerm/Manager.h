/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 **/

#ifndef SOURCETERM_MANAGER_H_
#define SOURCETERM_MANAGER_H_

#include "typedefs.hpp"
#include "NRF.h"

#include <Initializer/tree/Lut.hpp>
#include <Solver/time_stepping/TimeManager.h>
#include <Geometry/MeshReader.h>
#include <inttypes.h>

namespace seissol {
  namespace sourceterm {
    void transformNRFSourceToInternalSource( glm::dvec3 const&         centre,
                                             unsigned                 element,
                                             Subfault const&          subfault,
                                             Offsets const&           offsets,
                                             Offsets const&           nextOffsets,
                                             double *const            sliprates[3],
                                             seissol::model::Material material,
                                             PointSources&            pointSources,
                                             unsigned                 index );
    class Manager;
  }
}

class seissol::sourceterm::Manager {
private:
  ClusterMapping* cmps;
  PointSources* sources;

  void freeSources();

public:
  Manager() : cmps(NULL), sources(NULL) {}
  ~Manager() { freeSources(); }
  
  void mapPointSourcesToClusters( unsigned const*                 meshIds,
                                  unsigned                        numberOfSources,
                                  seissol::initializers::LTSTree* ltsTree,
                                  seissol::initializers::LTS*     lts,
                                  seissol::initializers::Lut*     ltsLut );

  void loadSourcesFromFSRM( double const*                   momentTensor,
                            double const*                   velocityComponent,
                            int                             numberOfSources,
                            double const*                   centres,
                            double const*                   strikes,
                            double const*                   dips,
                            double const*                   rakes,
                            double const*                   onsets,
                            double const*                   areas,
                            double                          timestep,
                            int                             numberOfSamples,
                            double const*                   timeHistories,
                            MeshReader const&               mesh,
                            seissol::initializers::LTSTree* ltsTree,
                            seissol::initializers::LTS*     lts,
                            seissol::initializers::Lut*     ltsLut,
                            time_stepping::TimeManager&     timeManager);

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
  void loadSourcesFromNRF(  char const*                     fileName,
                            MeshReader const&               mesh,
                            seissol::initializers::LTSTree* ltsTree,
                            seissol::initializers::LTS*     lts,
                            seissol::initializers::Lut*     ltsLut,
                            time_stepping::TimeManager&     timeManager );
#endif
};

#endif
