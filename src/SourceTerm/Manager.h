/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include "NRF.h"
#include "typedefs.hpp"
#include <cstdarg>

#include "Geometry/MeshReader.h"
#include "Initializer/tree/Lut.hpp"
#include "Kernels/PointSourceCluster.h"
#include "Solver/time_stepping/TimeManager.h"
#include <array>
#include <inttypes.h>
#include <memory>
#include <vector>

namespace seissol::sourceterm {

void computeMInvJInvPhisAtSources(
    const Eigen::Vector3d& centre,
    AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>& mInvJInvPhisAtSources,
    unsigned meshId,
    const seissol::geometry::MeshReader& mesh);
void transformNRFSourceToInternalSource(const Eigen::Vector3d& centre,
                                        unsigned meshId,
                                        const seissol::geometry::MeshReader& mesh,
                                        const Subfault& subfault,
                                        const Offsets& offsets,
                                        const Offsets& nextOffsets,
                                        const std::array<std::vector<double>, 3>& sliprates,
                                        seissol::model::Material* material,
                                        PointSources& pointSources,
                                        unsigned index,
                                        const AllocatorT& alloc);
class Manager;
} // namespace seissol::sourceterm

class seissol::sourceterm::Manager {
  public:
  Manager() = default;
  ~Manager() = default;

  void loadSources(seissol::initializer::parameters::PointSourceType sourceType,
                   const char* fileName,
                   const seissol::geometry::MeshReader& mesh,
                   seissol::initializer::LTSTree* ltsTree,
                   seissol::initializer::LTS* lts,
                   seissol::initializer::Lut* ltsLut,
                   time_stepping::TimeManager& timeManager);

  private:
  auto mapPointSourcesToClusters(const unsigned* meshIds,
                                 unsigned numberOfSources,
                                 seissol::initializer::LTSTree* ltsTree,
                                 seissol::initializer::LTS* lts,
                                 seissol::initializer::Lut* ltsLut,
                                 const AllocatorT& alloc)
      -> std::unordered_map<LayerType, std::vector<ClusterMapping>>;

  auto makePointSourceCluster(ClusterMapping mapping,
                              PointSources sources) -> std::unique_ptr<kernels::PointSourceCluster>;

  auto loadSourcesFromFSRM(const char* fileName,
                           const seissol::geometry::MeshReader& mesh,
                           seissol::initializer::LTSTree* ltsTree,
                           seissol::initializer::LTS* lts,
                           seissol::initializer::Lut* ltsLut,
                           const AllocatorT& alloc)
      -> std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>;

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
  auto loadSourcesFromNRF(const char* fileName,
                          const seissol::geometry::MeshReader& mesh,
                          seissol::initializer::LTSTree* ltsTree,
                          seissol::initializer::LTS* lts,
                          seissol::initializer::Lut* ltsLut,
                          const AllocatorT& alloc)
      -> std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>;
#endif
};

#endif
