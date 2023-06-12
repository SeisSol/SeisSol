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

#include "typedefs.hpp"
#include "NRF.h"
#include <cstdarg>

#include <Initializer/tree/Lut.hpp>
#include <Kernels/PointSourceCluster.h>
#include <Solver/time_stepping/TimeManager.h>
#include <Geometry/MeshReader.h>
#include <inttypes.h>
#include <memory>
#include <array>
#include <vector>

namespace seissol {
namespace sourceterm {
enum class SourceType : int { None = 0, NrfSource = 42, FsrmSource = 50 };

void computeMInvJInvPhisAtSources(
    Eigen::Vector3d const& centre,
    AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>& mInvJInvPhisAtSources,
    unsigned meshId,
    seissol::geometry::MeshReader const& mesh);
void transformNRFSourceToInternalSource(Eigen::Vector3d const& centre,
                                        unsigned meshId,
                                        seissol::geometry::MeshReader const& mesh,
                                        Subfault const& subfault,
                                        Offsets const& offsets,
                                        Offsets const& nextOffsets,
                                        std::array<std::vector<double>, 3> const& sliprates,
                                        seissol::model::Material* material,
                                        PointSources& pointSources,
                                        unsigned index,
                                        AllocatorT const& alloc);
class Manager;
} // namespace sourceterm
} // namespace seissol

class seissol::sourceterm::Manager {
  public:
  Manager() = default;
  ~Manager() = default;

  void loadSources(SourceType sourceType,
                   char const* fileName,
                   seissol::geometry::MeshReader const& mesh,
                   seissol::initializers::LTSTree* ltsTree,
                   seissol::initializers::LTS* lts,
                   seissol::initializers::Lut* ltsLut,
                   time_stepping::TimeManager& timeManager);

  private:
  auto mapPointSourcesToClusters(const unsigned* meshIds,
                                 unsigned numberOfSources,
                                 seissol::initializers::LTSTree* ltsTree,
                                 seissol::initializers::LTS* lts,
                                 seissol::initializers::Lut* ltsLut,
                                 AllocatorT const& alloc)
      -> std::unordered_map<LayerType, std::vector<ClusterMapping>>;

  auto makePointSourceCluster(ClusterMapping mapping, PointSources sources)
      -> std::unique_ptr<kernels::PointSourceCluster>;

  auto loadSourcesFromFSRM(char const* fileName,
                           seissol::geometry::MeshReader const& mesh,
                           seissol::initializers::LTSTree* ltsTree,
                           seissol::initializers::LTS* lts,
                           seissol::initializers::Lut* ltsLut,
                           AllocatorT const& alloc)
      -> std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>;

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
  auto loadSourcesFromNRF(char const* fileName,
                          seissol::geometry::MeshReader const& mesh,
                          seissol::initializers::LTSTree* ltsTree,
                          seissol::initializers::LTS* lts,
                          seissol::initializers::Lut* ltsLut,
                          AllocatorT const& alloc)
      -> std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>;
#endif
};

#endif
