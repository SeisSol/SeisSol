/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
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

#include "Parallel/MPI.h"

#include "FSRMReader.h"
#include "Initializer/PointMapper.h"
#include "Kernels/PointSourceClusterOnHost.h"
#include "Manager.h"
#include "Numerical/Transformation.h"
#include "PointSource.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <Common/Constants.h>
#include <Geometry/MeshReader.h>
#include <Geometry/MeshTools.h>
#include <Initializer/LTS.h>
#include <Initializer/MemoryAllocator.h>
#include <Initializer/Parameters/SourceParameters.h>
#include <Initializer/Tree/LTSTree.h>
#include <Initializer/Tree/Layer.h>
#include <Initializer/Tree/Lut.h>
#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Model/CommonDatastructures.h>
#include <Numerical/BasisFunction.h>
#include <Solver/time_stepping/TimeManager.h>
#include <SourceTerm/NRF.h>
#include <SourceTerm/Typedefs.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

#include <Model/Datastructures.h> // IWYU pragma: keep

#ifdef USE_NETCDF
#include "NRFReader.h"
#include <mpi.h>
#endif

#ifdef ACL_DEVICE
#include "Kernels/PointSourceClusterOnDevice.h"
#include <Parallel/Helper.h>
#endif

namespace {

using namespace seissol::sourceterm;

/**
 * Computes mInvJInvPhisAtSources[i] = |J|^-1 * M_ii^-1 * phi_i(xi, eta, zeta),
 * where xi, eta, zeta is the point in the reference tetrahedron corresponding to x, y, z.
 */
void computeMInvJInvPhisAtSources(
    const Eigen::Vector3d& centre,
    seissol::memory::AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>&
        mInvJInvPhisAtSources,
    unsigned meshId,
    const seissol::geometry::MeshReader& mesh) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  const double* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[elements[meshId].vertices[v]].coords;
  }
  const auto xiEtaZeta = transformations::tetrahedronGlobalToReference(
      coords[0], coords[1], coords[2], coords[3], centre);
  const auto basisFunctionsAtPoint = basisFunction::SampledBasisFunctions<real>(
      ConvergenceOrder, xiEtaZeta(0), xiEtaZeta(1), xiEtaZeta(2));

  const double volume = MeshTools::volume(elements[meshId], vertices);
  const double jInv = 1.0 / (6.0 * volume);

  kernel::computeMInvJInvPhisAtSources krnl;
  krnl.basisFunctionsAtPoint = basisFunctionsAtPoint.m_data.data();
  krnl.M3inv = init::M3inv::Values;
  krnl.mInvJInvPhisAtSources = mInvJInvPhisAtSources.data();
  krnl.JInv = jInv;
  krnl.execute();
}

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
                                        seissol::memory::Memkind memkind) {
  computeMInvJInvPhisAtSources(centre, pointSources.mInvJInvPhisAtSources[index], meshId, mesh);

  auto& faultBasis = pointSources.tensor[index];
  faultBasis[0] = subfault.tan1(0);
  faultBasis[1] = subfault.tan1(1);
  faultBasis[2] = subfault.tan1(2);
  faultBasis[3] = subfault.tan2(0);
  faultBasis[4] = subfault.tan2(1);
  faultBasis[5] = subfault.tan2(2);
  faultBasis[6] = subfault.normal(0);
  faultBasis[7] = subfault.normal(1);
  faultBasis[8] = subfault.normal(2);

  pointSources.A[index] = subfault.area;
  std::array<double, 81> stiffnessTensor{};
  switch (material->getMaterialType()) {
  case seissol::model::MaterialType::Anisotropic:
    [[fallthrough]];
  case seissol::model::MaterialType::Poroelastic:
    if (subfault.mu != 0) {
      logError() << "There are specific fault parameters for the fault. This is only compatible "
                    "with isotropic (visco)elastic materials.";
    }
    material->getFullStiffnessTensor(stiffnessTensor);
    break;
  default:
    seissol::model::ElasticMaterial em = *dynamic_cast<seissol::model::ElasticMaterial*>(material);
    em.mu = (subfault.mu == 0.0) ? em.mu : subfault.mu;
    em.getFullStiffnessTensor(stiffnessTensor);
    break;
  }
  std::copy(
      stiffnessTensor.begin(), stiffnessTensor.end(), pointSources.stiffnessTensor[index].begin());
  pointSources.onsetTime[index] = subfault.tinit;
  pointSources.samplingInterval[index] = subfault.timestep;
  for (unsigned sr = 0; sr < Offsets().size(); ++sr) {
    std::copy(sliprates[sr].begin() + offsets[sr],
              sliprates[sr].begin() + nextOffsets[sr],
              pointSources.sample[sr].data() + pointSources.sampleOffsets[sr][index]);
    pointSources.sampleOffsets[sr][index + 1] =
        pointSources.sampleOffsets[sr][index] + nextOffsets[sr] - offsets[sr];
  }
}

auto mapClusterToMesh(ClusterMapping& clusterMapping,
                      const unsigned* meshIds,
                      seissol::initializer::LTSTree* ltsTree,
                      seissol::initializer::LTS* lts,
                      seissol::initializer::Lut* ltsLut,
                      seissol::initializer::AllocationPlace place) {
  unsigned clusterSource = 0;
  unsigned mapping = 0;
  while (clusterSource < clusterMapping.sources.size()) {
    const unsigned meshId = meshIds[clusterMapping.sources[clusterSource]];
    unsigned next = clusterSource + 1;
    while (next < clusterMapping.sources.size() &&
           meshIds[clusterMapping.sources[next]] == meshId) {
      ++next;
    }

    for (unsigned ltsId = 0, dup = 0; dup < seissol::initializer::Lut::MaxDuplicates &&
                                      (ltsId = ltsLut->ltsId(lts->dofs.mask, meshId, dup)) !=
                                          std::numeric_limits<unsigned>::max();
         ++dup) {
      clusterMapping.cellToSources[mapping].dofs = &ltsTree->var(lts->dofs, place)[ltsId];
      clusterMapping.cellToSources[mapping].pointSourcesOffset = clusterSource;
      clusterMapping.cellToSources[mapping].numberOfPointSources = next - clusterSource;
      ++mapping;
    }

    clusterSource = next;
  }
  assert(mapping == clusterMapping.cellToSources.size());
}

auto mapPointSourcesToClusters(const unsigned* meshIds,
                               unsigned numberOfSources,
                               seissol::initializer::LTSTree* ltsTree,
                               seissol::initializer::LTS* lts,
                               seissol::initializer::Lut* ltsLut,
                               seissol::memory::Memkind memkind)
    -> std::unordered_map<LayerType, std::vector<ClusterMapping>> {
  auto layerClusterToPointSources =
      std::unordered_map<LayerType, std::vector<std::vector<unsigned>>>{};
  layerClusterToPointSources[Copy].resize(ltsTree->numChildren());
  layerClusterToPointSources[Interior].resize(ltsTree->numChildren());
  auto layerClusterToMeshIds = std::unordered_map<LayerType, std::vector<std::vector<unsigned>>>{};
  layerClusterToMeshIds[Copy].resize(ltsTree->numChildren());
  layerClusterToMeshIds[Interior].resize(ltsTree->numChildren());

  for (unsigned source = 0; source < numberOfSources; ++source) {
    const unsigned meshId = meshIds[source];
    const unsigned cluster = ltsLut->cluster(meshId);
    const LayerType layer = ltsLut->layer(meshId);
    assert(layer != Ghost);
    layerClusterToPointSources[layer][cluster].push_back(source);
    layerClusterToMeshIds[layer][cluster].push_back(meshId);
  }

  std::unordered_map<LayerType, std::vector<ClusterMapping>> layeredClusterMapping;
  for (auto layer : {Copy, Interior}) {
    layeredClusterMapping[layer].resize(ltsTree->numChildren(), ClusterMapping(memkind));
    auto& clusterToMeshIds = layerClusterToMeshIds[layer];
    auto& clusterToPointSources = layerClusterToPointSources[layer];
    auto& clusterMappings = layeredClusterMapping[layer];
    for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
      // Determine number of mappings by counting unique mesh Ids
      std::sort(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
      auto last = std::unique(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
      unsigned numberOfMappings = 0;
      for (auto it = clusterToMeshIds[cluster].begin(); it != last; ++it) {
        const unsigned meshId = *it;
        for (unsigned dup = 0;
             dup < seissol::initializer::Lut::MaxDuplicates &&
             ltsLut->ltsId(lts->dofs.mask, meshId, dup) != std::numeric_limits<unsigned>::max();
             ++dup) {
          ++numberOfMappings;
        }
      }

      clusterMappings[cluster].sources.resize(clusterToPointSources[cluster].size());
      clusterMappings[cluster].cellToSources.resize(numberOfMappings);

      for (unsigned source = 0; source < clusterToPointSources[cluster].size(); ++source) {
        clusterMappings[cluster].sources[source] = clusterToPointSources[cluster][source];
      }
      std::sort(clusterMappings[cluster].sources.begin(),
                clusterMappings[cluster].sources.end(),
                [&](unsigned i, unsigned j) { return meshIds[i] < meshIds[j]; });

      mapClusterToMesh(clusterMappings[cluster],
                       meshIds,
                       ltsTree,
                       lts,
                       ltsLut,
                       seissol::initializer::AllocationPlace::Host);
    }
  }

  return layeredClusterMapping;
}

auto makePointSourceCluster(const ClusterMapping& mapping,
                            const PointSources& sources,
                            const unsigned* meshIds,
                            seissol::initializer::LTSTree* ltsTree,
                            seissol::initializer::LTS* lts,
                            seissol::initializer::Lut* ltsLut)
    -> seissol::kernels::PointSourceClusterPair {
  auto hostData = std::pair<std::shared_ptr<ClusterMapping>, std::shared_ptr<PointSources>>(
      std::make_shared<ClusterMapping>(mapping), std::make_shared<PointSources>(sources));

#if defined(ACL_DEVICE) && !defined(MULTIPLE_SIMULATIONS)
  using GpuImpl = seissol::kernels::PointSourceClusterOnDevice;

  auto deviceData =
      [&]() -> std::pair<std::shared_ptr<ClusterMapping>, std::shared_ptr<PointSources>> {
    if (useUSM()) {
      return hostData;
    } else {
      constexpr auto GpuMemkind = seissol::memory::Memkind::DeviceGlobalMemory;
      auto predeviceClusterMapping = mapping;
      mapClusterToMesh(predeviceClusterMapping,
                       meshIds,
                       ltsTree,
                       lts,
                       ltsLut,
                       seissol::initializer::AllocationPlace::Device);
      auto deviceClusterMapping =
          std::make_shared<ClusterMapping>(predeviceClusterMapping, GpuMemkind);
      auto devicePointSources = std::make_shared<PointSources>(sources, GpuMemkind);
      return {deviceClusterMapping, devicePointSources};
    }
  }();
#else
  using GpuImpl = seissol::kernels::PointSourceClusterOnHost;
  auto deviceData = hostData;
#endif

  return seissol::kernels::PointSourceClusterPair{
      std::make_unique<kernels::PointSourceClusterOnHost>(hostData.first, hostData.second),
      std::make_unique<GpuImpl>(deviceData.first, deviceData.second)};
}

auto loadSourcesFromFSRM(const char* fileName,
                         const seissol::geometry::MeshReader& mesh,
                         seissol::initializer::LTSTree* ltsTree,
                         seissol::initializer::LTS* lts,
                         seissol::initializer::Lut* ltsLut,
                         seissol::memory::Memkind memkind)
    -> std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>> {
  // until further rewrite, we'll leave most of the raw pointers/arrays in here.

  const int rank = seissol::MPI::mpi.rank();

  seissol::sourceterm::FSRMSource fsrm;
  fsrm.read(std::string(fileName));

  logInfo(rank) << "Finding meshIds for point sources...";

  auto contained = std::vector<short>(fsrm.numberOfSources);
  auto meshIds = std::vector<unsigned>(fsrm.numberOfSources);

  initializer::findMeshIds(
      fsrm.centers.data(), mesh, fsrm.numberOfSources, contained.data(), meshIds.data());

#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring point sources for MPI...";
  initializer::cleanDoubles(contained.data(), fsrm.numberOfSources);
#endif

  auto originalIndex = std::vector<unsigned>(fsrm.numberOfSources);
  unsigned numSources = 0;
  for (unsigned source = 0; source < fsrm.numberOfSources; ++source) {
    originalIndex[numSources] = source;
    meshIds[numSources] = meshIds[source];
    numSources += contained[source];
  }

  logInfo(rank) << "Mapping point sources to LTS cells...";
  auto layeredClusterMapping =
      mapPointSourcesToClusters(meshIds.data(), numSources, ltsTree, lts, ltsLut, memkind);
  std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>>
      layeredSourceClusters;

  for (auto layer : {Interior, Copy}) {
    auto& sourceCluster = layeredSourceClusters[layer];
    sourceCluster.resize(ltsTree->numChildren());
    auto& clusterMappings = layeredClusterMapping[layer];
    for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
      auto numberOfSources = clusterMappings[cluster].sources.size();
      auto sources = PointSources{memkind};
      sources.mode = PointSourceMode::Fsrm;
      sources.numberOfSources = numberOfSources;
      sources.mInvJInvPhisAtSources.resize(numberOfSources);
      sources.tensor.resize(numberOfSources);
      sources.onsetTime.resize(numberOfSources);
      sources.samplingInterval.resize(numberOfSources);
      sources.sampleOffsets[0].resize(numberOfSources + 1);
      sources.sampleOffsets[0][0] = 0;
      sources.sample[0].resize(fsrm.numberOfSamples * numberOfSources);

      real* ptrVelocityComponent = fsrm.solidVelocityComponent;

      for (unsigned clusterSource = 0; clusterSource < numberOfSources; ++clusterSource) {
        const unsigned sourceIndex = clusterMappings[cluster].sources[clusterSource];
        const unsigned fsrmIndex = originalIndex[sourceIndex];

        computeMInvJInvPhisAtSources(fsrm.centers[fsrmIndex],
                                     sources.mInvJInvPhisAtSources[clusterSource],
                                     meshIds[sourceIndex],
                                     mesh);
        if (fsrm.solidVel3cVector[fsrmIndex][0] + fsrm.solidVel3cVector[fsrmIndex][1] +
                fsrm.solidVel3cVector[fsrmIndex][2] >
            1e-5) {
          ptrVelocityComponent = fsrm.solidVel3cVector[fsrmIndex].data();
        }
        transformMomentTensor(fsrm.momentTensor,
                              ptrVelocityComponent,
                              fsrm.pressureComponent,
                              fsrm.fluidVelocityComponent,
                              fsrm.strikes[fsrmIndex],
                              fsrm.dips[fsrmIndex],
                              fsrm.rakes[fsrmIndex],
                              sources.tensor[clusterSource]);

        for (unsigned i = 0; i < PointSources::TensorSize; ++i) {
          sources.tensor[clusterSource][i] *= fsrm.areas[fsrmIndex];
        }
#ifndef USE_POROELASTIC
        const seissol::model::Material& material =
            ltsLut->lookup(lts->material, meshIds[sourceIndex] - 1).local;
        for (unsigned i = 0; i < 3; ++i) {
          sources.tensor[clusterSource][6 + i] /= material.rho;
        }
#else
        logWarning() << "The poroelastic equation does not scale the force components with the "
                        "density. For the definition of the sources in poroelastic media, we refer "
                        "to the documentation of SeisSol.";
#endif

        sources.onsetTime[clusterSource] = fsrm.onsets[fsrmIndex];
        sources.samplingInterval[clusterSource] = fsrm.timestep;
        std::copy(std::begin(fsrm.timeHistories[fsrmIndex]),
                  std::end(fsrm.timeHistories[fsrmIndex]),
                  sources.sample[0].data() + sources.sampleOffsets[0][clusterSource]);
        sources.sampleOffsets[0][clusterSource + 1] =
            sources.sampleOffsets[0][clusterSource] + fsrm.timeHistories[fsrmIndex].size();
      }

      sourceCluster[cluster] = makePointSourceCluster(
          clusterMappings[cluster], sources, meshIds.data(), ltsTree, lts, ltsLut);
    }
  }

  logInfo(rank) << ".. finished point source initialization.";

  return layeredSourceClusters;
}

// TODO Add support for passive netCDF
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
auto loadSourcesFromNRF(const char* fileName,
                        const seissol::geometry::MeshReader& mesh,
                        seissol::initializer::LTSTree* ltsTree,
                        seissol::initializer::LTS* lts,
                        seissol::initializer::Lut* ltsLut,
                        seissol::memory::Memkind memkind)
    -> std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>> {
  const int rank = seissol::MPI::mpi.rank();

  logInfo(rank) << "Reading" << fileName;
  NRF nrf;
  readNRF(fileName, nrf);

  auto contained = std::vector<short>(nrf.size());
  auto meshIds = std::vector<unsigned>(nrf.size());

  logInfo(rank) << "Finding meshIds for point sources...";
  initializer::findMeshIds(nrf.centres.data(), mesh, nrf.size(), contained.data(), meshIds.data());

#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring point sources for MPI...";
  initializer::cleanDoubles(contained.data(), nrf.size());
#endif

  auto originalIndex = std::vector<unsigned>(nrf.size());
  unsigned numSources = 0;
  for (unsigned source = 0; source < nrf.size(); ++source) {
    originalIndex[numSources] = source;
    meshIds[numSources] = meshIds[source];
    numSources += contained[source];
  }

  // Checking that all sources are within the domain
  unsigned globalnumSources = numSources;
#ifdef USE_MPI
  MPI_Reduce(&numSources, &globalnumSources, 1, MPI_UNSIGNED, MPI_SUM, 0, seissol::MPI::mpi.comm());
#endif

  if (rank == 0) {
    const int numSourceOutside = nrf.size() - globalnumSources;
    if (numSourceOutside > 0) {
      logError() << nrf.size() - globalnumSources << " point sources are outside the domain.";
    }
  }

  logInfo(rank) << "Mapping point sources to LTS cells...";

  auto layeredClusterMapping =
      mapPointSourcesToClusters(meshIds.data(), numSources, ltsTree, lts, ltsLut, memkind);
  std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>>
      layeredSourceClusters;

  for (auto layer : {Interior, Copy}) {
    auto& sourceCluster = layeredSourceClusters[layer];
    sourceCluster.resize(ltsTree->numChildren());
    auto& clusterMappings = layeredClusterMapping[layer];
    for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
      auto numberOfSources = clusterMappings[cluster].sources.size();
      auto sources = PointSources{memkind};
      sources.mode = PointSourceMode::Nrf;
      sources.numberOfSources = numberOfSources;
      sources.mInvJInvPhisAtSources.resize(numberOfSources);
      sources.tensor.resize(numberOfSources);
      sources.A.resize(numberOfSources);
      sources.stiffnessTensor.resize(numberOfSources);
      sources.onsetTime.resize(numberOfSources);
      sources.samplingInterval.resize(numberOfSources);
      for (auto& so : sources.sampleOffsets) {
        so.resize(numberOfSources + 1);
        so[0] = 0;
      }

      for (std::size_t i = 0; i < Offsets().size(); ++i) {
        std::size_t sampleSize = 0;
        for (unsigned clusterSource = 0; clusterSource < numberOfSources; ++clusterSource) {
          const unsigned sourceIndex = clusterMappings[cluster].sources[clusterSource];
          const unsigned nrfIndex = originalIndex[sourceIndex];
          sampleSize += nrf.sroffsets[nrfIndex + 1][i] - nrf.sroffsets[nrfIndex][i];
        }
        sources.sample[i].resize(sampleSize);
      }

      for (unsigned clusterSource = 0; clusterSource < numberOfSources; ++clusterSource) {
        const unsigned sourceIndex = clusterMappings[cluster].sources[clusterSource];
        const unsigned nrfIndex = originalIndex[sourceIndex];
        transformNRFSourceToInternalSource(
            nrf.centres[nrfIndex],
            meshIds[sourceIndex],
            mesh,
            nrf.subfaults[nrfIndex],
            nrf.sroffsets[nrfIndex],
            nrf.sroffsets[nrfIndex + 1],
            nrf.sliprates,
            &ltsLut->lookup(lts->material, meshIds[sourceIndex]).local,
            sources,
            clusterSource,
            memkind);
      }
      sourceCluster[cluster] = makePointSourceCluster(
          clusterMappings[cluster], sources, meshIds.data(), ltsTree, lts, ltsLut);
    }
  }

  logInfo(rank) << ".. finished point source initialization.";

  return layeredSourceClusters;
}
#endif // defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)

} // namespace

namespace seissol::sourceterm {

void Manager::loadSources(seissol::initializer::parameters::PointSourceType sourceType,
                          const char* fileName,
                          const seissol::geometry::MeshReader& mesh,
                          seissol::initializer::LTSTree* ltsTree,
                          seissol::initializer::LTS* lts,
                          seissol::initializer::Lut* ltsLut,
                          time_stepping::TimeManager& timeManager) {
#ifdef ACL_DEVICE
  auto memkind = useUSM() ? seissol::memory::DeviceUnifiedMemory : seissol::memory::Standard;
#else
  auto memkind = seissol::memory::Standard;
#endif
  auto sourceClusters =
      std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>>{};
  if (sourceType == seissol::initializer::parameters::PointSourceType::NrfSource) {
    logInfo(seissol::MPI::mpi.rank()) << "Reading an NRF source (type 42).";
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
    sourceClusters = loadSourcesFromNRF(fileName, mesh, ltsTree, lts, ltsLut, memkind);
#else
    logError() << "NRF sources (type 42) need SeisSol to be linked with an (active) Netcdf "
                  "library. However, this is not the case for this build.";
#endif
  } else if (sourceType == seissol::initializer::parameters::PointSourceType::FsrmSource) {
    logInfo(seissol::MPI::mpi.rank()) << "Reading an FSRM source (type 50).";
    sourceClusters = loadSourcesFromFSRM(fileName, mesh, ltsTree, lts, ltsLut, memkind);
  } else if (sourceType == seissol::initializer::parameters::PointSourceType::None) {
    logInfo(seissol::MPI::mpi.rank()) << "No source term specified.";
  } else {
    logError() << "The source type" << static_cast<int>(sourceType)
               << "has been defined, but not yet been implemented in SeisSol.";
  }
  // otherwise, we do not have any source term.

  timeManager.setPointSourcesForClusters(std::move(sourceClusters));
}

} // namespace seissol::sourceterm
