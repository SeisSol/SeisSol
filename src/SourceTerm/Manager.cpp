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

#include "Manager.h"
#include "FSRMReader.h"
#include "NRFReader.h"
#include "PointSource.h"
#include "Numerical_aux/Transformation.h"
#include "generated_code/kernel.h"
#include "generated_code/init.h"
#include "generated_code/tensor.h"
#include "Parallel/MPI.h"

#include <Initializer/PointMapper.h>
#include <Kernels/PointSourceClusterOnHost.h>
#include <utils/logger.h>
#include <string>
#include <cstring>

#ifdef ACL_DEVICE
#include <Kernels/PointSourceClusterOnDevice.h>
#include <Parallel/AcceleratorDevice.h>
#include "Device/UsmAllocator.h"
#endif

/**
 * Computes mInvJInvPhisAtSources[i] = |J|^-1 * M_ii^-1 * phi_i(xi, eta, zeta),
 * where xi, eta, zeta is the point in the reference tetrahedron corresponding to x, y, z.
 */
void seissol::sourceterm::computeMInvJInvPhisAtSources(
    Eigen::Vector3d const& centre,
    AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>& mInvJInvPhisAtSources,
    unsigned meshId,
    seissol::geometry::MeshReader const& mesh) {
  auto const& elements = mesh.getElements();
  auto const& vertices = mesh.getVertices();

  double const* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[elements[meshId].vertices[v]].coords;
  }
  auto const xiEtaZeta = transformations::tetrahedronGlobalToReference(
      coords[0], coords[1], coords[2], coords[3], centre);
  auto const basisFunctionsAtPoint = basisFunction::SampledBasisFunctions<real>(
      CONVERGENCE_ORDER, xiEtaZeta(0), xiEtaZeta(1), xiEtaZeta(2));

  double volume = MeshTools::volume(elements[meshId], vertices);
  double JInv = 1.0 / (6.0 * volume);

  kernel::computeMInvJInvPhisAtSources krnl;
  krnl.basisFunctionsAtPoint = basisFunctionsAtPoint.m_data.data();
  krnl.M3inv = init::M3inv::Values;
  krnl.mInvJInvPhisAtSources = mInvJInvPhisAtSources.data();
  krnl.JInv = JInv;
  krnl.execute();
}

void seissol::sourceterm::transformNRFSourceToInternalSource(
    Eigen::Vector3d const& centre,
    unsigned meshId,
    seissol::geometry::MeshReader const& mesh,
    Subfault const& subfault,
    Offsets const& offsets,
    Offsets const& nextOffsets,
    std::array<std::vector<double>, 3> const& sliprates,
    seissol::model::Material* material,
    PointSources& pointSources,
    unsigned index,
    AllocatorT const& alloc) {
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
  switch (material->getMaterialType()) {
  case seissol::model::MaterialType::anisotropic:
    [[fallthrough]];
  case seissol::model::MaterialType::poroelastic:
    if (subfault.mu != 0) {
      logError() << "There are specific fault parameters for the fault. This is only compatible "
                    "with isotropic (visco)elastic materials.";
    }
    material->getFullStiffnessTensor(pointSources.stiffnessTensor[index]);
    break;
  default:
    seissol::model::ElasticMaterial em = *dynamic_cast<seissol::model::ElasticMaterial*>(material);
    em.mu = (subfault.mu == 0.0) ? em.mu : subfault.mu;
    em.getFullStiffnessTensor(pointSources.stiffnessTensor[index]);
    break;
  }

  for (unsigned sr = 0; sr < pointSources.slipRates.size(); ++sr) {
    unsigned numSamples = nextOffsets[sr] - offsets[sr];
    double const* samples = (numSamples > 0) ? &sliprates[sr][offsets[sr]] : NULL;
    pointSources.slipRates[sr][index] =
        PiecewiseLinearFunction1D(samples, numSamples, subfault.tinit, subfault.timestep, alloc);
  }
}

auto seissol::sourceterm::Manager::mapPointSourcesToClusters(
    const unsigned* meshIds,
    unsigned numberOfSources,
    seissol::initializers::LTSTree* ltsTree,
    seissol::initializers::LTS* lts,
    seissol::initializers::Lut* ltsLut,
    AllocatorT const& alloc) -> std::unordered_map<LayerType, std::vector<ClusterMapping>> {
  auto layerClusterToPointSources =
      std::unordered_map<LayerType, std::vector<std::vector<unsigned>>>{};
  layerClusterToPointSources[Copy].resize(ltsTree->numChildren());
  layerClusterToPointSources[Interior].resize(ltsTree->numChildren());
  auto layerClusterToMeshIds = std::unordered_map<LayerType, std::vector<std::vector<unsigned>>>{};
  layerClusterToMeshIds[Copy].resize(ltsTree->numChildren());
  layerClusterToMeshIds[Interior].resize(ltsTree->numChildren());

  for (unsigned source = 0; source < numberOfSources; ++source) {
    unsigned meshId = meshIds[source];
    unsigned cluster = ltsLut->cluster(meshId);
    LayerType layer = ltsLut->layer(meshId);
    assert(layer != Ghost);
    layerClusterToPointSources[layer][cluster].push_back(source);
    layerClusterToMeshIds[layer][cluster].push_back(meshId);
  }

  std::unordered_map<LayerType, std::vector<ClusterMapping>> layeredClusterMapping;
  for (auto layer : {Copy, Interior}) {
    layeredClusterMapping[layer].resize(ltsTree->numChildren(), ClusterMapping(alloc));
    auto& clusterToMeshIds = layerClusterToMeshIds[layer];
    auto& clusterToPointSources = layerClusterToPointSources[layer];
    auto& clusterMappings = layeredClusterMapping[layer];
    for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
      // Determine number of mappings by counting unique mesh Ids
      std::sort(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
      auto last = std::unique(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
      unsigned numberOfMappings = 0;
      for (auto it = clusterToMeshIds[cluster].begin(); it != last; ++it) {
        unsigned meshId = *it;
        for (unsigned dup = 0;
             dup < seissol::initializers::Lut::MaxDuplicates &&
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

      unsigned clusterSource = 0;
      unsigned mapping = 0;
      while (clusterSource < clusterMappings[cluster].sources.size()) {
        unsigned meshId = meshIds[clusterMappings[cluster].sources[clusterSource]];
        unsigned next = clusterSource + 1;
        while (next < clusterMappings[cluster].sources.size() &&
               meshIds[clusterMappings[cluster].sources[next]] == meshId) {
          ++next;
        }

        for (unsigned ltsId, dup = 0; dup < seissol::initializers::Lut::MaxDuplicates &&
                                      (ltsId = ltsLut->ltsId(lts->dofs.mask, meshId, dup)) !=
                                          std::numeric_limits<unsigned>::max();
             ++dup) {
          clusterMappings[cluster].cellToSources[mapping].dofs = &ltsTree->var(lts->dofs)[ltsId];
          clusterMappings[cluster].cellToSources[mapping].pointSourcesOffset = clusterSource;
          clusterMappings[cluster].cellToSources[mapping].numberOfPointSources =
              next - clusterSource;
          ++mapping;
        }

        clusterSource = next;
      }
      assert(mapping == clusterMappings[cluster].cellToSources.size());
    }
  }

  return layeredClusterMapping;
}

void seissol::sourceterm::Manager::loadSources(SourceType sourceType,
                                               char const* fileName,
                                               seissol::geometry::MeshReader const& mesh,
                                               seissol::initializers::LTSTree* ltsTree,
                                               seissol::initializers::LTS* lts,
                                               seissol::initializers::Lut* ltsLut,
                                               time_stepping::TimeManager& timeManager) {
#ifdef ACL_DEVICE
  auto& instance = device::DeviceInstance::getInstance();
  auto alloc = device::UsmAllocator<real>(instance);
#else
  auto alloc = AllocatorT();
#endif
  auto sourceClusters =
      std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>{};
  if (sourceType == SourceType::NrfSource) {
    logInfo(seissol::MPI::mpi.rank()) << "Reading an NRF source (type 42).";
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
    sourceClusters = loadSourcesFromNRF(fileName, mesh, ltsTree, lts, ltsLut, alloc);
#else
    logError() << "NRF sources (type 42) need SeisSol to be linked with an (active) Netcdf "
                  "library. However, this is not the case for this build.";
#endif
  } else if (sourceType == SourceType::FsrmSource) {
    logInfo(seissol::MPI::mpi.rank()) << "Reading an FSRM source (type 50).";
    sourceClusters = loadSourcesFromFSRM(fileName, mesh, ltsTree, lts, ltsLut, alloc);
  } else if (sourceType == SourceType::None) {
    logInfo(seissol::MPI::mpi.rank()) << "No source term specified.";
  } else {
    logError() << "The source type" << static_cast<int>(sourceType)
               << "has been defined, but not yet been implemented in SeisSol.";
  }
  // otherwise, we do not have any source term.

  timeManager.setPointSourcesForClusters(std::move(sourceClusters));
}

auto seissol::sourceterm::Manager::makePointSourceCluster(ClusterMapping mapping,
                                                          PointSources sources)
    -> std::unique_ptr<kernels::PointSourceCluster> {
#if defined(ACL_DEVICE) && !defined(MULTIPLE_SIMULATIONS)
  using Impl = kernels::PointSourceClusterOnDevice;
#else
  using Impl = kernels::PointSourceClusterOnHost;
#endif
  return std::make_unique<Impl>(std::move(mapping), std::move(sources));
}

auto seissol::sourceterm::Manager::loadSourcesFromFSRM(char const* fileName,
                                                       seissol::geometry::MeshReader const& mesh,
                                                       seissol::initializers::LTSTree* ltsTree,
                                                       seissol::initializers::LTS* lts,
                                                       seissol::initializers::Lut* ltsLut,
                                                       AllocatorT const& alloc)
    -> std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>> {
  // until further rewrite, we'll leave most of the raw pointers/arrays in here.

  int rank = seissol::MPI::mpi.rank();

  seissol::sourceterm::FSRMSource fsrm;
  fsrm.read(std::string(fileName));

  logInfo(rank) << "Finding meshIds for point sources...";

  auto contained = std::vector<short>(fsrm.numberOfSources);
  auto meshIds = std::vector<unsigned>(fsrm.numberOfSources);

  initializers::findMeshIds(
      fsrm.centers.data(), mesh, fsrm.numberOfSources, contained.data(), meshIds.data());

#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring point sources for MPI...";
  initializers::cleanDoubles(contained.data(), fsrm.numberOfSources);
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
      mapPointSourcesToClusters(meshIds.data(), numSources, ltsTree, lts, ltsLut, alloc);
  std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>
      layeredSourceClusters;

  for (auto layer : {Interior, Copy}) {
    auto& sourceCluster = layeredSourceClusters[layer];
    sourceCluster.resize(ltsTree->numChildren());
    auto& clusterMappings = layeredClusterMapping[layer];
    for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
      auto numberOfSources = clusterMappings[cluster].sources.size();
      auto sources = PointSources{alloc};
      sources.mode = PointSources::FSRM;
      sources.numberOfSources = numberOfSources;
      sources.mInvJInvPhisAtSources.resize(numberOfSources);
      sources.tensor.resize(numberOfSources);
      sources.slipRates[0].resize(numberOfSources, PiecewiseLinearFunction1D{alloc});

      for (unsigned clusterSource = 0; clusterSource < numberOfSources; ++clusterSource) {
        unsigned sourceIndex = clusterMappings[cluster].sources[clusterSource];
        unsigned fsrmIndex = originalIndex[sourceIndex];

        computeMInvJInvPhisAtSources(fsrm.centers[fsrmIndex],
                                     sources.mInvJInvPhisAtSources[clusterSource],
                                     meshIds[sourceIndex],
                                     mesh);
        transformMomentTensor(fsrm.momentTensor,
                              fsrm.solidVelocityComponent,
                              fsrm.pressureComponent,
                              fsrm.fluidVelocityComponent,
                              fsrm.strikes[fsrmIndex],
                              fsrm.dips[fsrmIndex],
                              fsrm.rakes[fsrmIndex],
                              sources.tensor[clusterSource]);

        for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
          sources.tensor[clusterSource][i] *= fsrm.areas[fsrmIndex];
        }
#ifndef USE_POROELASTIC
        seissol::model::Material& material =
            ltsLut->lookup(lts->material, meshIds[sourceIndex] - 1).local;
        for (unsigned i = 0; i < 3; ++i) {
          sources.tensor[clusterSource][6 + i] /= material.rho;
        }
#else
        logWarning() << "The poroelastic equation does not scale the force components with the "
                        "density. For the definition of the sources in poroelastic media, we refer "
                        "to the documentation of SeisSol.";
#endif

        sources.slipRates[0][clusterSource] =
            PiecewiseLinearFunction1D(fsrm.timeHistories[fsrmIndex].data(),
                                      fsrm.numberOfSamples,
                                      fsrm.onsets[fsrmIndex],
                                      fsrm.timestep,
                                      alloc);
      }

      sourceCluster[cluster] =
          makePointSourceCluster(std::move(clusterMappings[cluster]), std::move(sources));
    }
  }

  logInfo(rank) << ".. finished point source initialization.";

  return layeredSourceClusters;
}

// TODO Add support for passive netCDF
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
auto seissol::sourceterm::Manager::loadSourcesFromNRF(char const* fileName,
                                                      seissol::geometry::MeshReader const& mesh,
                                                      seissol::initializers::LTSTree* ltsTree,
                                                      seissol::initializers::LTS* lts,
                                                      seissol::initializers::Lut* ltsLut,
                                                      AllocatorT const& alloc)
    -> std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>> {
  int rank = seissol::MPI::mpi.rank();

  logInfo(rank) << "<--------------------------------------------------------->";
  logInfo(rank) << "<                      Point sources                      >";
  logInfo(rank) << "<--------------------------------------------------------->";

  logInfo(rank) << "Reading" << fileName;
  NRF nrf;
  readNRF(fileName, nrf);

  auto contained = std::vector<short>(nrf.size());
  auto meshIds = std::vector<unsigned>(nrf.size());

  logInfo(rank) << "Finding meshIds for point sources...";
  initializers::findMeshIds(nrf.centres.data(), mesh, nrf.size(), contained.data(), meshIds.data());

#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring point sources for MPI...";
  initializers::cleanDoubles(contained.data(), nrf.size());
#endif

  auto originalIndex = std::vector<unsigned>(nrf.size());
  unsigned numSources = 0;
  for (unsigned source = 0; source < nrf.size(); ++source) {
    originalIndex[numSources] = source;
    meshIds[numSources] = meshIds[source];
    numSources += contained[source];
  }

  // Checking that all sources are within the domain
  int globalnumSources = numSources;
#ifdef USE_MPI
  MPI_Reduce(&numSources, &globalnumSources, 1, MPI_INT, MPI_SUM, 0, seissol::MPI::mpi.comm());
#endif

  if (rank == 0) {
    int numSourceOutside = nrf.size() - globalnumSources;
    if (numSourceOutside > 0) {
      logError() << nrf.size() - globalnumSources << " point sources are outside the domain.";
    }
  }

  logInfo(rank) << "Mapping point sources to LTS cells...";

  auto layeredClusterMapping =
      mapPointSourcesToClusters(meshIds.data(), numSources, ltsTree, lts, ltsLut, alloc);
  std::unordered_map<LayerType, std::vector<std::unique_ptr<kernels::PointSourceCluster>>>
      layeredSourceClusters;

  for (auto layer : {Interior, Copy}) {
    auto& sourceCluster = layeredSourceClusters[layer];
    sourceCluster.resize(ltsTree->numChildren());
    auto& clusterMappings = layeredClusterMapping[layer];
    for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
      auto numberOfSources = clusterMappings[cluster].sources.size();
      auto sources = PointSources{alloc};
      sources.mode = PointSources::NRF;
      sources.numberOfSources = numberOfSources;
      sources.mInvJInvPhisAtSources.resize(numberOfSources);
      sources.tensor.resize(numberOfSources);
      sources.A.resize(numberOfSources);
      sources.stiffnessTensor.resize(numberOfSources);
      for (auto& sr : sources.slipRates) {
        sr.resize(numberOfSources, PiecewiseLinearFunction1D{alloc});
      }

      for (unsigned clusterSource = 0; clusterSource < numberOfSources; ++clusterSource) {
        unsigned sourceIndex = clusterMappings[cluster].sources[clusterSource];
        unsigned nrfIndex = originalIndex[sourceIndex];
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
            alloc);
      }
      sourceCluster[cluster] =
          makePointSourceCluster(std::move(clusterMappings[cluster]), std::move(sources));
    }
  }

  logInfo(rank) << ".. finished point source initialization.";

  return layeredSourceClusters;
}
#endif // defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
