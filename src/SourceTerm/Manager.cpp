// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)
// SPDX-FileContributor: Sebastian Rettenberger

#include "Manager.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Equations/Datastructures.h"
#include "FSRMReader.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Initializer/Parameters/SourceParameters.h"
#include "Initializer/PointMapper.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/PointSourceClusterOnHost.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Model/Datastructures.h" // IWYU pragma: keep
#include "Numerical/BasisFunction.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "PointSource.h"
#include "Solver/MultipleSimulations.h"
#include "Solver/TimeStepping/TimeManager.h"
#include "SourceTerm/Typedefs.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef USE_NETCDF
#include "NRFReader.h"
#include "SourceTerm/NRF.h"

#include <mpi.h>
#endif

#ifdef ACL_DEVICE
#include "Kernels/PointSourceClusterOnDevice.h"
#endif

namespace seissol::sourceterm {

namespace {

/**
 * Computes mInvJInvPhisAtSources[i] = |J|^-1 * M_ii^-1 * phi_i(xi, eta, zeta),
 * where xi, eta, zeta is the point in the reference tetrahedron corresponding to x, y, z.
 */
void computeMInvJInvPhisAtSources(
    const Eigen::Vector3d& center,
    seissol::memory::AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>&
        mInvJInvPhisAtSources,
    std::size_t meshId,
    const seissol::geometry::MeshReader& mesh) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  const auto transform = seissol::geometry::AffineTransform::fromMeshCell(meshId, mesh);
  const auto xiEtaZeta = transform.spaceToRef(center);
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

struct SourceFile {
  std::vector<std::size_t> originalIndex;
  std::vector<std::size_t> meshIds;
};

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
void transformNRFSourceToInternalSource(const Subfault& subfault,
                                        const Offsets& offsets,
                                        const Offsets& nextOffsets,
                                        const std::array<std::vector<double>, 3>& sliprates,
                                        const seissol::model::Material* material,
                                        PointSources& pointSources,
                                        std::size_t index,
                                        std::size_t tensorIndex) {
  std::array<real, 9> faultBasis{};
  faultBasis[0] = subfault.tan1(0);
  faultBasis[1] = subfault.tan1(1);
  faultBasis[2] = subfault.tan1(2);
  faultBasis[3] = subfault.tan2(0);
  faultBasis[4] = subfault.tan2(1);
  faultBasis[5] = subfault.tan2(2);
  faultBasis[6] = subfault.normal(0);
  faultBasis[7] = subfault.normal(1);
  faultBasis[8] = subfault.normal(2);

  std::array<double, 81> stiffnessTensor{};
  switch (material->getMaterialType()) {
  case seissol::model::MaterialType::Anisotropic:
    [[fallthrough]];
  case seissol::model::MaterialType::Acoustic:
    logError() << "NRF sources are only compatible with isotropic (visco)elastic and "
                  "poroelastic materials.";
    break;
  case seissol::model::MaterialType::Poroelastic:
    if (subfault.mu != 0) {
      logError() << "There are specific fault parameters for the fault. This is only compatible "
                    "with isotropic (visco)elastic materials.";
    }
    material->getFullStiffnessTensor(stiffnessTensor);
    break;
  default:
    seissol::model::ElasticMaterial em =
        *dynamic_cast<const seissol::model::ElasticMaterial*>(material);
    em.mu = (subfault.mu == 0.0) ? em.mu : subfault.mu;
    em.getFullStiffnessTensor(stiffnessTensor);
    break;
  }

  std::array<real, 81> stiffnessTensorReal{};
  std::copy(stiffnessTensor.begin(), stiffnessTensor.end(), stiffnessTensorReal.begin());

  kernel::transformNRF transformKernel;
  transformKernel.mArea = -subfault.area;
  transformKernel.mNormal = faultBasis.data() + 6;
  transformKernel.stiffnessTensor = stiffnessTensorReal.data();
  transformKernel.momentToNRF = init::momentToNRF::Values;
  transformKernel.rotateNRF = faultBasis.data();
  transformKernel.tensorNRF = pointSources.tensor.data() + tensorIndex * tensor::update::Size;

  transformKernel.execute();

  pointSources.onsetTime[index] = subfault.tinit;
  pointSources.samplingInterval[index] = subfault.timestep;
  for (std::size_t sr = 0; sr < Offsets().size(); ++sr) {
    std::copy(sliprates[sr].begin() + offsets[sr],
              sliprates[sr].begin() + nextOffsets[sr],
              pointSources.sample.data() + pointSources.sampleOffsets[tensorIndex + sr]);
    pointSources.sampleOffsets[tensorIndex + sr + 1] =
        pointSources.sampleOffsets[tensorIndex + sr] + nextOffsets[sr] - offsets[sr];
  }
}

struct NrfFile : public SourceFile {
  NRF nrf;
  void read(const std::string& file) { readNRF(file.c_str(), nrf); }

  [[nodiscard]] const std::vector<Eigen::Vector3d>& points() const { return nrf.centres; }

  [[nodiscard]] std::size_t dataSources(std::size_t /*sourceIndex*/) const { return 3; }

  [[nodiscard]] std::size_t sampleCount(std::size_t sourceIndex) const {
    const std::size_t nrfIndex = originalIndex[sourceIndex];
    std::size_t sampleSize = 0;
    for (std::size_t i = 0; i < Offsets().size(); ++i) {
      sampleSize += nrf.sroffsets[nrfIndex + 1][i] - nrf.sroffsets[nrfIndex][i];
    }
    return sampleSize;
  }

  void transform(PointSources& sources,
                 std::size_t sourceIndex,
                 std::size_t index,
                 const seissol::model::Material& material) {
    const std::size_t nrfIndex = originalIndex[sourceIndex];
    transformNRFSourceToInternalSource(nrf.subfaults[nrfIndex],
                                       nrf.sroffsets[nrfIndex],
                                       nrf.sroffsets[nrfIndex + 1],
                                       nrf.sliprates,
                                       &material,
                                       sources,
                                       index,
                                       sources.sampleRange[index]);
  }
};
#endif // defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)

struct FsrmFile : public SourceFile {
  FSRMSource fsrm;
  void read(const std::string& file) { fsrm.read(file); }

  [[nodiscard]] const std::vector<Eigen::Vector3d>& points() const { return fsrm.centers; }

  [[nodiscard]] std::size_t dataSources(std::size_t /*sourceIndex*/) const { return 1; }

  [[nodiscard]] std::size_t sampleCount(std::size_t /*sourceIndex*/) const {
    return fsrm.numberOfSamples;
  }

  void transform(PointSources& sources,
                 std::size_t sourceIndex,
                 std::size_t index,
                 const seissol::model::Material& material) {
    const std::size_t fsrmIndex = originalIndex[sourceIndex];

    auto* tensor = sources.tensor.data() + sources.sampleRange[index] * tensor::update::Size;
    transformMomentTensor(fsrm.momentTensor,
                          fsrm.solidVelocityComponent,
                          fsrm.pressureComponent,
                          fsrm.fluidVelocityComponent,
                          fsrm.strikes[fsrmIndex],
                          fsrm.dips[fsrmIndex],
                          fsrm.rakes[fsrmIndex],
                          tensor);

    for (std::size_t i = 0; i < tensor::update::Size; ++i) {
      tensor[i] *= fsrm.areas[fsrmIndex];
    }
    if (model::MaterialT::Type != model::MaterialType::Poroelastic) {
      for (std::size_t i = 0; i < Cell::Dim; ++i) {
        tensor[model::MaterialT::TractionQuantities + i] /= material.rho;
      }
    } else {
      logWarning() << "The poroelastic equation does not scale the force components with the "
                      "density. For the definition of the sources in poroelastic media, we refer "
                      "to the documentation of SeisSol.";
    }

    sources.onsetTime[index] = fsrm.onsets[fsrmIndex];
    sources.samplingInterval[index] = fsrm.timestep;
    std::copy(std::begin(fsrm.timeHistories[fsrmIndex]),
              std::end(fsrm.timeHistories[fsrmIndex]),
              sources.sample.data() + sources.sampleOffsets[index]);
    sources.sampleOffsets[index + 1] =
        sources.sampleOffsets[index] + fsrm.timeHistories[fsrmIndex].size();
  }
};

auto mapClusterToMesh(ClusterMapping& clusterMapping,
                      const std::size_t* meshIds,
                      LTS::Storage& ltsStorage,
                      LTS::Backmap& backmap,
                      seissol::initializer::AllocationPlace place) {
  std::size_t clusterSource = 0;
  std::size_t mapping = 0;
  while (clusterSource < clusterMapping.sources.size()) {
    const std::size_t meshId = meshIds[clusterMapping.sources[clusterSource]];
    std::size_t next = clusterSource + 1;
    while (next < clusterMapping.sources.size() &&
           meshIds[clusterMapping.sources[next]] == meshId) {
      ++next;
    }

    for (std::size_t dup = 0; dup < LTS::Backmap::MaxDuplicates; ++dup) {
      const auto position = backmap.getDup(meshId, dup);
      if (position.has_value()) {
        clusterMapping.cellToSources[mapping].dofs =
            ltsStorage.lookup<LTS::Dofs>(position.value(), place);
        clusterMapping.cellToSources[mapping].pointSourcesOffset = clusterSource;
        clusterMapping.cellToSources[mapping].numberOfPointSources = next - clusterSource;
        ++mapping;
      }
    }

    clusterSource = next;
  }
  assert(mapping == clusterMapping.cellToSources.size());
}

auto mapPointSourcesToClusters(const std::size_t* meshIds,
                               std::size_t numberOfSources,
                               LTS::Storage& ltsStorage,
                               LTS::Backmap& backmap,
                               seissol::memory::Memkind memkind) -> std::vector<ClusterMapping> {
  auto clusterToPointSources = std::vector<std::vector<std::size_t>>(ltsStorage.numChildren());
  auto clusterToMeshIds = std::vector<std::vector<std::size_t>>(ltsStorage.numChildren());

  for (std::size_t source = 0; source < numberOfSources; ++source) {
    const auto meshId = meshIds[source];
    const auto position = backmap.get(meshId);
    const auto id = position.color;
    clusterToPointSources[id].push_back(source);
    clusterToMeshIds[id].push_back(meshId);
  }

  std::vector<ClusterMapping> clusterMappings(ltsStorage.numChildren(), ClusterMapping(memkind));
  for (std::size_t cluster = 0; cluster < ltsStorage.numChildren(); ++cluster) {
    // Determine number of mappings by counting unique mesh Ids
    std::sort(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
    auto last = std::unique(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
    std::size_t numberOfMappings = 0;
    for (auto it = clusterToMeshIds[cluster].begin(); it != last; ++it) {
      const auto meshId = *it;

      for (std::size_t dup = 0; dup < LTS::Backmap::MaxDuplicates; ++dup) {
        const auto position = backmap.getDup(meshId, dup);
        if (position.has_value()) {
          ++numberOfMappings;
        }
      }
    }

    clusterMappings[cluster].sources.resize(clusterToPointSources[cluster].size());
    clusterMappings[cluster].cellToSources.resize(numberOfMappings);

    for (std::size_t source = 0; source < clusterToPointSources[cluster].size(); ++source) {
      clusterMappings[cluster].sources[source] = clusterToPointSources[cluster][source];
    }
    std::sort(clusterMappings[cluster].sources.begin(),
              clusterMappings[cluster].sources.end(),
              [&](std::size_t i, std::size_t j) { return meshIds[i] < meshIds[j]; });

    mapClusterToMesh(clusterMappings[cluster],
                     meshIds,
                     ltsStorage,
                     backmap,
                     seissol::initializer::AllocationPlace::Host);
  }

  return clusterMappings;
}

auto makePointSourceCluster(const ClusterMapping& mapping,
                            const PointSources& sources,
                            SEISSOL_GPU_PARAM const std::size_t* meshIds,
                            SEISSOL_GPU_PARAM LTS::Storage& ltsStorage,
                            SEISSOL_GPU_PARAM LTS::Backmap& backmap)
    -> seissol::kernels::PointSourceClusterPair {
  auto hostData = std::pair<std::shared_ptr<ClusterMapping>, std::shared_ptr<PointSources>>(
      std::make_shared<ClusterMapping>(mapping), std::make_shared<PointSources>(sources));

#ifdef ACL_DEVICE
  using GpuImpl = seissol::kernels::PointSourceClusterOnDevice;

  const auto deviceData =
      [&]() -> std::pair<std::shared_ptr<ClusterMapping>, std::shared_ptr<PointSources>> {
    if (useUSM()) {
      return hostData;
    } else {
      constexpr auto GpuMemkind = seissol::memory::Memkind::DeviceGlobalMemory;
      auto predeviceClusterMapping = mapping;
      mapClusterToMesh(predeviceClusterMapping,
                       meshIds,
                       ltsStorage,
                       backmap,
                       seissol::initializer::AllocationPlace::Device);
      const auto deviceClusterMapping =
          std::make_shared<ClusterMapping>(predeviceClusterMapping, GpuMemkind);
      const auto devicePointSources = std::make_shared<PointSources>(sources, GpuMemkind);
      return {deviceClusterMapping, devicePointSources};
    }
  }();
#else
  using GpuImpl = seissol::kernels::PointSourceClusterOnHost;
  const auto deviceData = hostData;
#endif

  return seissol::kernels::PointSourceClusterPair{
      std::make_unique<kernels::PointSourceClusterOnHost>(hostData.first, hostData.second),
      std::make_unique<GpuImpl>(deviceData.first, deviceData.second)};
}

template <typename SourceFileT>
auto loadSourceFile(const char* fileName,
                    const seissol::geometry::MeshReader& mesh,
                    LTS::Storage& ltsStorage,
                    LTS::Backmap& backmap,
                    seissol::memory::Memkind memkind)
    -> std::vector<seissol::kernels::PointSourceClusterPair> {

  logInfo() << "Reading point source file" << fileName << "...";

  SourceFileT file;
  file.read(fileName);

  logInfo() << "Finding mesh IDs for point sources...";

  const auto points = file.points();

  auto contained = std::vector<short>(points.size());
  auto meshIds = std::vector<std::size_t>(points.size());

  initializer::findMeshIds(points.data(), mesh, points.size(), contained.data(), meshIds.data());

  logInfo() << "Cleaning possible double occurring point sources in multi-rank setups...";
  initializer::cleanDoubles(contained.data(), points.size());

  auto originalIndex = std::vector<std::size_t>(points.size());
  std::size_t numSources = 0;
  for (std::size_t source = 0; source < points.size(); ++source) {
    originalIndex[numSources] = source;
    meshIds[numSources] = meshIds[source];
    numSources += contained[source] != 0 ? 1 : 0;
  }

  file.originalIndex = originalIndex;
  file.meshIds = meshIds;

  // Checking that all sources are within the domain
  std::size_t globalnumSources = numSources;
  MPI_Reduce(&numSources,
             &globalnumSources,
             1,
             Mpi::castToMpiType<std::size_t>(),
             MPI_SUM,
             0,
             seissol::Mpi::mpi.comm());

  logInfo() << "Found" << globalnumSources << "point sources.";
  const int rank = seissol::Mpi::mpi.rank();
  if (rank == 0 && points.size() > globalnumSources) {
    logError() << (points.size() - globalnumSources) << " point sources are outside the domain.";
  }

  logInfo() << "Mapping point sources to LTS cells...";
  auto clusterMappings =
      mapPointSourcesToClusters(meshIds.data(), numSources, ltsStorage, backmap, memkind);
  std::vector<seissol::kernels::PointSourceClusterPair> sourceCluster(ltsStorage.numChildren());

  for (std::size_t cluster = 0; cluster < ltsStorage.numChildren(); ++cluster) {
    auto numberOfSources = clusterMappings[cluster].sources.size();

    auto sources = PointSources{memkind};
    sources.numberOfSources = numberOfSources;
    sources.mInvJInvPhisAtSources.resize(numberOfSources);

    std::size_t dataSourceCount = 0;
    std::size_t sampleCount = 0;

    for (std::size_t i = 0; i < numberOfSources; ++i) {
      const auto sourceIndex = clusterMappings[cluster].sources[i];
      dataSourceCount += file.dataSources(sourceIndex);
      sampleCount += file.sampleCount(sourceIndex);
    }

    sources.tensor.resize(dataSourceCount * tensor::update::Size);
    sources.onsetTime.resize(numberOfSources);
    sources.samplingInterval.resize(numberOfSources);
    sources.sampleOffsets.resize(dataSourceCount + 1);
    sources.sampleOffsets[0] = 0;
    sources.sample.resize(sampleCount);

    sources.sampleRange.resize(numberOfSources + 1);
    sources.sampleRange[0] = 0;
    for (std::size_t i = 0; i < numberOfSources; ++i) {
      const auto sourceIndex = clusterMappings[cluster].sources[i];
      sources.sampleRange[i + 1] = sources.sampleRange[i] + file.dataSources(sourceIndex);
    }

    sources.simulationIndex.resize(numberOfSources);

    for (std::size_t clusterSource = 0; clusterSource < numberOfSources; ++clusterSource) {
      const std::size_t sourceIndex = clusterMappings[cluster].sources[clusterSource];
      const auto fileIndex = originalIndex[sourceIndex];

      sources.simulationIndex[clusterSource] = fileIndex % multisim::NumSimulations;
      computeMInvJInvPhisAtSources(points[fileIndex],
                                   sources.mInvJInvPhisAtSources[clusterSource],
                                   meshIds[sourceIndex],
                                   mesh);

      const auto position = backmap.get(meshIds[sourceIndex]);
      const auto& material = *ltsStorage.lookup<LTS::Material>(position).local;
      file.transform(sources, sourceIndex, clusterSource, material);
    }

    sourceCluster[cluster] = makePointSourceCluster(
        clusterMappings[cluster], sources, meshIds.data(), ltsStorage, backmap);
  }

  logInfo() << ".. finished point source initialization.";

  return sourceCluster;
}

} // namespace

void Manager::loadSources(seissol::initializer::parameters::PointSourceType sourceType,
                          const char* fileName,
                          const seissol::geometry::MeshReader& mesh,
                          LTS::Storage& ltsStorage,
                          LTS::Backmap& backmap,
                          time_stepping::TimeManager& timeManager) {
  const auto memkind =
      useUSM() ? seissol::memory::Memkind::DeviceUnifiedMemory : seissol::memory::Memkind::Standard;
  auto sourceClusters = std::vector<seissol::kernels::PointSourceClusterPair>();
  if (sourceType == seissol::initializer::parameters::PointSourceType::NrfSource) {
    logInfo() << "Reading an NRF source (type 42).";
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
    sourceClusters = loadSourceFile<NrfFile>(fileName, mesh, ltsStorage, backmap, memkind);
#else
    logError() << "NRF sources (type 42) need SeisSol to be linked with an (active) Netcdf "
                  "library. However, this is not the case for this build.";
#endif
  } else if (sourceType == seissol::initializer::parameters::PointSourceType::FsrmSource) {
    logInfo() << "Reading an FSRM source (type 50).";
    sourceClusters = loadSourceFile<FsrmFile>(fileName, mesh, ltsStorage, backmap, memkind);
  } else if (sourceType == seissol::initializer::parameters::PointSourceType::None) {
    logInfo() << "No source term specified.";
    sourceClusters =
        std::vector<seissol::kernels::PointSourceClusterPair>(ltsStorage.numChildren());
  } else {
    logError() << "The source type" << static_cast<int>(sourceType)
               << "has been defined, but not yet been implemented in SeisSol.";
  }
  // otherwise, we do not have any source term.

  timeManager.setPointSourcesForClusters(std::move(sourceClusters));
}

} // namespace seissol::sourceterm
