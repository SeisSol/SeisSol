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
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
#include "NRFReader.h"
#include "PointSource.h"
#include "Numerical_aux/Transformation.h"
#include "generated_code/kernel.h"
#include "generated_code/init.h"
#include "generated_code/tensor.h"

#include <Initializer/PointMapper.h>
#include <Solver/Interoperability.h>
#include <utils/logger.h>
#include <cstring>

template<typename T>
class index_sort_by_value
{
private:
    T const* value;
public:
    index_sort_by_value(T const* value) : value(value) {}
    inline bool operator()(unsigned i, unsigned j) const {
        return value[i] < value[j];
    }
};

/**
 * Computes mInvJInvPhisAtSources[i] = |J|^-1 * M_ii^-1 * phi_i(xi, eta, zeta),
 * where xi, eta, zeta is the point in the reference tetrahedron corresponding to x, y, z.
 */
void seissol::sourceterm::computeMInvJInvPhisAtSources(Eigen::Vector3d const& centre,
                                                       real* mInvJInvPhisAtSources,
                                                       unsigned meshId,
                                                       MeshReader const& mesh) {
  auto const& elements = mesh.getElements();
  auto const& vertices = mesh.getVertices();

  double const* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[ elements[meshId].vertices[v] ].coords;
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
  krnl.mInvJInvPhisAtSources = mInvJInvPhisAtSources;
  krnl.JInv = JInv;
  krnl.execute();
}

void seissol::sourceterm::transformNRFSourceToInternalSource( Eigen::Vector3d const&    centre,
                                                              unsigned                  meshId,
                                                              MeshReader const&         mesh,
                                                              Subfault const&           subfault,
                                                              Offsets const&            offsets,
                                                              Offsets const&            nextOffsets,
                                                              double *const             sliprates[3],
                                                              seissol::model::Material* material,
                                                              PointSources&             pointSources,
                                                              unsigned                  index )
{
  computeMInvJInvPhisAtSources(centre, pointSources.mInvJInvPhisAtSources[index], meshId, mesh);

  real* faultBasis = pointSources.tensor[index];
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
  switch(material->getMaterialType()) {
    case seissol::model::MaterialType::anisotropic:
      if (subfault.mu != 0) {
        logError() << "There are specific fault parameters for the fault. This version of SeisSol was compiled for anisotropic materials. This is only compatible if the material around the source is actually isotropic.";
      }
      dynamic_cast<seissol::model::AnisotropicMaterial*>(material)->getFullStiffnessTensor(pointSources.stiffnessTensor[index]);
      break;
    default:
      seissol::model::ElasticMaterial em = *dynamic_cast<seissol::model::ElasticMaterial*>(material);
      em.mu = (subfault.mu == 0.0) ? em.mu : subfault.mu;
      em.getFullStiffnessTensor(pointSources.stiffnessTensor[index]);
      break;
  }
 
  for (unsigned sr = 0; sr < 3; ++sr) {
    unsigned numSamples = nextOffsets[sr] - offsets[sr];
    double const* samples = (numSamples > 0) ? &sliprates[sr][ offsets[sr] ] : NULL;
    samplesToPiecewiseLinearFunction1D( samples,
                                        numSamples,
                                        subfault.tinit,
                                        subfault.timestep,
                                        &pointSources.slipRates[index][sr] );
  }
}

void seissol::sourceterm::Manager::freeSources()
{
  delete[] cmps;
  delete[] sources;
  cmps = NULL;
  sources = NULL;
}

void seissol::sourceterm::Manager::mapPointSourcesToClusters( unsigned const*                 meshIds,
                                                              unsigned                        numberOfSources,
                                                              seissol::initializers::LTSTree* ltsTree,
                                                              seissol::initializers::LTS*     lts,
                                                              seissol::initializers::Lut*     ltsLut )
{
  std::vector<std::vector<unsigned> > clusterToPointSources(ltsTree->numChildren());
  std::vector<std::vector<unsigned> > clusterToMeshIds(ltsTree->numChildren());

  for (unsigned source = 0; source < numberOfSources; ++source) {
    unsigned meshId = meshIds[source];
    unsigned cluster = ltsLut->cluster(meshId);
    clusterToPointSources[cluster].push_back(source);
    clusterToMeshIds[cluster].push_back(meshId);
  }

  cmps = new ClusterMapping[ltsTree->numChildren()];
  for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
    // Determine number of mappings by counting unique mesh Ids
    std::sort(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
    std::vector<unsigned>::iterator last = std::unique(clusterToMeshIds[cluster].begin(), clusterToMeshIds[cluster].end());
    unsigned numberOfMappings = 0;
    for (std::vector<unsigned>::iterator it = clusterToMeshIds[cluster].begin(); it != last; ++it) {
      unsigned meshId = *it;
      for (unsigned dup = 0; dup < seissol::initializers::Lut::MaxDuplicates && ltsLut->ltsId(lts->dofs.mask, meshId, dup) != std::numeric_limits<unsigned>::max(); ++dup) {
        ++numberOfMappings;
      } 
    }
    
    cmps[cluster].sources           = new unsigned[ clusterToPointSources[cluster].size() ];
    cmps[cluster].numberOfSources   = clusterToPointSources[cluster].size();
    cmps[cluster].cellToSources     = new CellToPointSourcesMapping[ numberOfMappings ];
    cmps[cluster].numberOfMappings  = numberOfMappings;
    
    for (unsigned source = 0; source < clusterToPointSources[cluster].size(); ++source) {
      cmps[cluster].sources[source] = clusterToPointSources[cluster][source];
    }
    std::sort(cmps[cluster].sources, cmps[cluster].sources + cmps[cluster].numberOfSources, index_sort_by_value<unsigned>(meshIds));
    
    unsigned clusterSource = 0;
    unsigned mapping = 0;
    while (clusterSource < cmps[cluster].numberOfSources) {
      unsigned meshId = meshIds[ cmps[cluster].sources[clusterSource] ];
      unsigned next = clusterSource + 1;
      while (next < cmps[cluster].numberOfSources && meshIds[ cmps[cluster].sources[next] ] == meshId) {
        ++next;
      }
      
      for (unsigned ltsId, dup = 0; dup < seissol::initializers::Lut::MaxDuplicates && (ltsId = ltsLut->ltsId(lts->dofs.mask, meshId, dup)) != std::numeric_limits<unsigned>::max(); ++dup) {
        cmps[cluster].cellToSources[mapping].dofs = &ltsTree->var(lts->dofs)[ltsId];
        cmps[cluster].cellToSources[mapping].pointSourcesOffset = clusterSource;
        cmps[cluster].cellToSources[mapping].numberOfPointSources = next - clusterSource;
        ++mapping;
      }      
      
      clusterSource = next;
    }
    assert(mapping == cmps[cluster].numberOfMappings);
  }
}

void seissol::sourceterm::Manager::loadSourcesFromFSRM( double const*                   momentTensor,
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
                                                        time_stepping::TimeManager&     timeManager )
{
  freeSources();

  int rank = seissol::MPI::mpi.rank();
  
  logInfo(rank) << "<--------------------------------------------------------->";
  logInfo(rank) << "<                      Point sources                      >";
  logInfo(rank) << "<--------------------------------------------------------->";

  short* contained = new short[numberOfSources];
  unsigned* meshIds = new unsigned[numberOfSources];
  Eigen::Vector3d* centres3 = new Eigen::Vector3d[numberOfSources];

  for (int source = 0; source < numberOfSources; ++source) {
    centres3[source](0) = centres[3*source];
    centres3[source](1) = centres[3*source + 1];
    centres3[source](2) = centres[3*source + 2];
  }

  logInfo(rank) << "Finding meshIds for point sources...";

  initializers::findMeshIds(centres3, mesh, numberOfSources, contained, meshIds);

#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring point sources for MPI...";
  initializers::cleanDoubles(contained, numberOfSources);
#endif

  unsigned* originalIndex = new unsigned[numberOfSources];
  unsigned numSources = 0;
  for (int source = 0; source < numberOfSources; ++source) {
    originalIndex[numSources] = source;
    meshIds[numSources] = meshIds[source];
    numSources += contained[source];
  }
  delete[] contained;

  logInfo(rank) << "Mapping point sources to LTS cells...";
  mapPointSourcesToClusters(meshIds, numSources, ltsTree, lts, ltsLut);

  real localMomentTensor[3][3];
  for (unsigned i = 0; i < 9; ++i) {
    *(&localMomentTensor[0][0] + i) = momentTensor[i];
  }
  real localVelocityComponent[3];
  for (unsigned i = 0; i < 3; i++) {
    localVelocityComponent[i] = velocityComponent[i];
  }
  
  sources = new PointSources[ltsTree->numChildren()];
  for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
    sources[cluster].mode                  = PointSources::FSRM;
    sources[cluster].numberOfSources       = cmps[cluster].numberOfSources;
    sources[cluster].originalIndex.reserve(sources[cluster].numberOfSources);
    int error = posix_memalign(reinterpret_cast<void**>(&sources[cluster].mInvJInvPhisAtSources), ALIGNMENT, cmps[cluster].numberOfSources*tensor::mInvJInvPhisAtSources::size()*sizeof(real));
    if (error) {
      logError() << "posix_memalign failed in source term manager.";
    }
    error = posix_memalign(reinterpret_cast<void**>(&sources[cluster].tensor), ALIGNMENT, cmps[cluster].numberOfSources*PointSources::TensorSize*sizeof(real));
    if (error) {
      logError() << "posix_memalign failed in source term manager.";
    }
    sources[cluster].slipRates.resize(cmps[cluster].numberOfSources);

    for (unsigned clusterSource = 0; clusterSource < cmps[cluster].numberOfSources; ++clusterSource) {
      unsigned sourceIndex = cmps[cluster].sources[clusterSource];
      unsigned fsrmIndex = originalIndex[sourceIndex];
      sources[cluster].originalIndex[clusterSource] = fsrmIndex;
      computeMInvJInvPhisAtSources(centres3[fsrmIndex],
              sources[cluster].mInvJInvPhisAtSources[clusterSource],
              meshIds[sourceIndex], mesh);

      transformMomentTensor( localMomentTensor,
                             localVelocityComponent,
                             strikes[fsrmIndex],
                             dips[fsrmIndex],
                             rakes[fsrmIndex],
                             sources[cluster].tensor[clusterSource]);

      for (unsigned i = 0; i < 9; ++i) {
        sources[cluster].tensor[clusterSource][i] *= areas[fsrmIndex];
      }
      seissol::model::Material& material = ltsLut->lookup(lts->material, meshIds[sourceIndex] - 1).local;
      for (unsigned i = 0; i < 3; ++i) {
        sources[cluster].tensor[clusterSource][6+i] /= material.rho;
      }

      samplesToPiecewiseLinearFunction1D( &timeHistories[fsrmIndex * numberOfSamples],
                                          numberOfSamples,
                                          onsets[fsrmIndex],
                                          timestep,
                                          &sources[cluster].slipRates[clusterSource][0] );
    }
  }
  delete[] originalIndex;
  delete[] meshIds;
  delete[] centres3;

  timeManager.setPointSourcesForClusters(cmps, sources);
  
  logInfo(rank) << ".. finished point source initialization.";
}

// TODO Add support for passive netCDF
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
void seissol::sourceterm::Manager::loadSourcesFromNRF(  char const*                     fileName,
                                                        MeshReader const&               mesh,
                                                        seissol::initializers::LTSTree* ltsTree,
                                                        seissol::initializers::LTS*     lts,
                                                        seissol::initializers::Lut*     ltsLut,
                                                        time_stepping::TimeManager&     timeManager )
{
  freeSources();

  int rank = seissol::MPI::mpi.rank();

  logInfo(rank) << "<--------------------------------------------------------->";
  logInfo(rank) << "<                      Point sources                      >";
  logInfo(rank) << "<--------------------------------------------------------->";

  logInfo(rank) << "Reading" << fileName;
  NRF nrf;
  readNRF(fileName, nrf);

  short* contained = new short[nrf.source];
  unsigned* meshIds = new unsigned[nrf.source];

  logInfo(rank) << "Finding meshIds for point sources...";
  initializers::findMeshIds(nrf.centres, mesh, nrf.source, contained, meshIds);

#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring point sources for MPI...";
  initializers::cleanDoubles(contained, nrf.source);
#endif

  unsigned* originalIndex = new unsigned[nrf.source];
  unsigned numSources = 0;
  for (unsigned source = 0; source < nrf.source; ++source) {
    originalIndex[numSources] = source;
    meshIds[numSources] = meshIds[source];
    numSources += contained[source];
  }
  delete[] contained;

  logInfo(rank) << "Mapping point sources to LTS cells...";
  mapPointSourcesToClusters(meshIds, numSources, ltsTree, lts, ltsLut);
  
  sources = new PointSources[ltsTree->numChildren()];
  for (unsigned cluster = 0; cluster < ltsTree->numChildren(); ++cluster) {
    sources[cluster].mode                  = PointSources::NRF;
    sources[cluster].numberOfSources       = cmps[cluster].numberOfSources;
    sources[cluster].originalIndex.reserve(sources[cluster].numberOfSources);
    int error = posix_memalign(reinterpret_cast<void**>(&sources[cluster].mInvJInvPhisAtSources), ALIGNMENT, cmps[cluster].numberOfSources*tensor::mInvJInvPhisAtSources::size()*sizeof(real));
    if (error) {
      logError() << "posix_memalign failed in source term manager.";
    }
    error = posix_memalign(reinterpret_cast<void**>(&sources[cluster].tensor), ALIGNMENT, cmps[cluster].numberOfSources*PointSources::TensorSize*sizeof(real));
    if (error) {
      logError() << "posix_memalign failed in source term manager.";
    }
    sources[cluster].A.resize(cmps[cluster].numberOfSources);
    sources[cluster].stiffnessTensor.resize(cmps[cluster].numberOfSources);
    sources[cluster].slipRates.resize(cmps[cluster].numberOfSources);

    for (unsigned clusterSource = 0; clusterSource < cmps[cluster].numberOfSources; ++clusterSource) {
      unsigned sourceIndex = cmps[cluster].sources[clusterSource];
      unsigned nrfIndex = originalIndex[sourceIndex];
      sources[cluster].originalIndex[clusterSource] =  nrfIndex;
      transformNRFSourceToInternalSource( nrf.centres[nrfIndex],
                                          meshIds[sourceIndex],
                                          mesh,
                                          nrf.subfaults[nrfIndex],
                                          nrf.sroffsets[nrfIndex],
                                          nrf.sroffsets[nrfIndex+1],
                                          nrf.sliprates,
                                          &ltsLut->lookup(lts->material, meshIds[sourceIndex]).local,
                                          sources[cluster],
                                          clusterSource );
    }
  }
  delete[] originalIndex;
  delete[] meshIds;

  timeManager.setPointSourcesForClusters(cmps, sources);
  
  logInfo(rank) << ".. finished point source initialization.";
}
#endif // defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
