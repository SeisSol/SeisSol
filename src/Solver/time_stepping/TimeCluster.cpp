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
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2017, SeisSol Group
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
 * LTS cluster in SeisSol.
 **/

#include "Parallel/MPI.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SeisSol.h"
#include "TimeCluster.h"
#include <Solver/Interoperability.h>
#include <SourceTerm/PointSource.h>
#include <Kernels/TimeCommon.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Receiver.h>
#include <Monitoring/FlopCounter.hpp>

#include <cassert>
#include <cstring>

//! fortran interoperability
extern seissol::Interoperability e_interoperability;

seissol::time_stepping::TimeCluster::TimeCluster(unsigned int i_clusterId,
                                                 unsigned int i_globalClusterId,
                                                 LayerType layerType,
                                                 double maxTimeStepSize,
                                                 long timeStepRate,
                                                 double timeTolerance,
                                                 bool printProgress,
                                                 struct GlobalData *i_globalData,
                                                 seissol::initializers::Layer *i_clusterData,
                                                 seissol::initializers::Layer *i_dynRupClusterData,
                                                 seissol::initializers::LTS *i_lts,
                                                 seissol::initializers::DynamicRupture *i_dynRup,
                                                 LoopStatistics *i_loopStatistics) :
    AbstractTimeCluster(maxTimeStepSize, timeTolerance, timeStepRate),
    // cluster ids
    m_clusterId(i_clusterId),
    m_globalClusterId(i_globalClusterId),
    layerType(layerType),
    // global data
    m_globalData(i_globalData),
    m_clusterData(i_clusterData),
    m_dynRupClusterData(i_dynRupClusterData),
    m_lts(i_lts),
    m_dynRup(i_dynRup),
    // cells
    m_cellToPointSources(nullptr),
    m_numberOfCellToPointSourcesMappings(0),
    m_pointSources(nullptr),
    m_loopStatistics(i_loopStatistics),
    m_receiverCluster(nullptr),
    printProgress(printProgress)
{
    assert( m_globalData                               != nullptr);
    assert( m_clusterData                              != nullptr);

  // default: no updates are allowed
#ifdef USE_MPI
  //m_sendLtsBuffers                = false;
#endif
  // set timings to zero
  m_receiverTime                  = 0;

  m_dynamicRuptureFaces = i_dynRupClusterData->getNumberOfCells() > 0;

  m_timeKernel.setGlobalData(m_globalData);
  m_localKernel.setGlobalData(m_globalData);
  m_localKernel.setInitConds(&e_interoperability.getInitialConditions());
  m_neighborKernel.setGlobalData(m_globalData);
  m_dynamicRuptureKernel.setGlobalData(m_globalData);

  computeFlops();

  m_regionComputeLocalIntegration = m_loopStatistics->getRegion("computeLocalIntegration");
  m_regionComputeNeighboringIntegration = m_loopStatistics->getRegion("computeNeighboringIntegration");
  m_regionComputeDynamicRupture = m_loopStatistics->getRegion("computeDynamicRupture");
}

seissol::time_stepping::TimeCluster::~TimeCluster() {
#ifndef NDEBUG
  logInfo() << "#(time steps):" << numberOfTimeSteps;
#endif
}

void seissol::time_stepping::TimeCluster::setPointSources( sourceterm::CellToPointSourcesMapping const* i_cellToPointSources,
                                                           unsigned i_numberOfCellToPointSourcesMappings,
                                                           sourceterm::PointSources const* i_pointSources )
{
  m_cellToPointSources = i_cellToPointSources;
  m_numberOfCellToPointSourcesMappings = i_numberOfCellToPointSourcesMappings;
  m_pointSources = i_pointSources;
}

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (m_receiverCluster != nullptr) {
    m_receiverTime = m_receiverCluster->calcReceivers(m_receiverTime, ct.correctionTime, timeStepSize());
  }

}

void seissol::time_stepping::TimeCluster::computeSources() {
  SCOREP_USER_REGION( "computeSources", SCOREP_USER_REGION_TYPE_FUNCTION )

  // Return when point sources not initialised. This might happen if there
  // are no point sources on this rank.
  if (m_numberOfCellToPointSourcesMappings != 0) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned mapping = 0; mapping < m_numberOfCellToPointSourcesMappings; ++mapping) {
      unsigned startSource = m_cellToPointSources[mapping].pointSourcesOffset;
      unsigned endSource =
          m_cellToPointSources[mapping].pointSourcesOffset + m_cellToPointSources[mapping].numberOfPointSources;
      if (m_pointSources->mode == sourceterm::PointSources::NRF) {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceNRF(m_pointSources->mInvJInvPhisAtSources[source],
                                                      m_pointSources->tensor[source],
                                                      m_pointSources->A[source],
                                                      m_pointSources->stiffnessTensor[source],
                                                      m_pointSources->slipRates[source],
                                                      ct.correctionTime,
                                                      ct.correctionTime + timeStepSize(),
                                                      *m_cellToPointSources[mapping].dofs);
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM(m_pointSources->mInvJInvPhisAtSources[source],
                                                       m_pointSources->tensor[source],
                                                       m_pointSources->slipRates[source][0],
                                                       ct.correctionTime,
                                                       ct.correctionTime + timeStepSize(),
                                                       *m_cellToPointSources[mapping].dofs);
        }
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  DRFaceInformation*                    faceInformation                                                   = layerData.var(m_dynRup->faceInformation);
  DRGodunovData*                        godunovData                                                       = layerData.var(m_dynRup->godunovData);
  real**                                timeDerivativePlus                                                = layerData.var(m_dynRup->timeDerivativePlus);
  real**                                timeDerivativeMinus                                               = layerData.var(m_dynRup->timeDerivativeMinus);
  real                                (*imposedStatePlus)[tensor::QInterpolated::size()]                  = layerData.var(m_dynRup->imposedStatePlus);
  real                                (*imposedStateMinus)[tensor::QInterpolated::size()]                 = layerData.var(m_dynRup->imposedStateMinus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                                    = layerData.var(m_dynRup->waveSpeedsPlus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                                   = layerData.var(m_dynRup->waveSpeedsMinus);

  alignas(ALIGNMENT) real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  alignas(ALIGNMENT) real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];

  m_dynamicRuptureKernel.setTimeStepWidth(timeStepSize());
#ifdef _OPENMP
  #pragma omp parallel for schedule(static) private(QInterpolatedPlus,QInterpolatedMinus)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.spaceTimeInterpolation(  faceInformation[face],
                                                    m_globalData,
                                                   &godunovData[face],
                                                    timeDerivativePlus[face],
                                                    timeDerivativeMinus[face],
                                                    QInterpolatedPlus,
                                                    QInterpolatedMinus,
                                                    timeDerivativePlus[prefetchFace],
                                                    timeDerivativeMinus[prefetchFace] );

    e_interoperability.evaluateFrictionLaw( static_cast<int>(faceInformation[face].meshFace),
                                            QInterpolatedPlus,
                                            QInterpolatedMinus,
                                            imposedStatePlus[face],
                                            imposedStateMinus[face],
                                            ct.correctionTime,
                                            m_dynamicRuptureKernel.timePoints,
                                            m_dynamicRuptureKernel.timeWeights,
                                            waveSpeedsPlus[face],
                                            waveSpeedsMinus[face] );
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells());
}


void seissol::time_stepping::TimeCluster::computeDynamicRuptureFlops(seissol::initializers::Layer& layerData)
{
  auto& nonZeroFlops = m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& hardwareFlops = m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];
  nonZeroFlops = 0;
  hardwareFlops = 0;

  DRFaceInformation* faceInformation = layerData.var(m_dynRup->faceInformation);

  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    long long faceNonZeroFlops, faceHardwareFlops;
    m_dynamicRuptureKernel.flopsGodunovState(faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

void seissol::time_stepping::TimeCluster::computeLocalIntegration(seissol::initializers::Layer& i_layerData, bool resetBuffers ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  real l_integrationBuffer[tensor::I::size()] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real* l_bufferPointer;

  real** buffers = i_layerData.var(m_lts->buffers);
  real** derivatives = i_layerData.var(m_lts->derivatives);
  real** displacements = i_layerData.var(m_lts->displacements);
  CellMaterialData* materialData = i_layerData.var(m_lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp{};

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer, tmp) schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++) {
    auto data = loader.entry(l_cell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided = (data.cellInformation.ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers = buffersProvided && ( (data.cellInformation.ltsSetup >> 10) %2 == 0 || resetBuffers ); // they should be reset

    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[l_cell] != nullptr);

      l_bufferPointer = buffers[l_cell];
    } else {
      // work on local buffer
      l_bufferPointer = l_integrationBuffer;
    }

    m_timeKernel.computeAder(timeStepSize(),
                             data,
                             tmp,
                             l_bufferPointer,
                             derivatives[l_cell]);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping (*boundaryMapping)[4] = i_layerData.var(m_lts->boundaryMapping);
    m_localKernel.computeIntegral(l_bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[l_cell],
                                  &boundaryMapping[l_cell],
                                  ct.correctionTime,
                                  timeStepSize()
    );

    // Update displacement
    if (displacements[l_cell] != nullptr) {
      kernel::addVelocity krnl;
      krnl.I = l_bufferPointer;
      krnl.selectVelocity = init::selectVelocity::Values;
      krnl.displacement = displacements[l_cell];
      krnl.execute();
    }

    // TODO: Integrate this step into the kernel
    // We've used a temporary buffer -> need to accumulate update in
    // shared buffer.
    if (!resetMyBuffers && buffersProvided) {
      assert(buffers[l_cell] != nullptr);

      for (unsigned int l_dof = 0; l_dof < tensor::I::size(); ++l_dof) {
        buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
      }
    }
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells());
}

void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData, double subTimeStart ) {
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
  CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
#ifdef USE_PLASTICITY
  updateRelaxTime();
  PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
  real (*pstrain)[7] = i_layerData.var(m_lts->pstrain);
  unsigned numberOTetsWithPlasticYielding = 0;
#endif

  kernels::NeighborData::Loader loader;
  loader.load(*m_lts, i_layerData);

  real *l_timeIntegrated[4];
  real *l_faceNeighbors_prefetch[4];

  //std::cout << "dt_max=" << ct.maxTimeStepSize << " dt=" << timeStepSize() << " t_sub=" << subTimeStart << std::endl;

#ifdef _OPENMP
#ifdef USE_PLASTICITY
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch) reduction(+:numberOTetsWithPlasticYielding)
#else
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch)
#endif
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    auto data = loader.entry(l_cell);
    seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                                   data.cellInformation.ltsSetup,
                                                   data.cellInformation.faceTypes,
                                                   subTimeStart,
                                                   timeStepSize(),
                                                   faceNeighbors[l_cell],
#ifdef _OPENMP
                                                   *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalData->integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
#else
                                                   *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalData->integrationBufferLTS),
#endif
                                                   l_timeIntegrated);

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][1] :
      drMapping[l_cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][2] :
      drMapping[l_cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][3] :
      drMapping[l_cell][3].godunov;

    // fourth face's prefetches
    if (l_cell < (i_layerData.getNumberOfCells()-1) ) {
      l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != FaceType::dynamicRupture) ?
	faceNeighbors[l_cell+1][0] :
	drMapping[l_cell+1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
    }
#endif

    m_neighborKernel.computeNeighborsIntegral( data,
                                               drMapping[l_cell],
#ifdef ENABLE_MATRIX_PREFETCH
                                               l_timeIntegrated, l_faceNeighbors_prefetch
#else
                                               l_timeIntegrated
#endif
                                               );

#ifdef USE_PLASTICITY
  numberOTetsWithPlasticYielding += seissol::kernels::Plasticity::computePlasticity( m_relaxTime,
                                                                                     timeStepSize(),
                                                                                     m_globalData,
                                                                                     &plasticity[l_cell],
                                                                                     data.dofs,
                                                                                     pstrain[l_cell] );
#endif
#ifdef INTEGRATE_QUANTITIES
  seissol::SeisSol::main.postProcessor().integrateQuantities( timeStepSize(),
                                                              i_layerData,
                                                              l_cell,
                                                              dofs[l_cell] );
#endif // INTEGRATE_QUANTITIES
  }

  #ifdef USE_PLASTICITY
  g_SeisSolNonZeroFlopsPlasticity +=
    i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)]
    + numberOTetsWithPlasticYielding * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
  g_SeisSolHardwareFlopsPlasticity +=
    i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)]
    + numberOTetsWithPlasticYielding * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];
  #endif

  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells());
}

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = layerData.var(m_lts->cellInformation);
  for (unsigned cell = 0; cell < layerData.getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    // TODO(Lukas) Maybe include avg. displacement computation here at some point.
    m_timeKernel.flopsAder(cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    m_localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
  }
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLaw)];
  auto& drFlopsHardware = m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLaw)];
  flopsNonZero = 0;
  flopsHardware = 0;
  drFlopsNonZero = 0;
  drFlopsHardware = 0;

  auto* cellInformation = layerData.var(m_lts->cellInformation);
  auto* drMapping = layerData.var(m_lts->drMapping);
  for (unsigned cell = 0; cell < layerData.getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    long long cellDRNonZero, cellDRHardware;
    m_neighborKernel.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
                                            cellInformation[cell].faceRelations,
                                            drMapping[cell],
                                            cellNonZero,
                                            cellHardware,
                                            cellDRNonZero,
                                            cellDRHardware );
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    drFlopsNonZero += cellDRNonZero;
    drFlopsHardware += cellDRHardware;

    /// \todo add lts time integration
    /// \todo add plasticity
  }
}

void seissol::time_stepping::TimeCluster::computeFlops() {
  computeLocalIntegrationFlops(*m_clusterData);
  computeNeighborIntegrationFlops(*m_clusterData);
  computeDynamicRuptureFlops(*m_dynRupClusterData);
  seissol::kernels::Plasticity::flopsPlasticity(
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)],
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)],
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)],
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]
          );
}

namespace seissol::time_stepping {
void TimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct.maxTimeStepSize) {
    lastSubTime = neighborCluster.ct.correctionTime;
  }
}
void TimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster&) {
  // Doesn't do anything
}
void TimeCluster::predict() {
  assert(state == ActorState::Corrected);
  bool resetBuffers = true;
  for (auto& neighbor : neighbors) {
      if (neighbor.ct.timeStepRate > ct.timeStepRate
          && ct.stepsSinceLastSync > neighbor.ct.stepsSinceLastSync) {
          resetBuffers = false;
        }
  }
  if (ct.stepsSinceLastSync == 0) resetBuffers = true;
  if (ct.stepsSinceLastSync == 0) {
    assert(resetBuffers);
  }

  // These methods compute the receivers/sources for both interior and copy cluster
  // and are called in actors for both copy AND interior.
  // TODO(Lukas) Make threadsafe later! (Concurrent copy/interior is unsafe due to shared state)
  writeReceivers();
  computeLocalIntegration(*m_clusterData, resetBuffers);
  computeSources();

  g_SeisSolNonZeroFlopsLocal += m_flops_nonZero[static_cast<int>(ComputePart::Local)];
  g_SeisSolHardwareFlopsLocal += m_flops_hardware[static_cast<int>(ComputePart::Local)];
}
void TimeCluster::correct() {
  assert(state == ActorState::Predicted);

  double subTimeStart = ct.correctionTime - lastSubTime;
  // Note, the following is likely wrong
  // TODO(Lukas) Check scheduling if this is correct
  if (m_dynamicRuptureFaces) {
    computeDynamicRupture(*m_dynRupClusterData);
    g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLaw)];
    g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLaw)];
  }

  computeNeighboringIntegration(*m_clusterData, subTimeStart);

  g_SeisSolNonZeroFlopsNeighbor += m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  g_SeisSolHardwareFlopsNeighbor += m_flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];

  // First cluster calls fault receiver output
  // TODO: Change from iteration based to time based
  // TODO(Lukas): Watch out that we use correct timestepSize here
  // TODO(Lukas) Are we onlyx calling this once?!
  if (m_clusterId == 0) {
    e_interoperability.faultOutput(ct.correctionTime + timeStepSize(), timeStepSize());
  }

  // TODO Use next correction time
  const auto nextCorrectionSteps = ct.nextCorrectionSteps();
  if constexpr (USE_MPI) {
    if (printProgress && ((nextCorrectionSteps % 100) == 0)) {
      const int rank = MPI::mpi.rank();
      logInfo(rank) << "#max-updates since sync: " << nextCorrectionSteps
                    << " @ " << ct.nextCorrectionTime(syncTime);

      }
  }

}

void TimeCluster::reset() {
    AbstractTimeCluster::reset();
    /*
    // TODO(Lukas) Do we need to do this in this way?
    for (unsigned mapping = 0; mapping < m_numberOfCellToPointSourcesMappings; ++mapping) {
        unsigned startSource = m_cellToPointSources[mapping].pointSourcesOffset;
        unsigned endSource =
                m_cellToPointSources[mapping].pointSourcesOffset + m_cellToPointSources[mapping].numberOfPointSources;
        for (unsigned source = startSource; source < endSource; ++source) {
            logInfo(MPI::mpi.rank()) << "Resetting " << source;
            m_pointSources->lastPredictionSteps[source] = -1;
        }
    }
     */
}

void TimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  const auto rank = MPI::mpi.rank();
  logWarning(rank)
      << "No update since " << timeSinceLastUpdate.count()
      << "[s] for global cluster " << m_globalClusterId
      << " with local cluster id " << m_clusterId
      << " at state " << actorStateToString(state)
      << " predTime = " << ct.predictionTime
      << " predictionsSinceSync = " << ct.predictionsSinceLastSync
      << " corrTime = " << ct.correctionTime
      << " correctionsSinceSync = " << ct.stepsSinceLastSync
      << " stepsTillSync = " << ct.stepsUntilSync
      << " mayPredict = " << mayPredict()
      << " mayCorrect = " << mayCorrect()
      << " maySync = " << maySync();
  for (auto& neighbor : neighbors) {
    logWarning(rank)
    << "Neighbor with rate = " << neighbor.ct.timeStepRate
    << "PredTime = " << neighbor.ct.predictionTime
    << "CorrTime = " << neighbor.ct.correctionTime
    << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
    << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }

}

}
