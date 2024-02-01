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

#include <Equations/Setup.h>

#include <Kernels/TimeCommon.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Receiver.h>
#include <Monitoring/FlopCounter.hpp>
#include <Monitoring/instrumentation.fpp>

#include <Numerical_aux/Quadrature.h>

#include <cassert>
#include <cstring>

#include <generated_code/kernel.h>

//! fortran interoperability
extern seissol::Interoperability e_interoperability;

void setStarMatrix( real* i_AT,
                    real* i_BT,
                    real* i_CT,
                    real  i_grad[3],
                    real* o_starMatrix );

seissol::time_stepping::TimeCluster::TimeCluster(unsigned int i_clusterId, unsigned int i_globalClusterId,
                                                 bool usePlasticity,
                                                 LayerType layerType, double maxTimeStepSize,
                                                 long timeStepRate, bool printProgress,
                                                 DynamicRuptureScheduler *dynamicRuptureScheduler,
                                                 CompoundGlobalData i_globalData,
                                                 seissol::initializers::Layer *i_clusterData,
                                                 seissol::initializers::Layer *dynRupInteriorData,
                                                 seissol::initializers::Layer *dynRupCopyData,
                                                 seissol::initializers::LTS *i_lts,
                                                 seissol::initializers::DynamicRupture* i_dynRup,
                                                 seissol::dr::friction_law::FrictionSolver* i_FrictionSolver,
                                                 dr::output::OutputManager* i_faultOutputManager,
                                                 LoopStatistics *i_loopStatistics,
                                                 ActorStateStatistics* actorStateStatistics) :
    AbstractTimeCluster(maxTimeStepSize, timeStepRate),
    // cluster ids
    usePlasticity(usePlasticity),
    m_globalDataOnHost( i_globalData.onHost ),
    m_globalDataOnDevice(i_globalData.onDevice ),
    m_clusterData(i_clusterData),
    // global data
    dynRupInteriorData(dynRupInteriorData),
    dynRupCopyData(dynRupCopyData),
    m_lts(i_lts),
    m_dynRup(i_dynRup),
    frictionSolver(i_FrictionSolver),
    faultOutputManager(i_faultOutputManager),
    m_cellToPointSources(nullptr),
    m_numberOfCellToPointSourcesMappings(0),
    m_pointSources(nullptr),
    // cells
    m_loopStatistics(i_loopStatistics),
    actorStateStatistics(actorStateStatistics),
    m_receiverCluster(nullptr),
    layerType(layerType),
    printProgress(printProgress),
    m_clusterId(i_clusterId),
    m_globalClusterId(i_globalClusterId),
    dynamicRuptureScheduler(dynamicRuptureScheduler)
{
    // assert all pointers are valid
    assert( m_clusterData                              != nullptr );
    assert( m_globalDataOnHost                         != nullptr );
    if constexpr (seissol::isDeviceOn()) {
        assert( m_globalDataOnDevice                   != nullptr );
    }

  // set timings to zero
  m_receiverTime                  = 0;

  m_timeKernel.setGlobalData(i_globalData);
  m_localKernel.setGlobalData(i_globalData);
  m_localKernel.setInitConds(&e_interoperability.getInitialConditions());
  m_neighborKernel.setGlobalData(i_globalData);
  m_dynamicRuptureKernel.setGlobalData(i_globalData);

  m_krnlNonlVolPrototype.kDivM = i_globalData.onHost->stiffnessMatrices;

  m_nonlinearInterpolation.V3mTo2n = i_globalData.onHost->faceToNodalMatrices;
  m_nonlSurfIntPrototype.V3mTo2nTWDivM = i_globalData.onHost->nodalFluxMatrices;

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
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
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
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
}

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  DRFaceInformation* faceInformation = layerData.var(m_dynRup->faceInformation);
  DRGodunovData* godunovData = layerData.var(m_dynRup->godunovData);
  DREnergyOutput* drEnergyOutput = layerData.var(m_dynRup->drEnergyOutput);
  real** timeDerivativePlus = layerData.var(m_dynRup->timeDerivativePlus);
  real** timeDerivativeMinus = layerData.var(m_dynRup->timeDerivativeMinus);
  auto* qInterpolatedPlus = layerData.var(m_dynRup->qInterpolatedPlus);
  auto* qInterpolatedMinus = layerData.var(m_dynRup->qInterpolatedMinus);

  m_dynamicRuptureKernel.setTimeStepWidth(timeStepSize());
  frictionSolver->computeDeltaT(m_dynamicRuptureKernel.timePoints);

#pragma omp parallel
  {
  LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation");
  }
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.spaceTimeInterpolation(faceInformation[face],
                                                  m_globalDataOnHost,
                                                  &godunovData[face],
                                                  &drEnergyOutput[face],
                                                  timeDerivativePlus[face],
                                                  timeDerivativeMinus[face],
                                                  qInterpolatedPlus[face],
                                                  qInterpolatedMinus[face],
                                                  timeDerivativePlus[prefetchFace],
                                                  timeDerivativeMinus[prefetchFace]);
  }
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureSpaceTimeInterpolation");
  LIKWID_MARKER_START("computeDynamicRuptureFrictionLaw");
  }

  SCOREP_USER_REGION_BEGIN(myRegionHandle, "computeDynamicRuptureFrictionLaw", SCOREP_USER_REGION_TYPE_COMMON )
  frictionSolver->evaluate(layerData,
                           m_dynRup,
                           ct.correctionTime,
                           m_dynamicRuptureKernel.timeWeights);
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_globalClusterId);
}
#else

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  if (layerData.getNumberOfCells() > 0) {
    // compute space time interpolation part

    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    m_dynamicRuptureKernel.setTimeStepWidth(timeStepSize());
    frictionSolver->computeDeltaT(m_dynamicRuptureKernel.timePoints);

    auto &table = layerData.getConditionalTable<inner_keys::Dr>();
    m_dynamicRuptureKernel.batchedSpaceTimeInterpolation(table);
    device.api->popLastProfilingMark();

    device.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    frictionSolver->evaluate(layerData,
                             m_dynRup,
                             ct.correctionTime,
                             m_dynamicRuptureKernel.timeWeights);
    device.api->popLastProfilingMark();
  }
  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_globalClusterId);
}
#endif


void seissol::time_stepping::TimeCluster::computeDynamicRuptureFlops( seissol::initializers::Layer& layerData,
                                                                      long long&                    nonZeroFlops,
                                                                      long long&                    hardwareFlops )
{
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

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration(seissol::initializers::Layer& i_layerData, bool resetBuffers ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  real l_integrationBuffer[tensor::I::size()] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real* l_bufferPointer;

  real** buffers = i_layerData.var(m_lts->buffers);
  real** derivatives = i_layerData.var(m_lts->derivatives);
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
                             derivatives[l_cell],
                             ct.correctionTime,
                             true);

    #ifdef USE_DAMAGEDELASTIC

    // real epsInitxx = 4.63e-4; // eps_xx0
    // real epsInityy = -1.85e-3; // eps_yy0
    // real epsInitzz = 4.63e-4; // eps_zz0
    // real epsInitxy = 1.11e-3; // eps_xx0
    // real epsInityz = -0e-1; // eps_yy0
    // real epsInitzx = -0e-1; // eps_zz0

    // real epsInitxx = -9.26e-4; // eps_xx0
    // real epsInityy = -9.26e-4; // eps_yy0
    // real epsInitzz = -9.26e-4; // eps_zz0
    // real epsInitxy = 1.11e-3; // eps_xx0
    // real epsInityz = -0e-1; // eps_yy0
    // real epsInitzx = -0e-1; // eps_zz0

    // tpv 5
    // real epsInitxx = 3.73854e-4; // eps_xx0
    // real epsInityy = -1.4963e-3; // eps_yy0
    // real epsInitzz = 3.73854e-4; // eps_zz0
    // real epsInitxy = 1.0909e-3; // eps_xx0
    // real epsInityz = -0e-1; // eps_yy0
    // real epsInitzx = -0e-1; // eps_zz0

    real epsInitxx = -1.8738e-4; // eps_xx0
    real epsInityy = -1.1225e-3; // eps_yy0
    real epsInitzz = -1.8738e-4; // eps_zz0
    real epsInitxy = 1.0909e-3; // eps_xy0
    real epsInityz = -0e-1; // eps_yz0
    real epsInitzx = -0e-1; // eps_zx0
    real const damage_para1 = data.material.local.Cd; // 1.2e-4*2;
    real const break_coeff = 1e1*damage_para1;
    real const beta_alpha = 0.03;

    real aB0 = 7.43e9;
    real aB1 = -22.14e9;
    real aB2 = 20.93e9;
    real aB3 = -6.067e9;

    // std::cout << data.material.local.Cd << std::endl;
    // real const damage_para2 = 3e-6;
    // real const lambda0 = 9.71e10; // data.material.local.lambda0
    // real const mu0 = 8.27e10; // data.material.local.mu0
    // Compute the Q at quadrature points in space and time
    /// Get quadrature points in time
    double timePoints[CONVERGENCE_ORDER];
    double timeWeights[CONVERGENCE_ORDER];
    seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);
    for (unsigned int point = 0; point < CONVERGENCE_ORDER; ++point) {
      timePoints[point] = 0.5 * (timeStepSize() * timePoints[point] + timeStepSize());
      timeWeights[point] = 0.5 * timeStepSize() * timeWeights[point];
    }

    /// Get Q_{lp}(tau_z) at different time quadrature points
    alignas(PAGESIZE_STACK) real QInterpolatedBody[CONVERGENCE_ORDER][tensor::Q::size()]; // initializations?
    alignas(PAGESIZE_STACK) real* QInterpolatedBodyi;
    alignas(PAGESIZE_STACK) real QInterpolatedBodyNodal[CONVERGENCE_ORDER][tensor::QNodal::size()];
    alignas(PAGESIZE_STACK) real* QInterpolatedBodyNodali;

    kernel::damageConvertToNodal d_converToKrnl;
    for (unsigned int timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval){
      QInterpolatedBodyi = QInterpolatedBody[timeInterval];
      QInterpolatedBodyNodali = QInterpolatedBodyNodal[timeInterval];
      m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, derivatives[l_cell], QInterpolatedBodyi);
      /// Convert Q_{lp}(tau_z) in modal basis to QN_{ip}(tau_z) in nodal basis
      d_converToKrnl.v = init::v::Values;
      d_converToKrnl.QNodal = QInterpolatedBodyNodali;
      d_converToKrnl.Q = QInterpolatedBodyi;
      d_converToKrnl.execute();
    }

    alignas(PAGESIZE_STACK) real FInterpolatedBody[CONVERGENCE_ORDER][tensor::QNodal::size()] = {{0}}; // initializations?

    alignas(PAGESIZE_STACK) real FluxInterpolatedBodyX[CONVERGENCE_ORDER][tensor::QNodal::size()] = {{0}};
    alignas(PAGESIZE_STACK) real FluxInterpolatedBodyY[CONVERGENCE_ORDER][tensor::QNodal::size()] = {{0}};
    alignas(PAGESIZE_STACK) real FluxInterpolatedBodyZ[CONVERGENCE_ORDER][tensor::QNodal::size()] = {{0}};

    alignas(PAGESIZE_STACK) real sxxNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
    alignas(PAGESIZE_STACK) real syyNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
    alignas(PAGESIZE_STACK) real szzNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
    alignas(PAGESIZE_STACK) real sxyNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
    alignas(PAGESIZE_STACK) real syzNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
    alignas(PAGESIZE_STACK) real szxNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};

    for (unsigned int timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval){
      real* exxNodal = (QInterpolatedBodyNodal[timeInterval] + 0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* eyyNodal = (QInterpolatedBodyNodal[timeInterval] + 1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* ezzNodal = (QInterpolatedBodyNodal[timeInterval] + 2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* alphaNodal = (QInterpolatedBodyNodal[timeInterval] + 9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* breakNodal = (QInterpolatedBodyNodal[timeInterval] + 10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      // std::cout << exxNodal[0] << " " << solNData[0] << std::endl;

      real* exyNodal = (QInterpolatedBodyNodal[timeInterval] + 3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* eyzNodal = (QInterpolatedBodyNodal[timeInterval] + 4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* ezxNodal = (QInterpolatedBodyNodal[timeInterval] + 5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* vxNodal = (QInterpolatedBodyNodal[timeInterval] + 6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* vyNodal = (QInterpolatedBodyNodal[timeInterval] + 7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
      real* vzNodal = (QInterpolatedBodyNodal[timeInterval] + 8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);

      // real alpha_ave = 0.0;
      // real break_ave = 0.0;
      // real w_ave = 1.0/NUMBER_OF_ALIGNED_BASIS_FUNCTIONS;
      // for (unsigned int q = 0; q<NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q){
      //   break_ave += breakNodal[q] * w_ave;
      //   alpha_ave += alphaNodal[q] * w_ave;
      // }

      real alpha_ave = alphaNodal[0];
      real break_ave = breakNodal[0];
      for (unsigned int q = 0; q<NUMBER_OF_ALIGNED_BASIS_FUNCTIONS-1; ++q){
        break_ave = std::max(break_ave, breakNodal[q]);
        alpha_ave = std::max(alpha_ave, alphaNodal[q]);
      }

      for (unsigned int q = 0; q<NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q){
        real EspI = (exxNodal[q]+epsInitxx) + (eyyNodal[q]+epsInityy) + (ezzNodal[q]+epsInitzz);
        real EspII = (exxNodal[q]+epsInitxx)*(exxNodal[q]+epsInitxx)
          + (eyyNodal[q]+epsInityy)*(eyyNodal[q]+epsInityy)
          + (ezzNodal[q]+epsInitzz)*(ezzNodal[q]+epsInitzz)
          + 2*(exyNodal[q]+epsInitxy)*(exyNodal[q]+epsInitxy)
          + 2*(eyzNodal[q]+epsInityz)*(eyzNodal[q]+epsInityz)
          + 2*(ezxNodal[q]+epsInitzx)*(ezxNodal[q]+epsInitzx);
        real xi;
        if (EspII > 1e-30){
          xi = EspI / std::sqrt(EspII);
        } else{
          xi = 0.0;
        }

        // Compute alpha_{cr}
        real aCR = (3.0*xi*xi - 3.0)*data.material.local.gammaR*data.material.local.gammaR
        + 6.0*xi*data.material.local.gammaR*data.material.local.xi0*data.material.local.gammaR
        + 4.0*data.material.local.xi0*data.material.local.gammaR*data.material.local.xi0*data.material.local.gammaR;

        real bCR = - (8.0*data.material.local.mu0 + 6.0*data.material.local.lambda0) * data.material.local.xi0*data.material.local.gammaR
        - xi * (xi*xi* data.material.local.lambda0 + 6.0*data.material.local.mu0) * data.material.local.gammaR;

        real cCR = 4.0 * data.material.local.mu0 * data.material.local.mu0
        + 6.0 * data.material.local.mu0 * data.material.local.lambda0;

        real alphaCR1q = ( -bCR - std::sqrt(bCR*bCR - 4.0*aCR*cCR) )/(2.0*aCR);
        real alphaCR2q = 2.0*data.material.local.mu0
          /data.material.local.gammaR/(xi+2.0*data.material.local.xi0);

         real alphaCRq = 1.0;
        if (alphaCR1q > 0.0){
          if (alphaCR2q > 0.0){
            alphaCRq = std::min(1.0,
              std::min( alphaCR1q, alphaCR2q )
            );
          }
        }

        if (xi + data.material.local.xi0 > 0) {
          if (alpha_ave < 0.9 ){
            if (break_ave < 0.85 ){
              FInterpolatedBody[timeInterval][10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
                (1 - breakNodal[q]) * 1.0/(std::exp( (alphaCRq - alphaNodal[q])/beta_alpha ) + 1.0) * break_coeff
                  *data.material.local.gammaR * EspII * (xi + data.material.local.xi0);
              FInterpolatedBody[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
                (1 - breakNodal[q]) * damage_para1
                  *data.material.local.gammaR * EspII * (xi + data.material.local.xi0);
            } else {
              FInterpolatedBody[timeInterval][10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
              FInterpolatedBody[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
            }
          } else {
            FInterpolatedBody[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
            FInterpolatedBody[timeInterval][10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
          }
        } else if (alpha_ave > 5e-1 ) {
          FInterpolatedBody[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
            0.0*damage_para1
              *data.material.local.gammaR * EspII * (xi + data.material.local.xi0);
          FInterpolatedBody[timeInterval][10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
            0.0*damage_para1
              *data.material.local.gammaR * EspII * (xi + data.material.local.xi0);
        }
        else {
          FInterpolatedBody[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
          FInterpolatedBody[timeInterval][10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        }

        // 1.0e0/damage_para2*(damage_para1*(exxNodal[q] + eyyNodal[q] + ezzNodal[q])*(exxNodal[q] + eyyNodal[q] + ezzNodal[q]) - alphaNodal[q]);

        // Compute nonlinear flux term

        // damage stress
        real mu_eff = data.material.local.mu0 - alphaNodal[q]*data.material.local.gammaR*data.material.local.xi0
            - 0.5*alphaNodal[q]*data.material.local.gammaR*xi;
        real sxx_s = data.material.local.lambda0*EspI
                      - alphaNodal[q]*data.material.local.gammaR * std::sqrt(EspII)
                      + 2*mu_eff*(exxNodal[q]+epsInitxx);
        real syy_s = data.material.local.lambda0*EspI
                      - alphaNodal[q]*data.material.local.gammaR * std::sqrt(EspII)
                      + 2*mu_eff*(eyyNodal[q]+epsInityy);

        real szz_s = data.material.local.lambda0*EspI
                      - alphaNodal[q]*data.material.local.gammaR * std::sqrt(EspII)
                      + 2*mu_eff*(ezzNodal[q]+epsInitzz);

        real sxy_s = 2*mu_eff*(exyNodal[q]+epsInitxy);
        real syz_s = 2*mu_eff*(eyzNodal[q]+epsInityz);
        real szx_s = 2*mu_eff*(ezxNodal[q]+epsInitzx);

        // breakage stress
        real sxx_b = (2.0*aB2 + 3.0*xi*aB3)*EspI
                      + aB1 * std::sqrt(EspII)
                      + (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(exxNodal[q]+epsInitxx);
        real syy_b = (2.0*aB2 + 3.0*xi*aB3)*EspI
                      + aB1 * std::sqrt(EspII)
                      + (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(eyyNodal[q]+epsInityy);
        real szz_b = (2.0*aB2 + 3.0*xi*aB3)*EspI
                      + aB1 * std::sqrt(EspII)
                      + (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(ezzNodal[q]+epsInitzz);

        real sxy_b = (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(exyNodal[q]+epsInitxy);
        real syz_b = (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(eyzNodal[q]+epsInityz);
        real szx_b = (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(ezxNodal[q]+epsInitzx);

        sxxNodal[q] = (1-breakNodal[q])*sxx_s + breakNodal[q]*sxx_b;
        syyNodal[q] = (1-breakNodal[q])*syy_s + breakNodal[q]*syy_b;
        szzNodal[q] = (1-breakNodal[q])*szz_s + breakNodal[q]*szz_b;
        sxyNodal[q] = (1-breakNodal[q])*sxy_s + breakNodal[q]*sxy_b;
        syzNodal[q] = (1-breakNodal[q])*syz_s + breakNodal[q]*syz_b;
        szxNodal[q] = (1-breakNodal[q])*szx_s + breakNodal[q]*szx_b;

        // //--- x-dir
        FluxInterpolatedBodyX[timeInterval][0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -vxNodal[q];
        FluxInterpolatedBodyX[timeInterval][1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyX[timeInterval][2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyX[timeInterval][3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -0.5*vyNodal[q];
        FluxInterpolatedBodyX[timeInterval][4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyX[timeInterval][5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -0.5*vzNodal[q];
        FluxInterpolatedBodyX[timeInterval][6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -sxxNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyX[timeInterval][7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -sxyNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyX[timeInterval][8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -szxNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyX[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        //--- y-dir
        FluxInterpolatedBodyY[timeInterval][0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyY[timeInterval][1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -vyNodal[q];
        FluxInterpolatedBodyY[timeInterval][2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyY[timeInterval][3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -0.5*vxNodal[q];
        FluxInterpolatedBodyY[timeInterval][4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -0.5*vzNodal[q];
        FluxInterpolatedBodyY[timeInterval][5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyY[timeInterval][6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -sxyNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyY[timeInterval][7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -syyNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyY[timeInterval][8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -syzNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyY[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        //--- z-dir
        FluxInterpolatedBodyZ[timeInterval][0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyZ[timeInterval][1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyZ[timeInterval][2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -vzNodal[q];
        FluxInterpolatedBodyZ[timeInterval][3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FluxInterpolatedBodyZ[timeInterval][4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -0.5*vyNodal[q];
        FluxInterpolatedBodyZ[timeInterval][5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -0.5*vxNodal[q];
        FluxInterpolatedBodyZ[timeInterval][6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -szxNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyZ[timeInterval][7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -syzNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyZ[timeInterval][8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -szzNodal[q]/data.material.local.rho;
        FluxInterpolatedBodyZ[timeInterval][9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      }
    }

    /// Convert Q_{lp} at the initial time step from modal to nodal space
    //// Quadrature in time and nodal space and project back to modal space
    kernel::damageIntegration d_timeIntegration;
    for (unsigned int timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval){
      d_timeIntegration.vInv = init::vInv::Values;
      d_timeIntegration.QTNodal(timeInterval) = FInterpolatedBody[timeInterval];
      d_timeIntegration.Q = data.dofs;
      d_timeIntegration.Tweight = (timeWeights[timeInterval]);
      d_timeIntegration.execute(timeInterval);
    }

    // for (unsigned int timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval){
    //   real* exxNodal = (QInterpolatedBodyNodal[timeInterval] + 0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* eyyNodal = (QInterpolatedBodyNodal[timeInterval] + 1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* ezzNodal = (QInterpolatedBodyNodal[timeInterval] + 2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* alphaNodal = (QInterpolatedBodyNodal[timeInterval] + 9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* breakNodal = (QInterpolatedBodyNodal[timeInterval] + 10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   // std::cout << exxNodal[0] << " " << solNData[0] << std::endl;

    //   real* exyNodal = (QInterpolatedBodyNodal[timeInterval] + 3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* eyzNodal = (QInterpolatedBodyNodal[timeInterval] + 4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* ezxNodal = (QInterpolatedBodyNodal[timeInterval] + 5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* vxNodal = (QInterpolatedBodyNodal[timeInterval] + 6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* vyNodal = (QInterpolatedBodyNodal[timeInterval] + 7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   real* vzNodal = (QInterpolatedBodyNodal[timeInterval] + 8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    //   for (unsigned int q = 0; q<NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q){
    //     real EspI = (exxNodal[q]+epsInitxx) + (eyyNodal[q]+epsInityy) + (ezzNodal[q]+epsInitzz);
    //     real EspII = (exxNodal[q]+epsInitxx)*(exxNodal[q]+epsInitxx)
    //       + (eyyNodal[q]+epsInityy)*(eyyNodal[q]+epsInityy)
    //       + (ezzNodal[q]+epsInitzz)*(ezzNodal[q]+epsInitzz)
    //       + 2*(exyNodal[q]+epsInitxy)*(exyNodal[q]+epsInitxy)
    //       + 2*(eyzNodal[q]+epsInityz)*(eyzNodal[q]+epsInityz)
    //       + 2*(ezxNodal[q]+epsInitzx)*(ezxNodal[q]+epsInitzx);
    //     real xi;
    //     if (EspII > 1e-30){
    //       xi = EspI / std::sqrt(EspII);
    //     } else{
    //       xi = 0.0;
    //     }

    //     // Compute alpha_{cr}
    //     real aCR = (3.0*xi*xi - 3.0)*data.material.local.gammaR*data.material.local.gammaR
    //     + 6.0*xi*data.material.local.gammaR*data.material.local.xi0*data.material.local.gammaR
    //     + 4.0*data.material.local.xi0*data.material.local.gammaR*data.material.local.xi0*data.material.local.gammaR;

    //     real bCR = - (8.0*data.material.local.mu0 + 6.0*data.material.local.lambda0) * data.material.local.xi0*data.material.local.gammaR
    //     - xi * (xi*xi* data.material.local.lambda0 + 6.0*data.material.local.mu0) * data.material.local.gammaR;

    //     real cCR = 4.0 * data.material.local.mu0 * data.material.local.mu0
    //     + 6.0 * data.material.local.mu0 * data.material.local.lambda0;

    //     real alphaCR1q = ( -bCR - std::sqrt(bCR*bCR - 4.0*aCR*cCR) )/(2.0*aCR);
    //     real alphaCR2q = 2.0*data.material.local.mu0
    //       /data.material.local.gammaR/(xi+2.0*data.material.local.xi0);

    //     real alphaCRq = std::min(1.0,
    //       std::min(alphaCR1q,alphaCR2q)
    //     );

    //     if (data.dofs[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] > 0.9) {

    //       std::cout << (1 - breakNodal[q]) * 1.0/(std::exp( (alphaCRq - alphaNodal[q])/beta_alpha ) + 1.0) * break_coeff << std::endl;

    //       std::cout << q << " qs: " << data.dofs[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] << std::endl;

    //       std::cout << q << " q2s: " << breakNodal[q] << std::endl;
    //     }

    //   }
    // }

    // Do integration of the nonlinear volumetric flux here
    /// Integrate V^-1_li * F_ip^d * Theta_ed * K^e_kl in time (*TweightN)
    kernel::nonlinearVolumeIntegration d_nonlinearVolumeIntegration = m_krnlNonlVolPrototype;
    for (unsigned int timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval){
      /// Convert F^d_{lp}(tau_z) in nodal basis back to modal basis
      d_nonlinearVolumeIntegration.vInv = init::vInv::Values;
      d_nonlinearVolumeIntegration.gradXiEtaZetaX0 = data.localIntegration.gradXiEtaZeta[0][0]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaX1 = data.localIntegration.gradXiEtaZeta[0][1]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaX2 = data.localIntegration.gradXiEtaZeta[0][2]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaY0 = data.localIntegration.gradXiEtaZeta[1][0]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaY1 = data.localIntegration.gradXiEtaZeta[1][1]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaY2 = data.localIntegration.gradXiEtaZeta[1][2]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaZ0 = data.localIntegration.gradXiEtaZeta[2][0]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaZ1 = data.localIntegration.gradXiEtaZeta[2][1]*timeWeights[timeInterval];
      d_nonlinearVolumeIntegration.gradXiEtaZetaZ2 = data.localIntegration.gradXiEtaZeta[2][2]*timeWeights[timeInterval];

      d_nonlinearVolumeIntegration.Q = data.dofs;
      d_nonlinearVolumeIntegration.FluxVolX(timeInterval) = FluxInterpolatedBodyX[timeInterval];
      d_nonlinearVolumeIntegration.FluxVolY(timeInterval) = FluxInterpolatedBodyY[timeInterval];
      d_nonlinearVolumeIntegration.FluxVolZ(timeInterval) = FluxInterpolatedBodyZ[timeInterval];
      // d_nonlinearVolumeIntegration.TweightN = (timeWeights[timeInterval]);
      // std::cout << timeWeights[timeInterval] << std::endl;
      d_nonlinearVolumeIntegration.execute(timeInterval);
    }
    //================================END OF DAMAGE===============================
    #endif

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

    for (unsigned face = 0; face < 4; ++face) {
      auto& curFaceDisplacements = data.faceDisplacements[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr
          && data.cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = m_globalDataOnHost->V3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
        addVelocityKrnl.faceDisplacement = data.faceDisplacements[face];
        addVelocityKrnl.I = l_bufferPointer;
        addVelocityKrnl.execute(face);
      }
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

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration(
  seissol::initializers::Layer& i_layerData,
  bool resetBuffers) {

  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  std::cout<<'Using ACL_DEVICE==================\n';

  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  auto& dataTable = i_layerData.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = i_layerData.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = i_layerData.getConditionalTable<inner_keys::Indices>();

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp;

  m_timeKernel.computeBatchedAder(timeStepSize(),
                                  tmp,
                                  dataTable,
                                  materialTable,
                                  ct.correctionTime,
                                  true);

  m_localKernel.computeBatchedIntegral(dataTable,
                                       materialTable,
                                       indicesTable,
                                       loader,
                                       tmp,
                                       ct.correctionTime,
                                       timeStepSize());

  auto defaultStream = device.api->getDefaultStream();

  for (unsigned face = 0; face < 4; ++face) {
    ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
    if (dataTable.find(key) != dataTable.end()) {
      auto &entry = dataTable[key];
      // NOTE: integrated velocities have been computed implicitly, i.e
      // it is 6th, 7the and 8th columns of integrated dofs

      kernel::gpu_addVelocity displacementKrnl;
      displacementKrnl.faceDisplacement = entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
      displacementKrnl.integratedVelocities = const_cast<real const**>(entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
      displacementKrnl.V3mTo2nFace = m_globalDataOnDevice->V3mTo2nFace;

      // Note: this kernel doesn't require tmp. memory
      displacementKrnl.numElements = entry.get(inner_keys::Wp::Id::FaceDisplacement)->getSize();
      displacementKrnl.streamPtr = device.api->getDefaultStream();
      displacementKrnl.execute(face);
    }
  }

  ConditionalKey key = ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
  if (dataTable.find(key) != dataTable.end()) {
    auto &entry = dataTable[key];

    if (resetBuffers) {
      device.algorithms.streamBatchedData((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
                                          (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                                          tensor::I::Size,
                                          (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                                          defaultStream);
    }
    else {
      device.algorithms.accumulateBatchedData((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
                                              (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                                              tensor::I::Size,
                                              (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                                              defaultStream);
    }
  }

  device.api->synchDevice();
  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);

  device.api->popLastProfilingMark();
}
#endif // ACL_DEVICE

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegration(seissol::initializers::Layer& i_layerData,
                                                                        double subTimeStart) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(i_layerData, subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(i_layerData, subTimeStart);
  }
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData,
                                                                         double subTimeStart) {
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  auto& table = i_layerData.getConditionalTable<inner_keys::Wp>();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(m_timeKernel,
                                                        subTimeStart,
                                                        timeStepSize(),
                                                        table);
  m_neighborKernel.computeBatchedNeighborsIntegral(table);

  if (usePlasticity) {
    PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
    unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(m_oneMinusIntegratingFactor,
                                                                                      timeStepSize(),
                                                                                      m_tv,
                                                                                      m_globalDataOnDevice,
                                                                                      table,
                                                                                      plasticity);

    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsPlasticity(
        i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsPlasticity(
        i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]);
  }

  device.api->synchDevice();
  device.api->popLastProfilingMark();
  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);
}
#endif // ACL_DEVICE

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = layerData.var(m_lts->cellInformation);
  for (unsigned cell = 0; cell < layerData.getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    m_timeKernel.flopsAder(cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    m_localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    // Contribution from displacement/integrated displacement
    for (unsigned face = 0; face < 4; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::freeSurfaceGravity) {
        const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
        GravitationalFreeSurfaceBc::getFlopsDisplacementFace(face,
                                                             cellInformation[cell].faceTypes[face]);
        flopsNonZero += nonZeroFlopsDisplacement;
        flopsHardware += hardwareFlopsDisplacement;
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];
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
  computeDynamicRuptureFlops(*dynRupInteriorData,
                             m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)],
                             m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
  computeDynamicRuptureFlops(*dynRupCopyData,
                             m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)],
                             m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
  seissol::kernels::Plasticity::flopsPlasticity(
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)],
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)],
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)],
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]
          );
}

void seissol::time_stepping::TimeCluster::updateDerivatives() {
  // SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  // Access the neighboring solutions
  seissol::initializers::Layer& i_layerData = *m_clusterData;
  real** derivatives = i_layerData.var(m_lts->derivatives);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  // kernels::LocalTmp tmp{};

#ifdef _OPENMP
  #pragma omp parallel
    {
#endif

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++) {
      auto data = loader.entry(l_cell);
      if (derivatives[l_cell] != NULL) {
        for (unsigned dof = 0; dof < tensor::Q::size(); ++dof) {
          // zero time integration buffers
          derivatives[l_cell][dof] = data.dofs[dof];
          // std::cout << derivatives[l_cell][dof] << std::endl;
        }
      }
    }

#ifdef _OPENMP
    }
#endif

  // m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);
}

void seissol::time_stepping::TimeCluster::updateMaterialLocal(seissol::initializers::Layer& i_layerData) {
  // SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  // m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // / To recompute the matrices, information required includes:
  // mesh information, material input,

  MeshReader& meshReader = seissol::SeisSol::main.meshReader();
  std::vector<Element> const& elements = meshReader.getElements();
  std::vector<Vertex> const& vertices = meshReader.getVertices();

  // CellMaterialData* materialData = i_layerData.var(m_lts->material);
  CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);

  // Access the neighboring solutions
  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  // real** derivatives = i_layerData.var(m_lts->buffers);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  // kernels::LocalTmp tmp{};

#ifdef _OPENMP
  #pragma omp parallel
    {
#endif
    CellMaterialData* materialData = i_layerData.var(m_lts->material);
    real ATData[tensor::star::size(0)];
    real ATtildeData[tensor::star::size(0)];
    real BTData[tensor::star::size(1)];
    real CTData[tensor::star::size(2)];
    auto AT = init::star::view<0>::create(ATData);
    // AT with elastic parameters in local coordinate system, used for flux kernel
    auto ATtilde = init::star::view<0>::create(ATtildeData);
    auto BT = init::star::view<0>::create(BTData);
    auto CT = init::star::view<0>::create(CTData);

    real TData[seissol::tensor::T::size()];
    real TinvData[seissol::tensor::Tinv::size()];
    auto T = init::T::view::create(TData);
    auto Tinv = init::Tinv::view::create(TinvData);

    real QgodLocalData[tensor::QgodLocal::size()];
    real QgodNeighborData[tensor::QgodNeighbor::size()];
    auto QgodLocal = init::QgodLocal::view::create(QgodLocalData);
    auto QgodNeighbor = init::QgodNeighbor::view::create(QgodNeighborData);
    kernel::cellAve m_cellAverageKernel;
    // real Q_aveData[NUMBER_OF_QUANTITIES];
    real Q_aveData[tensor::QAve::size()];
    auto Q_ave = init::QAve::view::create(Q_aveData);

    // real Q_aveNeighborData[tensor::QAve::size()];
    // auto Q_aveNeighbor = init::QAve::view::create(Q_aveNeighborData);

    // seissol::model::DamagedElasticMaterial changed_materialLocal;
    // seissol::model::DamagedElasticMaterial changed_materialNeighbor[4];

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++) {
      auto data = loader.entry(l_cell);

      real * derivatives_neighbor[4];
      for (unsigned int i_nei=0; i_nei<4; i_nei++){
        derivatives_neighbor[i_nei] = faceNeighbors[l_cell][i_nei];
      }

      // TODO: compute cell-average of solutions and compute new material properties
      // real Q_ave[NUMBER_OF_QUANTITIES];
      m_cellAverageKernel.phiAve = init::phiAve::Values;
      m_cellAverageKernel.Q = data.dofs;
      m_cellAverageKernel.QAve = Q_aveData;
      m_cellAverageKernel.execute();

      // std::cout << data.dofs[0+0] << std::endl;
      // real epsInitxx = -0e-2; // eps_xx0
      // real epsInityy = -0e-1; // eps_yy0
      // real epsInitzz = -0e-1; // eps_zz0
      // real epsInitxy = -0e-1; // eps_xx0
      // real epsInityz = -0e-1; // eps_yy0
      // real epsInitzx = -0e-1; // eps_zz0

      // real epsInitxx = 4.63e-4; // eps_xx0
      // real epsInityy = -1.85e-3; // eps_yy0
      // real epsInitzz = 4.63e-4; // eps_zz0
      // real epsInitxy = 1.11e-3; // eps_xx0
      // real epsInityz = -0e-1; // eps_yy0
      // real epsInitzx = -0e-1; // eps_zz0

      // real epsInitxx = -9.26e-4; // eps_xx0
      // real epsInityy = -9.26e-4; // eps_yy0
      // real epsInitzz = -9.26e-4; // eps_zz0
      // real epsInitxy = 1.11e-3; // eps_xx0
      // real epsInityz = -0e-1; // eps_yy0
      // real epsInitzx = -0e-1; // eps_zz0

      // tpv 5
      // real epsInitxx = 3.73854e-4; // eps_xx0
      // real epsInityy = -1.4963e-3; // eps_yy0
      // real epsInitzz = 3.73854e-4; // eps_zz0
      // real epsInitxy = 1.0909e-3; // eps_xx0
      // real epsInityz = -0e-1; // eps_yy0
      // real epsInitzx = -0e-1; // eps_zz0

      real epsInitxx = -1.8738e-4; // eps_xx0
      real epsInityy = -1.1225e-3; // eps_yy0
      real epsInitzz = -1.8738e-4; // eps_zz0
      real epsInitxy = 1.0909e-3; // eps_xy0
      real epsInityz = -0e-1; // eps_yz0
      real epsInitzx = -0e-1; // eps_zx0

      // END TODO
      real EspI = (Q_aveData[0]+epsInitxx) + (Q_aveData[1]+epsInityy) + (Q_aveData[2]+epsInitzz);
      real EspII = (Q_aveData[0]+epsInitxx)*(Q_aveData[0]+epsInitxx)
        +  (Q_aveData[1]+epsInityy)*(Q_aveData[1]+epsInityy)
        +  (Q_aveData[2]+epsInitzz)*(Q_aveData[2]+epsInitzz)
        +  2*(Q_aveData[3]+epsInitxy)*(Q_aveData[3]+epsInitxy)
        +  2*(Q_aveData[4]+epsInityz)*(Q_aveData[4]+epsInityz)
        +  2*(Q_aveData[5]+epsInitzx)*(Q_aveData[5]+epsInitzx);

      real xi;
      // real xiInv;
      if (EspII > 1e-30){
        xi = EspI / std::sqrt(EspII);
      } else{
        xi = 0.0;
      }

      // if ( std::abs(xi) > 1e-1){
      //   xiInv = 1 / xi;
      // } else{
      //   xiInv = 0.0;
      // }

      real alphaAve = Q_aveData[9];
      real breakAve = Q_aveData[10];

      real lambda0 = materialData[l_cell].local.lambda0;
      real mu0 = materialData[l_cell].local.mu0;
      real beta_m = 0e2;

      // real aB0 = 4.95e9;
      // real aB1 = -18.89e9;
      // real aB2 = 23.96e9;
      // real aB3 = -10.112e9;

      real aB0 = 7.43e9;
      real aB1 = -22.14e9;
      real aB2 = 20.93e9;
      real aB3 = -6.067e9;

      unsigned int meshId = data.localIntegration.globalMeshId;

      // changed_materialLocal = materialData[l_cell].local;
      // changed_materialLocal.mu
      // = materialData[l_cell].local.mu*(1-1e6*EspI);


      // for (unsigned side = 0; side < 4; ++side){
      //   changed_materialNeighbor[side].lambda
      //   = materialData[l_cell].neighbor[side].lambda*0.9;
      // }

       materialData[l_cell].local.mu = (1-breakAve) * (mu0
          - alphaAve*materialData[l_cell].local.xi0*materialData[l_cell].local.gammaR
          - 0.5*alphaAve*materialData[l_cell].local.gammaR*xi)
          + breakAve * (
            (aB0 + 0.5*aB1*xi - 0.5*aB3*xi*xi*xi)
          );
       materialData[l_cell].local.lambda = (1-breakAve) * (lambda0
        - alphaAve*materialData[l_cell].local.gammaR*(Q_aveData[0]+epsInitxx)/std::sqrt(EspII))
        + breakAve * (
          (2.0*aB2 + 3.0*aB3*xi) + aB1*(Q_aveData[0]+epsInitxx)/std::sqrt(EspII)
        );
       materialData[l_cell].local.gamma = alphaAve*materialData[l_cell].local.gammaR;

       materialData[l_cell].local.epsxx_alpha = (Q_aveData[0]+epsInitxx);
       materialData[l_cell].local.epsyy_alpha = (Q_aveData[1]+epsInityy);
       materialData[l_cell].local.epszz_alpha = (Q_aveData[2]+epsInitzz);
       materialData[l_cell].local.epsxy_alpha = (Q_aveData[3]+epsInitxy);
       materialData[l_cell].local.epsyz_alpha = (Q_aveData[4]+epsInityz);
       materialData[l_cell].local.epszx_alpha = (Q_aveData[5]+epsInitzx);

      //  for (unsigned side = 0; side < 4; ++side){
      //   materialData[l_cell].neighbor[side].mu = 8.27e10*(1-1e7*EspI);
      // }

      // if (meshId > 1000 && meshId < 1002){
      //   std::cout << seissol::MPI::mpi.rank() << "=========================" << std::endl;
      //   // std::cout << materialData[l_cell].local.lambda << std::endl;
      //   // std::cout << changed_material.lambda << std::endl;
      //   // std::cout << changed_materialLocal.mu/materialData[l_cell].local.mu << std::endl;
      //   std::cout << materialData[l_cell].local.sigmaxx_alpha << std::endl;
      //   // std::cout << Q_ave.size() << std::endl;
      // }

      /// global coordinates of the vertices
      real x[4];
      real y[4];
      real z[4];
      real gradXi[3];
      real gradEta[3];
      real gradZeta[3];

      // Iterate over all 4 vertices of the tetrahedron
      for (unsigned vertex = 0; vertex < 4; ++vertex) {
        VrtxCoords const& coords = vertices[ elements[meshId].vertices[vertex] ].coords;
        x[vertex] = coords[0];
        y[vertex] = coords[1];
        z[vertex] = coords[2];
      }

      seissol::transformations::tetrahedronGlobalToReferenceJacobian( x, y, z, gradXi, gradEta, gradZeta );
      // seissol::model::getTransposedCoefficientMatrix( changed_materialLocal, 0, AT );
      // seissol::model::getTransposedCoefficientMatrix( changed_materialLocal, 1, BT );
      // seissol::model::getTransposedCoefficientMatrix( changed_materialLocal, 2, CT );
      seissol::model::getTransposedCoefficientMatrix( materialData[l_cell].local, 0, AT );
      seissol::model::getTransposedCoefficientMatrix( materialData[l_cell].local, 1, BT );
      seissol::model::getTransposedCoefficientMatrix( materialData[l_cell].local, 2, CT );
      setStarMatrix(ATData, BTData, CTData, gradXi, data.localIntegration.starMatrices[0]);
      setStarMatrix(ATData, BTData, CTData, gradEta, data.localIntegration.starMatrices[1]);
      setStarMatrix(ATData, BTData, CTData, gradZeta, data.localIntegration.starMatrices[2]);


      // if (meshId > 1000 && meshId < 1002){
      //   std::cout << seissol::MPI::mpi.rank() << "=========================" << std::endl;
      //   std::cout << CTData[26]*gradXi[2] << std::endl;
      // }

      double volume = MeshTools::volume(elements[meshId], vertices);

      auto ATtildeBC = init::star::view<0>::create(data.localIntegration.ATtildeBC);
      seissol::model::getTransposedCoefficientMatrix( materialData[l_cell].local, 0, ATtildeBC );

      // for (int i_out = 0; i_out<50; ++i_out){
      //   std::cout
      //             // << faceNeighbors[l_cell][0][20*i_out+0] << " "
      //             << derivatives[l_cell][20*i_out+0]
      //             << "s";
      // }
      // std::cout << faceNeighbors[l_cell][0][50] << " "<< std::endl;

      for (unsigned side = 0; side < 4; ++side) {
        // seissol::model::getTransposedGodunovState(  materialData[l_cell].local,
        //                                             materialData[l_cell].neighbor[side],
        //                                             cellInformation[l_cell].faceTypes[side],
        //                                             QgodLocal,
        //                                             QgodNeighbor );

        // // const real *deri_n_s = derivatives_neighbor[side];
        // // derivatives_neighbor[side] = faceNeighbors[l_cell][side];
        // if (cellInformation[l_cell].faceTypes[side] != FaceType::outflow &&
        // cellInformation[l_cell].faceTypes[side] != FaceType::dynamicRupture ) {
        //   // auto solutions = init::Q::view::create(derivatives_neighbor[side]);
        //   std::cout << derivatives_neighbor[side][0] << ' '  << std::endl;
        // }

        lambda0 = materialData[l_cell].neighbor[side].lambda0;
        mu0 = materialData[l_cell].neighbor[side].mu0;

        if (cellInformation[l_cell].faceTypes[side] != FaceType::outflow &&
        cellInformation[l_cell].faceTypes[side] != FaceType::dynamicRupture ) {
          m_cellAverageKernel.phiAve = init::phiAve::Values;
          m_cellAverageKernel.Q = derivatives_neighbor[side];
          m_cellAverageKernel.QAve = Q_aveData;
          m_cellAverageKernel.execute();

          // real EspINeigh = (Q_aveData[0] + Q_aveData[1] + Q_aveData[2])/timeStepSize();
          // real alphaAveNeigh = Q_aveData[9]/timeStepSize();

          real EspINeigh = ((Q_aveData[0]+epsInitxx) + (Q_aveData[1]+epsInityy) + (Q_aveData[2]+epsInitzz));
          real EspII = (Q_aveData[0]+epsInitxx)*(Q_aveData[0]+epsInitxx)
            +  (Q_aveData[1]+epsInityy)*(Q_aveData[1]+epsInityy)
            +  (Q_aveData[2]+epsInitzz)*(Q_aveData[2]+epsInitzz)
            +  2*(Q_aveData[3]+epsInitxy)*(Q_aveData[3]+epsInitxy)
            +  2*(Q_aveData[4]+epsInityz)*(Q_aveData[4]+epsInityz)
            +  2*(Q_aveData[5]+epsInitzx)*(Q_aveData[5]+epsInitzx);
          real alphaAveNeigh = Q_aveData[9];
          real breakAveNeigh = Q_aveData[10];
          real xi;
          real xiInv;
          if (EspII > 1e-30){
            xi = EspINeigh / std::sqrt(EspII);
          } else{
            xi = 0.0;
          }

          // if ( std::abs(xi) > 1e-1){
          //   xiInv = 1 / xi;
          // } else{
          //   xiInv = 0.0;
          // }

          // std::cout << xi << ", " << xiInv  << ", " << alphaAve << std::endl;
          // std::cout << materialData[l_cell].local.xi0 << ", "
          //   << materialData[l_cell].neighbor[side].gammaR  << ", "
          //   << materialData[l_cell].neighbor[side].mu << std::endl;
          materialData[l_cell].neighbor[side].mu = (1-breakAveNeigh) * (mu0
            - alphaAveNeigh*materialData[l_cell].neighbor[side].xi0*materialData[l_cell].neighbor[side].gammaR
            - 0.5*alphaAveNeigh*materialData[l_cell].neighbor[side].gammaR*xi)
            + breakAveNeigh * (
                (aB0 + 0.5*aB1*xi - 0.5*aB3*xi*xi*xi)
              );
          materialData[l_cell].neighbor[side].lambda = (1-breakAveNeigh) * (lambda0
            - alphaAveNeigh*materialData[l_cell].neighbor[side].gammaR*(Q_aveData[0]+epsInitxx)/std::sqrt(EspII))
            + breakAveNeigh * (
              (2.0*aB2 + 3.0*aB3*xi) + aB1*(Q_aveData[0]+epsInitxx)/std::sqrt(EspII)
            );
          materialData[l_cell].neighbor[side].gamma = alphaAveNeigh*materialData[l_cell].neighbor[side].gammaR;

          materialData[l_cell].neighbor[side].epsxx_alpha = (Q_aveData[0]+epsInitxx);
          materialData[l_cell].neighbor[side].epsyy_alpha = (Q_aveData[1]+epsInityy);
          materialData[l_cell].neighbor[side].epszz_alpha = (Q_aveData[2]+epsInitzz);
          materialData[l_cell].neighbor[side].epsxy_alpha = (Q_aveData[3]+epsInitxy);
          materialData[l_cell].neighbor[side].epsyz_alpha = (Q_aveData[4]+epsInityz);
          materialData[l_cell].neighbor[side].epszx_alpha = (Q_aveData[5]+epsInitzx);
        }


        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);
        double surface = MeshTools::surface(normal);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);

        real NLocalData[6*6];
        seissol::model::getBondMatrix(normal, tangent1, tangent2, NLocalData);
        if (materialData[l_cell].local.getMaterialType() == seissol::model::MaterialType::anisotropic) {
          seissol::model::getTransposedGodunovState(  seissol::model::getRotatedMaterialCoefficients(NLocalData, *dynamic_cast<seissol::model::AnisotropicMaterial*>(&materialData[l_cell].local)),
                                                      seissol::model::getRotatedMaterialCoefficients(NLocalData, *dynamic_cast<seissol::model::AnisotropicMaterial*>(&materialData[l_cell].neighbor[side])),
                                                      cellInformation[l_cell].faceTypes[side],
                                                      QgodLocal,
                                                      QgodNeighbor );
          seissol::model::getTransposedCoefficientMatrix( seissol::model::getRotatedMaterialCoefficients(NLocalData, *dynamic_cast<seissol::model::AnisotropicMaterial*>(&materialData[l_cell].local)), 0, ATtilde );
        } else {
          seissol::model::getTransposedGodunovState(  materialData[l_cell].local,
                                                      materialData[l_cell].neighbor[side],
                                                      cellInformation[l_cell].faceTypes[side],
                                                      QgodLocal,
                                                      QgodNeighbor );
          seissol::model::getTransposedCoefficientMatrix( materialData[l_cell].local, 0, ATtilde );
        }

        // if (meshId > 1000 && meshId < 1002){
        //   std::cout << seissol::MPI::mpi.rank() << "=========================" << std::endl;
        //   std::cout << materialData[l_cell].local.lambda << std::endl;
        //   std::cout << &materialData[l_cell].local.lambda << std::endl;
        // }



        // Calculate transposed T instead
        seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, T, Tinv);

        // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
        // must be subtracted.
        real fluxScale = -2.0 * surface / (6.0 * volume);

        kernel::computeFluxSolverLocal localKrnl;
        localKrnl.fluxScale = fluxScale;
        localKrnl.AplusT = data.localIntegration.nApNm1[side];
        localKrnl.QgodLocal = QgodLocalData;
        localKrnl.T = TData;
        localKrnl.Tinv = TinvData;
        localKrnl.star(0) = ATtildeData;
        localKrnl.execute();

        kernel::computeFluxSolverNeighbor neighKrnl;
        neighKrnl.fluxScale = fluxScale;
        neighKrnl.AminusT = data.neighboringIntegration.nAmNm1[side];
        neighKrnl.QgodNeighbor = QgodNeighborData;
        neighKrnl.T = TData;
        neighKrnl.Tinv = TinvData;
        neighKrnl.star(0) = ATtildeData;
        if (cellInformation[l_cell].faceTypes[side] == FaceType::dirichlet ||
            cellInformation[l_cell].faceTypes[side] == FaceType::freeSurfaceGravity) {
          // Already rotated!
          neighKrnl.Tinv = init::identityT::Values;
        }
        neighKrnl.execute();


        // auto ATtildeBC = init::star::view<0>::create(data.localIntegration.ATtildeBC);
        // seissol::model::getTransposedCoefficientMatrix( materialData[l_cell].local, 0, ATtildeBC );
      }


    }

#ifdef _OPENMP
    }
#endif

  // m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);
}

namespace seissol::time_stepping {
ActResult TimeCluster::act() {
  actorStateStatistics->enter(state);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics->enter(state);
  return result;
}

void TimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct.maxTimeStepSize) {
    lastSubTime = neighborCluster.ct.correctionTime;
  }
  #ifdef USE_DAMAGEDELASTIC
  else{
    lastSubTime = neighborCluster.ct.correctionTime;
  }
  #endif
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
  if (ct.stepsSinceLastSync == 0) {
    resetBuffers = true;
  }

  // 02.02.2023, change material matrices based on the solution at the last time-step
  // real aa;

  // kernels::LocalData::Loader loader;
  // loader.load(*m_lts, *m_clusterData);
  // auto data = loader.entry(10);
  // for (unsigned int ii = 26; ii<27; ++ii){
  //   std::cout << data.localIntegration.starMatrices[0][ii] << " ";
  //   aa = data.localIntegration.starMatrices[0][ii];
  // }
  // std::cout << std::endl;

  // std::cout << "called" << std::endl;

  // updateDerivates(*m_clusterData);

  updateMaterialLocal(*m_clusterData);

  // kernels::LocalData::Loader loader1;
  // loader1.load(*m_lts, *m_clusterData);
  // auto data1 = loader1.entry(10);
  // for (unsigned int ii = 26; ii<27; ++ii){
  //   std::cout << "after: " << data1.localIntegration.starMatrices[0][ii] - aa;
  // }
  // std::cout << std::endl;


  // These methods compute the receivers/sources for both interior and copy cluster
  // and are called in actors for both copy AND interior.
  writeReceivers();

  computeLocalIntegration(*m_clusterData, resetBuffers);
  computeSources();
  // updateMaterialLocal(*m_clusterData);

  seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsLocal(m_flops_nonZero[static_cast<int>(ComputePart::Local)]);
  seissol::SeisSol::main.flopCounter().incrementHardwareFlopsLocal(m_flops_hardware[static_cast<int>(ComputePart::Local)]);
}
void TimeCluster::correct() {
  assert(state == ActorState::Predicted);
  // std::cout << "layer type: " << layerType << std::endl;

  /* Sub start time of width respect to the next cluster; use 0 if not relevant, for example in GTS.
   * LTS requires to evaluate a partial time integration of the derivatives. The point zero in time refers to the derivation of the surrounding time derivatives, which
   * coincides with the last completed time step of the next cluster. The start/end of the time step is the start/end of this clusters time step relative to the zero point.
   *   Example:
   *                                              5 dt
   *   |-----------------------------------------------------------------------------------------| <<< Time stepping of the next cluster (Cn) (5x larger than the current).
   *   |                 |                 |                 |                 |                 |
   *   |*****************|*****************|+++++++++++++++++|                 |                 | <<< Status of the current cluster.
   *   |                 |                 |                 |                 |                 |
   *   |-----------------|-----------------|-----------------|-----------------|-----------------| <<< Time stepping of the current cluster (Cc).
   *   0                 dt               2dt               3dt               4dt               5dt
   *
   *   In the example above two clusters are illustrated: Cc and Cn. Cc is the current cluster under consideration and Cn the next cluster with respect to LTS terminology.
   *   Cn is currently at time 0 and provided Cc with derivatives valid until 5dt. Cc updated already twice and did its last full update to reach 2dt (== subTimeStart). Next
   *   computeNeighboringCopy is called to accomplish the next full update to reach 3dt (+++). Besides working on the buffers of own buffers and those of previous clusters,
   *   Cc needs to evaluate the time prediction of Cn in the interval [2dt, 3dt].
   */
  double subTimeStart = ct.correctionTime - lastSubTime;

  // Note, if this is a copy layer actor, we need the FL_Copy and the FL_Int.
  // Otherwise, this is an interior layer actor, and we need only the FL_Int.
  // We need to avoid computing it twice.
  if (dynamicRuptureScheduler->hasDynamicRuptureFaces()) {
    if (dynamicRuptureScheduler->mayComputeInterior(ct.stepsSinceStart)) {
      computeDynamicRupture(*dynRupInteriorData);
      seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      dynamicRuptureScheduler->setLastCorrectionStepsInterior(ct.stepsSinceStart);
    }
    if (layerType == Copy) {
      computeDynamicRupture(*dynRupCopyData);
      seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      dynamicRuptureScheduler->setLastCorrectionStepsCopy((ct.stepsSinceStart));
    }

  }
  computeNeighboringIntegration(*m_clusterData, subTimeStart);

  seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsNeighbor(m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)]);
  seissol::SeisSol::main.flopCounter().incrementHardwareFlopsNeighbor(m_flops_hardware[static_cast<int>(ComputePart::Neighbor)]);
  seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)]);
  seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)]);

  // First cluster calls fault receiver output
  // Call fault output only if both interior and copy parts of DR were computed
  // TODO: Change from iteration based to time based
  if (dynamicRuptureScheduler->isFirstClusterWithDynamicRuptureFaces()
      && dynamicRuptureScheduler->mayComputeFaultOutput(ct.stepsSinceStart)) {
    faultOutputManager->writePickpointOutput(ct.correctionTime + timeStepSize(), timeStepSize());
    dynamicRuptureScheduler->setLastFaultOutput(ct.stepsSinceStart);
  }

  // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
  const auto nextCorrectionSteps = ct.nextCorrectionSteps();
  if constexpr (USE_MPI) {
    if (printProgress && (((nextCorrectionSteps / timeStepRate) % 100) == 0)) {
      const int rank = MPI::mpi.rank();
      logInfo(rank) << "#max-updates since sync: " << nextCorrectionSteps
                    << " @ " << ct.nextCorrectionTime(syncTime);

      }
  }

}

void TimeCluster::reset() {
    AbstractTimeCluster::reset();
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

unsigned int TimeCluster::getClusterId() const {
  return m_clusterId;
}

unsigned int TimeCluster::getGlobalClusterId() const {
  return m_globalClusterId;
}

LayerType TimeCluster::getLayerType() const {
  return layerType;
}
void TimeCluster::setReceiverTime(double receiverTime) {
  m_receiverTime = receiverTime;
}

}

