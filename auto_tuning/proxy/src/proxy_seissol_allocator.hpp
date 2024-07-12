/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * Copyright (c) 2013-2014, SeisSol Group
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
 **/

#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/LTS.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/GlobalData.h"
#include "Solver/time_stepping/MiniSeisSol.cpp"
#include <yateto.h>
#include <unordered_set>
#include <Parallel/Runtime/Stream.hpp>

#ifdef ACL_DEVICE
#include <device.h>
#include <unordered_set>
#include "Initializer/BatchRecorders/Recorders.h"
#include "Solver/Pipeline/DrPipeline.h"
#endif

seissol::initializer::LTSTree               *m_ltsTree{nullptr};
seissol::initializer::LTS                   m_lts;
seissol::initializer::LTSTree               *m_dynRupTree{nullptr};
seissol::initializer::DynamicRupture        m_dynRup;

GlobalData m_globalDataOnHost;
GlobalData m_globalDataOnDevice;

real* m_fakeDerivatives = nullptr;
real* m_fakeDerivativesHost = nullptr;

seissol::kernels::Time      m_timeKernel;
seissol::kernels::Local     m_localKernel;
seissol::kernels::Neighbor  m_neighborKernel;
seissol::kernels::DynamicRupture m_dynRupKernel;

seissol::memory::ManagedAllocator *m_allocator{nullptr};

seissol::parallel::runtime::StreamRuntime* runtime;

namespace tensor = seissol::tensor;
seissol::SeisSol* seissolInstance = nullptr;

void initGlobalData() {
  seissol::initializer::GlobalDataInitializerOnHost::init(m_globalDataOnHost,
                                                           *m_allocator,
                                                           seissol::memory::Standard,
                                                           seissolInstance);

  CompoundGlobalData globalData{};
  globalData.onHost = &m_globalDataOnHost;
  globalData.onDevice = nullptr;
  if constexpr (seissol::isDeviceOn()) {
    seissol::initializer::GlobalDataInitializerOnDevice::init(m_globalDataOnDevice,
                                                               *m_allocator,
                                                               seissol::memory::DeviceGlobalMemory,
                                                               seissolInstance);
    globalData.onDevice = &m_globalDataOnDevice;
  }
  m_timeKernel.setGlobalData(globalData);
  m_localKernel.setGlobalData(globalData);
  m_neighborKernel.setGlobalData(globalData);
  m_dynRupKernel.setGlobalData(globalData);
}

unsigned int initDataStructures(unsigned int i_cells, bool enableDynamicRupture) {
  // init RNG
  srand48(i_cells);
  m_lts.addTo(*m_ltsTree, false); // proxy does not use plasticity
  m_ltsTree->setNumberOfTimeClusters(1);
  m_ltsTree->fixate();
  
  seissol::initializer::TimeCluster& cluster = m_ltsTree->child(0);
  cluster.child<Ghost>().setNumberOfCells(0);
  cluster.child<Copy>().setNumberOfCells(0);
  cluster.child<Interior>().setNumberOfCells(i_cells);
  
  seissol::initializer::Layer& layer = cluster.child<Interior>();
  layer.setBucketSize(m_lts.buffersDerivatives, sizeof(real) * tensor::I::size() * layer.getNumberOfCells());
  
  m_ltsTree->allocateVariables();
  m_ltsTree->touchVariables();
  m_ltsTree->allocateBuckets();
  
  if (enableDynamicRupture) {
    m_dynRup.addTo(*m_dynRupTree);
    m_dynRupTree->setNumberOfTimeClusters(1);
    m_dynRupTree->fixate();
    
    seissol::initializer::TimeCluster& cluster = m_dynRupTree->child(0);
    cluster.child<Ghost>().setNumberOfCells(0);
    cluster.child<Copy>().setNumberOfCells(0);
    cluster.child<Interior>().setNumberOfCells(4*i_cells); /// Every face is a potential dynamic rupture face
  
    m_dynRupTree->allocateVariables();
    m_dynRupTree->touchVariables();
    
    m_fakeDerivativesHost = (real*) m_allocator->allocateMemory(i_cells * yateto::computeFamilySize<tensor::dQ>() * sizeof(real), PagesizeHeap, seissol::memory::Standard);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned cell = 0; cell < i_cells; ++cell) {
      for (unsigned i = 0; i < yateto::computeFamilySize<tensor::dQ>(); i++) {
        m_fakeDerivativesHost[cell*yateto::computeFamilySize<tensor::dQ>() + i] = (real)drand48();
      }
    }

#ifdef ACL_DEVICE
    m_fakeDerivatives = (real*) m_allocator->allocateMemory(i_cells * yateto::computeFamilySize<tensor::dQ>() * sizeof(real), PagesizeHeap, seissol::memory::DeviceGlobalMemory);
    const auto& device = ::device::DeviceInstance::getInstance();
    device.api->copyTo(m_fakeDerivatives, m_fakeDerivativesHost, i_cells * yateto::computeFamilySize<tensor::dQ>() * sizeof(real));
#else
    m_fakeDerivatives = m_fakeDerivativesHost;
#endif
  }

  /* cell information and integration data*/
  seissol::fakeData(m_lts, layer, (enableDynamicRupture) ? FaceType::dynamicRupture : FaceType::regular);

  if (enableDynamicRupture) {
    // From lts tree
    CellDRMapping (*drMapping)[4] = m_ltsTree->var(m_lts.drMapping);

    // From dynamic rupture tree
    seissol::initializer::Layer& interior = m_dynRupTree->child(0).child<Interior>();
    real (*imposedStatePlus)[seissol::tensor::QInterpolated::size()] = interior.var(m_dynRup.imposedStatePlus);
    real (*fluxSolverPlus)[seissol::tensor::fluxSolver::size()]     = interior.var(m_dynRup.fluxSolverPlus);
    real** timeDerivativePlus = interior.var(m_dynRup.timeDerivativePlus);
    real** timeDerivativeMinus = interior.var(m_dynRup.timeDerivativeMinus);
    DRFaceInformation* faceInformation = interior.var(m_dynRup.faceInformation);
    
    /* init drMapping */
    for (unsigned cell = 0; cell < i_cells; ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        CellDRMapping& drm = drMapping[cell][face];
        unsigned side = (unsigned int)lrand48() % 4;
        unsigned orientation = (unsigned int)lrand48() % 3;
        unsigned drFace = (unsigned int)lrand48() % interior.getNumberOfCells();
        drm.side = side;
        drm.faceRelation = orientation;
        drm.godunov = imposedStatePlus[drFace];
        drm.fluxSolver = fluxSolverPlus[drFace];
      }
    }

    /* init dr godunov state */
    for (unsigned face = 0; face < interior.getNumberOfCells(); ++face) {
      unsigned plusCell = (unsigned int)lrand48() % i_cells;
      unsigned minusCell = (unsigned int)lrand48() % i_cells;
      timeDerivativePlus[face] = &m_fakeDerivatives[plusCell * yateto::computeFamilySize<tensor::dQ>()];
      timeDerivativeMinus[face] = &m_fakeDerivatives[minusCell * yateto::computeFamilySize<tensor::dQ>()];
      
      faceInformation[face].plusSide = (unsigned int)lrand48() % 4;
      faceInformation[face].minusSide = (unsigned int)lrand48() % 4;
      faceInformation[face].faceRelation = (unsigned int)lrand48() % 3;
    }
  }
  
  return i_cells;
}

#ifdef ACL_DEVICE
void initDataStructuresOnDevice(bool enableDynamicRupture) {
  const auto& device = ::device::DeviceInstance::getInstance();
  m_ltsTree->synchronizeTo(seissol::initializer::AllocationPlace::Device, device.api->getDefaultStream());
  device.api->syncDefaultStreamWithHost();

  seissol::initializer::TimeCluster& cluster = m_ltsTree->child(0);
  seissol::initializer::Layer& layer = cluster.child<Interior>();

  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(*m_ltsTree, m_lts);
  m_ltsTree->allocateScratchPads();

  seissol::initializer::recording::CompositeRecorder<seissol::initializer::LTS> recorder;
  recorder.addRecorder(new seissol::initializer::recording::LocalIntegrationRecorder);
  recorder.addRecorder(new seissol::initializer::recording::NeighIntegrationRecorder);

  recorder.addRecorder(new seissol::initializer::recording::PlasticityRecorder);
  recorder.record(m_lts, layer);
  if (enableDynamicRupture) {
    seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForDr(*m_dynRupTree, m_dynRup);
    m_dynRupTree->allocateScratchPads();

    CompositeRecorder <seissol::initializer::DynamicRupture> drRecorder;
    drRecorder.addRecorder(new DynamicRuptureRecorder);

    auto &drLayer = m_dynRupTree->child(0).child<Interior>();
    drRecorder.record(m_dynRup, drLayer);
  }
}
#endif // ACL_DEVICE
