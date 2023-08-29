/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * 
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 * 
 **/

#include "MiniSeisSol.h"

#include <Kernels/Time.h>
#include <Kernels/Local.h>
#include <Monitoring/Stopwatch.h>
#include "utils/env.h"

#ifdef ACL_DEVICE
#include <Initializer/BatchRecorders/Recorders.h>
#include "device.h"
#endif

namespace seissol::mini {
struct Config {
  int numRepeats{10};
  int numElements{50000};
};

Config getConfig() {
  const auto rank = seissol::MPI::mpi.rank();
  constexpr int numRepeats{10};
  constexpr int numElements{50000};

  Config config{};
  utils::Env env{};

  try {
    config.numRepeats = env.get("SEISSOL_MINI_NUM_REPEATS", numRepeats);
    if (config.numRepeats < 1) {
      throw std::runtime_error("expecting a positive integer number");
    }
  }
  catch (std::runtime_error& err) {
    logWarning(rank) << "failed to read `SEISSOL_MINI_NUM_REPEATS`," << err.what();
    config.numRepeats = numRepeats;
  }

  try {
    config.numElements = env.get("SEISSOL_MINI_NUM_ELEMENTS", numElements);
    if (config.numElements < 1) {
      throw std::runtime_error("expecting a positive integer number");
    }
  }
  catch (std::runtime_error& err) {
    logWarning(rank) << "failed to read `SEISSOL_MINI_NUM_ELEMENTS`," << err.what();
    config.numElements = numElements;
  }
  return config;
}
} // namespace seissol::mini


template<typename Config>
void seissol::localIntegration<Config>(GlobalData<Config>* globalData,
                               initializers::LTS<Config>& lts,
                               initializers::Layer& layer) {
  kernels::Local localKernel;
  localKernel.setHostGlobalData(globalData);
  kernels::Time  timeKernel;
  timeKernel.setHostGlobalData(globalData);

  typename Config::RealT**                buffers                       = layer.var(lts.buffers);

  typename kernels::LocalData<Config>::Loader loader;
  loader.load(lts, layer);
  kernels::LocalTmp<Config> tmp;

#ifdef _OPENMP
  #pragma omp parallel for private(tmp) schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    auto data = loader.entry(cell);
    timeKernel.computeAder(miniSeisSolTimeStep<Config>,
                           data,
                           tmp,
                           buffers[cell],
                           nullptr);
    localKernel.computeIntegral(buffers[cell],
                                data,
                                tmp,
                                nullptr,
                                nullptr,
                                0.0,
                                0.0);
  }
}

template<typename Config>
void seissol::fillWithStuff(  typename Config::RealT* buffer,
                              unsigned nValues) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned n = 0; n < nValues; ++n) {
    // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
    buffer[n] = static_cast<typename Config::RealT>((214013*n + 2531011) / 65536);
  }
}

template<typename Config>
void seissol::localIntegrationOnDevice<Config>(CompoundGlobalData<Config>& globalData,
                                       initializers::LTS<Config>& lts,
                                       initializers::Layer& layer) {
#ifdef ACL_DEVICE
  kernels::Time  timeKernel;
  timeKernel.setGlobalData(globalData);

  kernels::Local localKernel;
  localKernel.setGlobalData(globalData);

  const auto &device = device::DeviceInstance::getInstance();

  typename kernels::LocalData<Config>::Loader loader;
  loader.load(lts, layer);
  kernels::LocalTmp<Config> tmp;

  auto &dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto &materialTable = layer.getConditionalTable<inner_keys::Material>();
  auto &indicesTable = layer.getConditionalTable<inner_keys::Indices>();

  timeKernel.computeBatchedAder(miniSeisSolTimeStep, tmp, dataTable, materialTable, false);
  localKernel.computeBatchedIntegral(dataTable, materialTable, indicesTable, loader, tmp, 0.0);
#endif
}

template<typename Config>
void seissol::fakeData<Config>(initializers::LTS<Config>& lts,
                       initializers::Layer& layer,
                       FaceType faceTp) {
  using RealT = typename Config::RealT;
  RealT                      (*dofs)[Yateto<Config>::Tensor::Q::size()]      = layer.var(lts.dofs);
  RealT**                      buffers                       = layer.var(lts.buffers);
  RealT**                      derivatives                   = layer.var(lts.derivatives);
  RealT*                     (*faceNeighbors)[4]             = layer.var(lts.faceNeighbors);
  LocalIntegrationData*       localIntegration              = layer.var(lts.localIntegration);
  NeighboringIntegrationData* neighboringIntegration        = layer.var(lts.neighboringIntegration);
  CellLocalInformation*       cellInformation               = layer.var(lts.cellInformation);
  RealT*                       bucket                        = static_cast<RealT*>(layer.bucket(lts.buffersDerivatives));

  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    buffers[cell] = bucket + cell * Yateto<Config>::Tensor::I::size();
    derivatives[cell] = nullptr;

    for (unsigned f = 0; f < 4; ++f) {
      cellInformation[cell].faceTypes[f] = faceTp;
      cellInformation[cell].faceRelations[f][0] = ((unsigned int)lrand48() % 4);
      cellInformation[cell].faceRelations[f][1] = ((unsigned int)lrand48() % 3);
      cellInformation[cell].faceNeighborIds[f] =  ((unsigned int)lrand48() % layer.getNumberOfCells());
    }    
    cellInformation[cell].ltsSetup = 0;
  }

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {    
    for (unsigned f = 0; f < 4; ++f) {
      switch (faceTp) {
      case FaceType::freeSurface:
          faceNeighbors[cell][f] = buffers[cell];
          break;
      case FaceType::periodic:
      case FaceType::regular:
          faceNeighbors[cell][f] = buffers[ cellInformation[cell].faceNeighborIds[f] ];
          break;
        default:
          faceNeighbors[cell][f] = nullptr;
          break;
      }
    }
  }
  
  fillWithStuff(reinterpret_cast<RealT*>(dofs),   Yateto<Config>::Tensor::Q::size() * layer.getNumberOfCells());
  fillWithStuff(bucket, Yateto<Config>::Tensor::I::size() * layer.getNumberOfCells());
  fillWithStuff(reinterpret_cast<RealT*>(localIntegration), sizeof(LocalIntegrationData<Config>)/sizeof(RealT) * layer.getNumberOfCells());
  fillWithStuff(reinterpret_cast<RealT*>(neighboringIntegration), sizeof(NeighboringIntegrationData<Config>)/sizeof(RealT) * layer.getNumberOfCells());

#ifdef USE_POROELASTIC
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {    
    localIntegration[cell].specific.typicalTimeStepWidth = miniSeisSolTimeStep;
  }
#endif
}

template<typename Config>
double seissol::miniSeisSol<Config>(initializers::MemoryManager& memoryManager, bool usePlasticity) {
  initializers::LTSTree ltsTree;
  initializers::LTS<Config>     lts;

  lts.Plasticity = usePlasticity;
  lts.addTo(ltsTree);
  ltsTree.setNumberOfTimeClusters(1);
  ltsTree.fixate();

  auto config = mini::getConfig();
  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "miniSeisSol configured with"
                << config.numElements << "elements and"
                << config.numRepeats << "repeats per process";

  initializers::TimeCluster& cluster = ltsTree.child(0);
  cluster.child<Ghost>().setNumberOfCells(0);
  cluster.child<Copy>().setNumberOfCells(0);
  cluster.child<Interior>().setNumberOfCells(config.numElements);

  ltsTree.allocateVariables();
  ltsTree.touchVariables();
  
  initializers::Layer& layer = cluster.child<Interior>();
  
  layer.setBucketSize(lts.buffersDerivatives, sizeof(typename Config::RealT) * Yateto<Config>::Tensor::I::size() * layer.getNumberOfCells());
  ltsTree.allocateBuckets();

  fakeData(lts, layer);

#ifdef ACL_DEVICE
  seissol::initializers::MemoryManager::deriveRequiredScratchpadMemoryForWp(ltsTree, lts);
  ltsTree.allocateScratchPads();

  seissol::initializers::recording::CompositeRecorder<seissol::initializers::LTS> recorder;
  recorder.addRecorder(new seissol::initializers::recording::LocalIntegrationRecorder);
  recorder.addRecorder(new seissol::initializers::recording::NeighIntegrationRecorder);
  recorder.record(lts, layer);
  ltsTree.allocateScratchPads();

  auto globalData = memoryManager.getGlobalData();
  auto runBenchmark = [&globalData, &lts, &layer]() {
    localIntegrationOnDevice(globalData, lts, layer);
  };

  const auto &device = device::DeviceInstance::getInstance();
  auto syncBenchmark = [&device]() {
    device.api->syncDevice();
  };
#else
  auto* globalData = memoryManager.getGlobalDataOnHost();
  auto runBenchmark = [globalData, &lts, &layer]() {
    localIntegration(globalData, lts, layer);
  };
  auto syncBenchmark = []() {};
#endif

  runBenchmark();
  syncBenchmark();

  Stopwatch stopwatch;
  stopwatch.start();
  for (unsigned t = 0; t < config.numRepeats; ++t) {
    runBenchmark();
  }
  syncBenchmark();

  return stopwatch.stop();
}
