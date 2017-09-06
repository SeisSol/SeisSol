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

void seissol::localIntegration( struct GlobalData* globalData,
                                initializers::LTS& lts,
                                initializers::Layer& layer ) {
  kernels::Local localKernel;
  kernels::Time  timeKernel;

  real                (*dofs)[NUMBER_OF_ALIGNED_DOFS] = layer.var(lts.dofs);
  real**                buffers                       = layer.var(lts.buffers);
  LocalIntegrationData* localIntegration              = layer.var(lts.localIntegration);
  CellLocalInformation* cellInformation               = layer.var(lts.cellInformation);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    timeKernel.computeAder( 1.0,
                            globalData,
                            &localIntegration[cell],
                            dofs[cell],
                            buffers[cell],
                            nullptr );
    localKernel.computeIntegral( cellInformation[cell].faceTypes,
                                 globalData,
                                 &localIntegration[cell],
                                 buffers[cell],
                                 dofs[cell] );
  }
}

void seissol::fillWithStuff(  real* buffer,
                              unsigned nValues) {
  static unsigned long const x = 123456789, y = 362436069, z = 521288629;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned n = 0; n < nValues; ++n) {
    // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
    buffer[n] = static_cast<real>((214013*n + 2531011) / 65536);
  }
}

void seissol::fakeData( initializers::LTS& lts,
                        initializers::Layer& layer ) {
  real                (*dofs)[NUMBER_OF_ALIGNED_DOFS] = layer.var(lts.dofs);
  real**                buffers                       = layer.var(lts.buffers);
  LocalIntegrationData* localIntegration              = layer.var(lts.localIntegration);
  CellLocalInformation* cellInformation               = layer.var(lts.cellInformation);
  real*                 bucket                        = static_cast<real*>(layer.bucket(lts.buffersDerivatives));
  
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    buffers[cell] = bucket + cell * NUMBER_OF_ALIGNED_DOFS;
    
    for (unsigned f = 0; f < 4; ++f) {
      cellInformation[cell].faceTypes[f] = regular;
      cellInformation[cell].faceRelations[f][0] = ((unsigned int)lrand48() % 4);
      cellInformation[cell].faceRelations[f][1] = ((unsigned int)lrand48() % 3);
      cellInformation[cell].faceNeighborIds[f] =  ((unsigned int)lrand48() % layer.getNumberOfCells());
    }    
    cellInformation[cell].ltsSetup = 0;
  }
  
  fillWithStuff(reinterpret_cast<real*>(dofs),   NUMBER_OF_ALIGNED_DOFS * layer.getNumberOfCells());
  fillWithStuff(bucket, NUMBER_OF_ALIGNED_DOFS * layer.getNumberOfCells());
  fillWithStuff(reinterpret_cast<real*>(localIntegration), sizeof(LocalIntegrationData)/sizeof(real) * layer.getNumberOfCells());
}

double seissol::miniSeisSol(initializers::MemoryManager& memoryManager) {
  struct GlobalData* globalData = memoryManager.getGlobalData();

  initializers::LTSTree ltsTree;
  initializers::LTS     lts;
  
  lts.addTo(ltsTree);
  ltsTree.setNumberOfTimeClusters(1);
  ltsTree.fixate();
  
  initializers::TimeCluster& cluster = ltsTree.child(0);
  cluster.child<Ghost>().setNumberOfCells(0);
  cluster.child<Copy>().setNumberOfCells(0);
  cluster.child<Interior>().setNumberOfCells(50000);

  ltsTree.allocateVariables();
  ltsTree.touchVariables();
  
  initializers::Layer& layer = cluster.child<Interior>();
  
  layer.setBucketSize(lts.buffersDerivatives, sizeof(real) * NUMBER_OF_ALIGNED_DOFS * layer.getNumberOfCells());
  ltsTree.allocateBuckets();
  
  fakeData(lts, layer);
  
  localIntegration(globalData, lts, layer);
  
  Stopwatch stopwatch;
  stopwatch.start();
  for (unsigned t = 0; t < 10; ++t) {
    localIntegration(globalData, lts, layer);
  }
  return stopwatch.stop();
}
