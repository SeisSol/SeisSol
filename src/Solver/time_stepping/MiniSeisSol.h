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

#ifndef MINISEISSOL_H_
#define MINISEISSOL_H_

#include "MiniSeisSolTypes.h"
#include <Initializer/MemoryManager.h>

namespace seissol {
  void localIntegration(  struct GlobalData* globalData,
                          initializers::LTS& lts,
                          initializers::Layer& layer );
  
  void fillWithStuff( real* buffer,
                      unsigned nValues );

  void fakeData(  initializers::LTS& lts,
                  initializers::Layer& layer,
                  FaceType faceTp = FaceType::regular);

  inline void fakeData(ProxyData& data, FaceType faceType = FaceType::regular) {
    auto elementViewFactory = mneme::createViewFactory().withPlan(data.elementStoragePlan).withStorage(data.elementStorage);
    auto elementViewInterior = elementViewFactory.createDenseView<InteriorLayer>();

    auto displacementBucketViewFactory = mneme::createViewFactory().withPlan(data.buffersBucketPlan).withStorage(data.buffersBucket);
    auto displacementBucketViewInterior = displacementBucketViewFactory.withStride<tensor::I::size()>().createStridedView<InteriorLayer>();

    for (unsigned cell = 0; cell < elementViewInterior.size(); ++cell) {
      auto& curElement = elementViewInterior[cell];
      //auto& dofs = elementViewInterior[cell].get<dofs>();
      // TODO(Lukas) Is the offset for the bucket correct?
      curElement.get<buffer>() = &displacementBucketViewInterior[cell][0];
      auto& cellInformation = curElement.get<cellLocalInformation>();
      for (unsigned f = 0; f < 4; ++f) {
        cellInformation.faceTypes[f] = faceType;
        cellInformation.faceRelations[f][0] = ((unsigned int) lrand48() % 4);
        cellInformation.faceRelations[f][1] = ((unsigned int) lrand48() % 3);
        cellInformation.faceNeighborIds[f] = ((unsigned int) lrand48() % elementViewInterior.size());
      }
      cellInformation.ltsSetup = 0;
    }
#pragma omp parallel for schedule(static) default(none) shared(elementViewInterior, faceType)
    for (unsigned cell = 0; cell < elementViewInterior.size(); ++cell) {
      auto& curElement = elementViewInterior[cell];
      auto& curCellInformation = curElement.get<cellLocalInformation>();
      auto& curFaceNeighbors = curElement.get<faceNeighbors>();
      auto* curBuffer = curElement.get<buffer>();
      for (unsigned f = 0; f < 4; ++f) {
        switch (faceType) {
          case FaceType::freeSurface:
            curFaceNeighbors[f] = curBuffer;
            break;
          case FaceType::periodic:
            [[fallthrough]];
          case FaceType::regular:
            curFaceNeighbors[f] = elementViewInterior[curCellInformation.faceNeighborIds[f]].get<buffer>();
            break;
          default:
            curFaceNeighbors[f] = nullptr;
            break;
        }
      }
      fillWithStuff(curElement.get<dofs>().data(), tensor::Q::size());

      //fillWithStuff(reinterpret_cast<real*>(dofs),   tensor::Q::size() * layer.getNumberOfCells());

      fillWithStuff(curElement.get<buffer>(), tensor::I::size());
      //fillWithStuff(bucket, tensor::I::size() * layer.getNumberOfCells());

      fillWithStuff(reinterpret_cast<real*>(&curElement.get<localIntegrationData>()), sizeof(LocalIntegrationData)/sizeof(real));
      //fillWithStuff(reinterpret_cast<real*>(localIntegration), sizeof(LocalIntegrationData)/sizeof(real) * layer.getNumberOfCells());

      fillWithStuff(reinterpret_cast<real*>(&curElement.get<neighborIntegrationData>()), sizeof(NeighboringIntegrationData)/sizeof(real));
      //fillWithStuff(reinterpret_cast<real*>(neighboringIntegration), sizeof(NeighboringIntegrationData)/sizeof(real) * layer.getNumberOfCells());
    }
  }
  
  double miniSeisSol(initializers::MemoryManager& memoryManager);
}


#endif // MINISEISSOL_H_
