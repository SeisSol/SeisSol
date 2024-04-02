/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2019, SeisSol Group
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

#include "Monitoring/FlopCounter.hpp"
#include "Numerical_aux/BasisFunction.h"
#include "Parallel/DataCollector.h"
#include "Receiver.h"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <unordered_map>

#ifdef ACL_DEVICE
#include "device.h"
#endif

void seissol::kernels::ReceiverCluster::addReceiver(  unsigned                          meshId,
                                                      unsigned                          pointId,
                                                      Eigen::Vector3d const&            point,
                                                      seissol::geometry::MeshReader const&                 mesh,
                                                      seissol::initializer::Lut const& ltsLut,
                                                      seissol::initializer::LTS const& lts ) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  double const* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[ elements[meshId].vertices[v] ].coords;
  }

  // (time + number of quantities) * number of samples until sync point
  size_t reserved = ncols() * (m_syncPointInterval / m_samplingInterval + 1);
  m_receivers.emplace_back( pointId,
                            point,
                            coords,
                            kernels::LocalData::lookup(lts, ltsLut, meshId),
                            reserved);
}

double seissol::kernels::ReceiverCluster::calcReceivers(  double time,
                                                          double expansionPoint,
                                                          double timeStepWidth ) {
  alignas(Alignment) real timeEvaluated[tensor::Q::size()];
  alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint::size()];
  alignas(Alignment) real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint::size()];
#ifdef USE_STP
  alignas(PagesizeStack) real stp[tensor::spaceTimePredictor::size()];
  kernel::evaluateDOFSAtPointSTP krnl;
  krnl.QAtPoint = timeEvaluatedAtPoint;
  krnl.spaceTimePredictor = stp;
  kernel::evaluateDerivativeDOFSAtPointSTP derivativeKrnl;
  derivativeKrnl.QDerivativeAtPoint = timeEvaluatedDerivativesAtPoint;
  derivativeKrnl.spaceTimePredictor = stp;
#else
  alignas(Alignment) real timeDerivatives[yateto::computeFamilySize<tensor::dQ>()];
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

  kernel::evaluateDOFSAtPoint krnl;
  krnl.QAtPoint = timeEvaluatedAtPoint;
  krnl.Q = timeEvaluated;
  kernel::evaluateDerivativeDOFSAtPoint derivativeKrnl;
  derivativeKrnl.QDerivativeAtPoint = timeEvaluatedDerivativesAtPoint;
  derivativeKrnl.Q = timeEvaluated;
#endif


  auto qAtPoint = init::QAtPoint::view::create(timeEvaluatedAtPoint);
  auto qDerivativeAtPoint = init::QDerivativeAtPoint::view::create(timeEvaluatedDerivativesAtPoint);

#ifdef ACL_DEVICE
  deviceCollector->gatherToHost(device::DeviceInstance::getInstance().api->getDefaultStream());
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif

  double receiverTime = time;
  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    for (size_t i = 0; i < m_receivers.size(); ++i) {
      auto& receiver = m_receivers[i];
      krnl.basisFunctionsAtPoint = receiver.basisFunctions.m_data.data();
      derivativeKrnl.basisFunctionDerivativesAtPoint = receiver.basisFunctionDerivatives.m_data.data();

      // Copy DOFs from device to host.
      LocalData tmpReceiverData { receiver.data };
#ifdef ACL_DEVICE
      tmpReceiverData.dofs_ptr = reinterpret_cast<decltype(tmpReceiverData.dofs_ptr)>(deviceCollector->get(deviceIndices[i]));
#endif
      

#ifdef USE_STP
      m_timeKernel.executeSTP(timeStepWidth, tmpReceiverData, timeEvaluated, stp);
#else
      m_timeKernel.computeAder( timeStepWidth,
                                tmpReceiverData,
                                tmp,
                                timeEvaluated, // useless but the interface requires it
                                timeDerivatives );
#endif
      seissolInstance.flopCounter().incrementNonZeroFlopsOther(m_nonZeroFlops);
      seissolInstance.flopCounter().incrementHardwareFlopsOther(m_hardwareFlops);

      receiverTime = time;
      while (receiverTime < expansionPoint + timeStepWidth) {
#ifdef USE_STP
        //eval time basis
        double tau = (time - expansionPoint) / timeStepWidth;
        seissol::basisFunction::SampledTimeBasisFunctions<real> timeBasisFunctions(ConvergenceOrder, tau);
        krnl.timeBasisFunctionsAtPoint = timeBasisFunctions.m_data.data();
        derivativeKrnl.timeBasisFunctionsAtPoint = timeBasisFunctions.m_data.data();
#else
        m_timeKernel.computeTaylorExpansion(receiverTime, expansionPoint, timeDerivatives, timeEvaluated);
#endif

        krnl.execute();
        derivativeKrnl.execute();

        receiver.output.push_back(receiverTime);
#ifdef MULTIPLE_SIMULATIONS
        for (unsigned sim = init::QAtPoint::Start[0]; sim < init::QAtPoint::Stop[0]; ++sim) {
          for (auto quantity : m_quantities) {
           if (!std::isfinite(qAtPoint(sim, quantity))) {
             logError()
                 << "Detected Inf/NaN in receiver output at"
                 << receiver.coordinates[0] << ","
                 << receiver.coordinates[1] << ","
                 << receiver.coordinates[2] << "."
                 << "Aborting.";
          }
            receiver.output.push_back(qAtPoint(sim, quantity));
          }
          if (m_computeRotation) {
            receiver.output.push_back(qDerivativeAtPoint(sim, 8, 1) - qDerivativeAtPoint(sim, 7, 2));
            receiver.output.push_back(qDerivativeAtPoint(sim, 6, 2) - qDerivativeAtPoint(sim, 8, 0));
            receiver.output.push_back(qDerivativeAtPoint(sim, 7, 0) - qDerivativeAtPoint(sim, 6, 1));
          }
        }
#else //MULTIPLE_SIMULATIONS
        for (auto quantity : m_quantities) {
          if (!std::isfinite(qAtPoint(quantity))) {
            logError()
                << "Detected Inf/NaN in receiver output at"
                << receiver.position[0] << ","
                << receiver.position[1] << ","
                << receiver.position[2] << "."
                << "Aborting.";
          }
          receiver.output.push_back(qAtPoint(quantity));
        }
        if (m_computeRotation) {
          receiver.output.push_back(qDerivativeAtPoint(8, 1) - qDerivativeAtPoint(7, 2));
          receiver.output.push_back(qDerivativeAtPoint(6, 2) - qDerivativeAtPoint(8, 0));
          receiver.output.push_back(qDerivativeAtPoint(7, 0) - qDerivativeAtPoint(6, 1));
        }
#endif //MULTITPLE_SIMULATIONS

        receiverTime += m_samplingInterval;
      }
    }
  }
  return receiverTime;
}

void seissol::kernels::ReceiverCluster::allocateData() {
#ifdef ACL_DEVICE
  // collect all data pointers to transfer. If we have multiple receivers on the same cell, we make sure to only transfer the related data once (hence, we use the `indexMap` here)
  deviceIndices.resize(m_receivers.size());
  std::vector<real*> dofs;
  std::unordered_map<real*, size_t> indexMap;
  for (size_t i = 0; i < m_receivers.size(); ++i) {
    real* currentDofs = m_receivers[i].data.dofs();
    if (indexMap.find(currentDofs) == indexMap.end()) {
      // point to the current array end
      indexMap[currentDofs] = dofs.size();
      dofs.push_back(currentDofs);
    }
    deviceIndices[i] = indexMap.at(currentDofs);
  }
  deviceCollector = std::make_unique<seissol::parallel::DataCollector>(dofs, tensor::Q::size());
#endif
}
void seissol::kernels::ReceiverCluster::freeData() {
#ifdef ACL_DEVICE
  deviceCollector.reset(nullptr);
#endif
}
