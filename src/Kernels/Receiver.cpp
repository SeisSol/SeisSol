/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include "Receiver.h"
#include "Monitoring/FlopCounter.hpp"
#include "Numerical_aux/BasisFunction.h"
#include "Parallel/DataCollector.h"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <Common/Executor.hpp>
#include <Initializer/tree/Layer.hpp>
#include <Kernels/Plasticity.h>
#include <Kernels/common.hpp>
#include <Kernels/Plasticity.h>
#include <cstddef>
#include <cstring>
#include <init.h>
#include <omp.h>
#include <string>
#include <tensor.h>
#include <unordered_map>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace {
#ifdef MULTIPLE_SIMULATIONS
template <typename T, typename F, typename... Args>
T multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), sim, std::forward<Args>(args)...);
}
constexpr size_t MultisimStart = init::QAtPoint::Start[0];
constexpr size_t MultisimEnd = init::QAtPoint::Stop[0];
#else
template <typename F, typename... Args>
auto multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), std::forward<Args>(args)...);
}
constexpr size_t MultisimStart = 0;
constexpr size_t MultisimEnd = 1;
#endif
} // namespace

namespace seissol::kernels {

Receiver::Receiver(unsigned pointId,
                   Eigen::Vector3d position,
                   const double* elementCoords[4],
                   kernels::LocalData dataHost,
                   kernels::LocalData dataDevice,
                   size_t reserved)
    : pointId(pointId), position(std::move(position)), dataHost(dataHost), dataDevice(dataDevice) {
  output.reserve(reserved);

  auto xiEtaZeta = seissol::transformations::tetrahedronGlobalToReference(
      elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], position);
  basisFunctions = basisFunction::SampledBasisFunctions<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives = basisFunction::SampledBasisFunctionDerivatives<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives.transformToGlobalCoordinates(elementCoords);
}

ReceiverCluster::ReceiverCluster(seissol::SeisSol& seissolInstance)
    : m_nonZeroFlops(0), m_hardwareFlops(0), m_samplingInterval(1.0e99), m_syncPointInterval(0.0),
      seissolInstance(seissolInstance) {}

ReceiverCluster::ReceiverCluster(
    const GlobalData* global,
    const std::vector<unsigned>& quantities,
    double samplingInterval,
    double syncPointInterval,
    const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
    seissol::SeisSol& seissolInstance)
    : m_quantities(quantities), m_samplingInterval(samplingInterval),
      m_syncPointInterval(syncPointInterval), derivedQuantities(derivedQuantities),
      seissolInstance(seissolInstance),
      m_globalData(global) {
  m_timeKernel.setHostGlobalData(global);
  m_timeKernel.flopsAder(m_nonZeroFlops, m_hardwareFlops);
}

void ReceiverCluster::addReceiver(unsigned meshId,
                                  unsigned pointId,
                                  const Eigen::Vector3d& point,
                                  const seissol::geometry::MeshReader& mesh,
                                  const seissol::initializer::Lut& ltsLut,
                                  const seissol::initializer::LTS& lts) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  const double* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[elements[meshId].vertices[v]].coords;
  }

  // (time + number of quantities) * number of samples until sync point
  size_t reserved = ncols() * (m_syncPointInterval / m_samplingInterval + 1);
  m_receivers.emplace_back(
      pointId,
      point,
      coords,
      kernels::LocalData::lookup(lts, ltsLut, meshId, initializer::AllocationPlace::Host),
      kernels::LocalData::lookup(lts,
                                 ltsLut,
                                 meshId,
                                 isDeviceOn() ? initializer::AllocationPlace::Device
                                              : initializer::AllocationPlace::Host),
      reserved);
}

double ReceiverCluster::calcReceivers(
    double time, double expansionPoint, double timeStepWidth, Executor executor, void* stream) {
  std::size_t ncols = this->ncols();

  double outReceiverTime = time;
  while (outReceiverTime < expansionPoint + timeStepWidth) {
    outReceiverTime += m_samplingInterval;
  }

#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    deviceCollector->gatherToHost(device::DeviceInstance::getInstance().api->getDefaultStream());
    deviceCollectorPlasticity->gatherToHost(
        device::DeviceInstance::getInstance().api->getDefaultStream());
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif
  auto tv = seissolInstance.getSeisSolParameters().model.tv;
  auto oneMinusIntegratingFactor = (tv > 0.0) ? 1.0 - exp(-timeStepWidth / tv) : 1.0;
  auto globalData = seissolInstance.getMemoryManager().getGlobalData().onHost;

  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    // heuristic; to avoid the overhead from the parallel region
    auto threshold = std::max(1000, omp_get_num_threads() * 100);
    auto recvCount = m_receivers.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (recvCount >= threshold)
#endif
    for (size_t i = 0; i < recvCount; ++i) {
      alignas(Alignment) real timeEvaluated[tensor::Q::size()];
      alignas(Alignment) real timeEvaluatedPrev[tensor::Q::size()];
      alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint::size()];
      alignas(Alignment) real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint::size()];
      alignas(Alignment) real QEtaModal[tensor::QEtaModal::size()];
      alignas(Alignment) real QEtaNodal[tensor::QEtaNodal::size()];
      alignas(Alignment) real dudt_pstrain[tensor::QStress::size()];
      alignas(Alignment) real QStressNodal[tensor::QStressNodal::size()];
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
      auto qDerivativeAtPoint =
          init::QDerivativeAtPoint::view::create(timeEvaluatedDerivativesAtPoint);

      auto& receiver = m_receivers[i];
      krnl.basisFunctionsAtPoint = receiver.basisFunctions.m_data.data();
      derivativeKrnl.basisFunctionDerivativesAtPoint =
          receiver.basisFunctionDerivatives.m_data.data();

      // Copy DOFs from device to host.
      LocalData tmpReceiverData{receiver.dataHost};

      // we always need a pstrain copy here (the dofs will get overridden by the predictor again,
      // but that's not true for the pstrain)
      alignas(Alignment)
          real pstrain[seissol::tensor::QStress::size() + seissol::tensor::QEtaModal::size()];
          
      if (executor == Executor::Device) {
#ifdef ACL_DEVICE
        tmpReceiverData.dofs_ptr = reinterpret_cast<decltype(tmpReceiverData.dofs_ptr)>(
            deviceCollector->get(deviceIndices[i]));
        std::memcpy(pstrain,
                    reinterpret_cast<decltype(tmpReceiverData.pstrain_ptr)>(
                        deviceCollectorPlasticity->get(deviceIndices[i])),
                    sizeof(pstrain));
#endif
      } else {
        std::memcpy(pstrain, receiver.dataHost.pstrain_ptr, sizeof(pstrain));
      }

#ifdef USE_STP
      m_timeKernel.executeSTP(timeStepWidth, tmpReceiverData, timeEvaluated, stp);
#else
      m_timeKernel.computeAder(m_samplingInterval,
                               tmpReceiverData,
                               tmp,
                               timeEvaluated, // useless but the interface requires it
                               timeDerivatives);
#endif

      seissolInstance.flopCounter().incrementNonZeroFlopsOther(m_nonZeroFlops);
      seissolInstance.flopCounter().incrementHardwareFlopsOther(m_hardwareFlops);

      double receiverTime = time;
      kernel::evaluatePstrainAtPoint pstrainKrnl;
      pstrainKrnl.basisFunctionsAtPoint = receiver.basisFunctions.m_data.data();
      pstrainKrnl.Pstrain = pstrain;
      alignas(Alignment) real pStrainAtPoint[tensor::PstrainAtPoint::size()];
      auto pAtPoint = init::QAtPoint::view::create(pStrainAtPoint);
      pstrainKrnl.PstrainAtPoint = pStrainAtPoint;

      // m_timeKernel.computeTaylorExpansion(
      //       expansionPoint, expansionPoint, timeDerivatives, timeEvaluatedPrev);

      while (receiverTime < expansionPoint + timeStepWidth) {
#ifdef USE_STP
        // eval time basis
        double tau = (time - expansionPoint) / timeStepWidth;
        seissol::basisFunction::SampledTimeBasisFunctions<real> timeBasisFunctions(ConvergenceOrder,
                                                                                   tau);
        krnl.timeBasisFunctionsAtPoint = timeBasisFunctions.m_data.data();
        derivativeKrnl.timeBasisFunctionsAtPoint = timeBasisFunctions.m_data.data();
#else
        m_timeKernel.computeTaylorExpansion(
            receiverTime, expansionPoint, timeDerivatives, timeEvaluated);
#endif

        Plasticity::computePlasticity(oneMinusIntegratingFactor,
                                      m_samplingInterval,
                                      tv,
                                      globalData,
                                      &tmpReceiverData.plasticity(),
                                      timeEvaluated,
                                      pstrain);
      // for (unsigned q = 0; q < tensor::QStress::size(); ++q) {
      //   /**
      //    * Equation (10) from Wollherr et al.:
      //    *
      //    * d/dt strain_{ij} = (sigma_{ij} + sigma0_{ij} - P_{ij}(sigma)) / (2mu T_v)
      //    *
      //    * where (11)
      //    *
      //    * P_{ij}(sigma) = { tau_c/tau s_{ij} + m delta_{ij}         if     tau >= taulim
      //    *                 { sigma_{ij} + sigma0_{ij}                else
      //    *
      //    * Thus,
      //    *
      //    * d/dt strain_{ij} = { (1 - tau_c/tau) / (2mu T_v) s_{ij}   if     tau >= taulim
      //    *                    { 0                                    else
      //    *
      //    * Consider tau >= taulim first. We have (1 - tau_c/tau) = -yield / r. Therefore,
      //    *
      //    * d/dt strain_{ij} = -1 / (2mu T_v r) yield s_{ij}
      //    *                  = -1 / (2mu T_v r) (sigmaNew_{ij} - sigma_{ij})
      //    *                  = (sigma_{ij} - sigmaNew_{ij}) / (2mu T_v r)
      //    *
      //    * If tau < taulim, then sigma_{ij} - sigmaNew_{ij} = 0.
      //    */
      //   real factor = tmpReceiverData.plasticity().mufactor / (tv * oneMinusIntegratingFactor);
      //   dudt_pstrain[q] = factor * (timeEvaluatedPrev[q] - timeEvaluated[q]);
      //   // Integrate with explicit Euler
      //   pstrain[q] += m_samplingInterval * dudt_pstrain[q];
      // }

      // kernel::plConvertToNodalNoLoading m2nKrnl_dudt_pstrain;
      // m2nKrnl_dudt_pstrain.v = m_globalData->vandermondeMatrix;
      // m2nKrnl_dudt_pstrain.QStress = dudt_pstrain;
      // m2nKrnl_dudt_pstrain.QStressNodal = QStressNodal;
      // m2nKrnl_dudt_pstrain.execute();

      // for (unsigned q = 0; q < tensor::QEtaModal::size(); ++q) {
      //   QEtaModal[q] = pstrain[tensor::QStress::size() + q];
      // }

      // kernel::plConvertEtaModal2Nodal m2n_eta_Krnl;
      // m2n_eta_Krnl.v = m_globalData->vandermondeMatrix;
      // m2n_eta_Krnl.QEtaModal = QEtaModal;
      // m2n_eta_Krnl.QEtaNodal = QEtaNodal;
      // m2n_eta_Krnl.execute();

      // auto QStressNodalView = init::QStressNodal::view::create(QStressNodal);
      // unsigned numNodes = QStressNodalView.shape(0);
      // for (unsigned i = 0; i < numNodes; ++i) {
      //   // eta := int_0^t sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt) dt
      //   // Approximate with eta += timeStepWidth * sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt)
      //   QEtaNodal[i] = std::max((real) 0.0, QEtaNodal[i]) + 
      //                  timeStepWidth * sqrt(0.5 * (QStressNodalView(i, 0) * QStressNodalView(i, 0)  + QStressNodalView(i, 1) * QStressNodalView(i, 1)
      //                                             + QStressNodalView(i, 2) * QStressNodalView(i, 2)  + QStressNodalView(i, 3) * QStressNodalView(i, 3)
      //                                             + QStressNodalView(i, 4) * QStressNodalView(i, 4)  + QStressNodalView(i, 5) * QStressNodalView(i, 5)));
      // }
      // kernel::plConvertEtaNodal2Modal n2m_eta_Krnl;
      // n2m_eta_Krnl.vInv = m_globalData->vandermondeMatrixInverse;
      // n2m_eta_Krnl.QEtaNodal = QEtaNodal;
      // n2m_eta_Krnl.QEtaModal = QEtaModal;
      // n2m_eta_Krnl.execute();

      // for (unsigned q = 0; q < tensor::QEtaModal::size(); ++q) {
      //   pstrain[tensor::QStress::size() + q] = QEtaModal[q];
      // }

        krnl.execute();
        derivativeKrnl.execute();
        pstrainKrnl.execute();

        // note: necessary receiver space is reserved in advance
        receiver.output.push_back(receiverTime);
        for (unsigned sim = MultisimStart; sim < MultisimEnd; ++sim) {
          for (auto quantity : m_quantities) {
            if (!std::isfinite(multisimWrap(qAtPoint, sim, quantity))) {
              logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0] << ","
                         << receiver.position[1] << "," << receiver.position[2] << " in simulation"
                         << sim << "."
                         << "Aborting.";
            }
            receiver.output.push_back(multisimWrap(qAtPoint, sim, quantity));
          }
          for (const auto& derived : derivedQuantities) {
            derived->compute(sim, receiver.output, qAtPoint, qDerivativeAtPoint);
          }
          for (unsigned int i=0; i < 7; i++){
            receiver.output.push_back(multisimWrap(pAtPoint, sim, i));
          }
        }

        receiverTime += m_samplingInterval;
        std::memcpy(timeEvaluatedPrev, timeEvaluated, sizeof(timeEvaluatedPrev));
      }
    }
  }
  return outReceiverTime;
}

void ReceiverCluster::allocateData() {
#ifdef ACL_DEVICE
  // collect all data pointers to transfer. If we have multiple receivers on the same cell, we make
  // sure to only transfer the related data once (hence, we use the `indexMap` here)
  deviceIndices.resize(m_receivers.size());
  std::vector<real*> dofs;
  std::vector<real*> pstrain;
  std::unordered_map<real*, size_t> indexMap;
  for (size_t i = 0; i < m_receivers.size(); ++i) {
    real* currentDofs = m_receivers[i].dataDevice.dofs();
    if (indexMap.find(currentDofs) == indexMap.end()) {
      // point to the current array end
      indexMap[currentDofs] = dofs.size();
      dofs.push_back(currentDofs);
      pstrain.push_back(m_receivers[i].dataDevice.pstrain());
    }
    deviceIndices[i] = indexMap.at(currentDofs);
  }
  deviceCollector =
      std::make_unique<seissol::parallel::DataCollector>(dofs, tensor::Q::size(), useUSM());
  deviceCollectorPlasticity = std::make_unique<seissol::parallel::DataCollector>(
      pstrain, seissol::tensor::QStress::size() + seissol::tensor::QEtaModal::size(), useUSM());
#endif
}
void ReceiverCluster::freeData() {
#ifdef ACL_DEVICE
  deviceCollector.reset(nullptr);
  deviceCollectorPlasticity.reset(nullptr);
#endif
}

size_t ReceiverCluster::ncols() const {
  size_t ncols = m_quantities.size();
  for (const auto& derived : derivedQuantities) {
    ncols += derived->quantities().size();
  }
  ncols *= MultisimEnd - MultisimStart;
  return 1 + ncols;
}

std::vector<std::string> ReceiverRotation::quantities() const { return {"rot1", "rot2", "rot3"}; }
void ReceiverRotation::compute(size_t sim,
                               std::vector<real>& output,
                               seissol::init::QAtPoint::view::type& qAtPoint,
                               seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint) {
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 8, 1) -
                   multisimWrap(qDerivativeAtPoint, sim, 7, 2));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6, 2) -
                   multisimWrap(qDerivativeAtPoint, sim, 8, 0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 7, 0) -
                   multisimWrap(qDerivativeAtPoint, sim, 6, 1));
}

std::vector<std::string> ReceiverStrain::quantities() const {
  return {"epsxx", "epsxy", "epsxz", "epsyy", "epsyz", "epszz"};
}

std::vector<std::string> ReceiverPlasticStrain::quantities() const {
  return {"plasticepsxx", "plasticepsxy", "plasticepsxz", "plasticepsyy", "plasticepsyz", "plasticepszz"};
}

void ReceiverPlasticStrain::compute(size_t sim, std::vector<real>& output, seissol::init::QAtPoint::view::type& qAtPoint,
seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint){
  // maintaing 6 as of now, will need to extend once we get confirmation that it is more
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6,0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6,0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6,0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6,0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6,0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6,0));

}

void ReceiverStrain::compute(size_t sim,
                             std::vector<real>& output,
                             seissol::init::QAtPoint::view::type& qAtPoint,
                             seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint) {
  // actually 9 quantities; 3 removed due to symmetry

  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6, 0));
  output.push_back(
      (multisimWrap(qDerivativeAtPoint, sim, 6, 1) + multisimWrap(qDerivativeAtPoint, sim, 7, 0)) /
      2);
  output.push_back(
      (multisimWrap(qDerivativeAtPoint, sim, 6, 2) + multisimWrap(qDerivativeAtPoint, sim, 8, 0)) /
      2);
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 7, 1));
  output.push_back(
      (multisimWrap(qDerivativeAtPoint, sim, 7, 2) + multisimWrap(qDerivativeAtPoint, sim, 8, 1)) /
      2);
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 8, 2));
}

} // namespace seissol::kernels
