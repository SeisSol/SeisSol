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

#include "Receiver.h"

#include <Initializer/PointMapper.h>
#include <Numerical_aux/Transformation.h>
#include <Parallel/MPI.h>
#include <Monitoring/FlopCounter.hpp>
#include <generated_code/kernel.h>

void seissol::kernels::ReceiverCluster::addReceiver(  unsigned                          meshId,
                                                      unsigned                          pointId,
                                                      Eigen::Vector3d const&            point,
                                                      MeshReader const&                 mesh,
                                                      seissol::initializers::Lut const& ltsLut,
                                                      seissol::initializers::LTS const& lts ) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  double const* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[ elements[meshId].vertices[v] ].coords;
  }
  auto xiEtaZeta = seissol::transformations::tetrahedronGlobalToReference(coords[0], coords[1], coords[2], coords[3], point);

  // (time + number of quantities) * number of samples until sync point
  size_t reserved = ncols() * (m_syncPointInterval / m_samplingInterval + 1);
  m_receivers.emplace_back( pointId,
                            xiEtaZeta[0],
                            xiEtaZeta[1],
                            xiEtaZeta[2],
                            kernels::LocalData::lookup(lts, ltsLut, meshId),
                            reserved);
}

double seissol::kernels::ReceiverCluster::calcReceivers(  double time,
                                                          double expansionPoint,
                                                          double timeStepWidth ) {
  real timeEvaluated[tensor::Q::size()] __attribute__((aligned(ALIGNMENT)));
  real timeDerivatives[yateto::computeFamilySize<tensor::dQ>()] __attribute__((aligned(ALIGNMENT)));
  real timeEvaluatedAtPoint[tensor::QAtPoint::size()] __attribute__((aligned(ALIGNMENT)));

  kernels::LocalTmp tmp;

  kernel::evaluateDOFSAtPoint krnl;
  krnl.QAtPoint = timeEvaluatedAtPoint;
  krnl.Q = timeEvaluated;

  auto qAtPoint = init::QAtPoint::view::create(timeEvaluatedAtPoint);

  double receiverTime = time;
  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    for (auto& receiver : m_receivers) {
      krnl.basisFunctionsAtPoint = receiver.basisFunctions.m_data.data();

      m_timeKernel.computeAder( timeStepWidth,
                                receiver.data,
                                tmp,
                                timeEvaluated, // useless but the interface requires it
                                timeDerivatives );
      g_SeisSolNonZeroFlopsOther += m_nonZeroFlops;
      g_SeisSolHardwareFlopsOther += m_hardwareFlops;

      receiverTime = time;
      while (receiverTime < expansionPoint + timeStepWidth) {
        m_timeKernel.computeTaylorExpansion(receiverTime, expansionPoint, timeDerivatives, timeEvaluated);
        krnl.execute();

        receiver.output.push_back(receiverTime);
#ifdef MULTIPLE_SIMULATIONS
        for (unsigned sim = init::QAtPoint::Start[0]; sim < init::QAtPoint::Stop[0]; ++sim) {
          for (auto quantity : m_quantities) {
           if (!std::isfinite(qAtPoint(sim, quantity))) {
            logError() << "Detected Inf/NaN in receiver output. Aborting.";
          }
            receiver.output.push_back(qAtPoint(sim, quantity));
          }
        }
#else
        for (auto quantity : m_quantities) {
          if (!std::isfinite(qAtPoint(quantity))) {
            logError() << "Detected Inf/NaN in receiver output. Aborting.";
          }
          receiver.output.push_back(qAtPoint(quantity));
        }
#endif

        receiverTime += m_samplingInterval;
      }
    }
  }
  return receiverTime;
}

