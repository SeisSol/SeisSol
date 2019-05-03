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

void seissol::kernels::ReceiverCluster::addReceiver(  unsigned                          meshId,
                                                      unsigned                          pointId,
                                                      glm::dvec3 const&                 point,
                                                      MeshReader const&                 mesh,
                                                      seissol::initializers::Lut const& ltsLut,
                                                      seissol::initializers::LTS const& lts ) {
  auto const elements = mesh.getElements();
  auto const vertices = mesh.getVertices();

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
                            &ltsLut.lookup(lts.dofs, meshId)[0],
                            &ltsLut.lookup(lts.localIntegration, meshId),
                            reserved);
}

double seissol::kernels::ReceiverCluster::calcReceivers(  double time,
                                                          double expansionPoint,
                                                          double timeStepWidth ) {
  real timeEvaluated[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));
  real timeDerivatives[NUMBER_OF_ALIGNED_DERS] __attribute__((aligned(ALIGNMENT)));

  assert(m_global != nullptr);

  double receiverTime = time;
  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    for (auto& receiver : m_receivers) {
      m_timeKernel.computeAder( 0,
                                m_global,
                                receiver.local,
                                receiver.dofs,
                                timeEvaluated, // useless but the interface requires it
                                timeDerivatives );
      g_SeisSolNonZeroFlopsOther += m_nonZeroFlops;
      g_SeisSolHardwareFlopsOther += m_hardwareFlops;

      receiverTime = time;
      while (receiverTime < expansionPoint + timeStepWidth) {
        m_timeKernel.computeTaylorExpansion(receiverTime, expansionPoint, timeDerivatives, timeEvaluated);

        receiver.output.push_back(receiverTime);
        for (auto quantity : m_quantities) {
          auto value = receiver.basisFunctions.evalWithCoeffs(timeEvaluated + quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          receiver.output.push_back(value);
        }

        receiverTime += m_samplingInterval;
      }
    }
  }
  return receiverTime;
}

