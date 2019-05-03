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

#ifndef KERNELS_RECEIVER_H_
#define KERNELS_RECEIVER_H_

#include <vector>
#include <glm/vec3.hpp>
#include <Geometry/MeshReader.h>
#include <Numerical_aux/BasisFunction.h>
#include <Initializer/tree/Lut.hpp>
#include <Initializer/LTS.h>
#include <Initializer/PointMapper.h>
#include <Kernels/Time.h>
#include <Kernels/Interface.hpp>
#include <generated_code/init.h>

class GlobalData;
namespace seissol {
  namespace kernels {
    struct Receiver {
      Receiver(unsigned pointId, double xi, double eta, double zeta, kernels::LocalData data, size_t reserved)
        : pointId(pointId),
          basisFunctions(CONVERGENCE_ORDER, xi, eta, zeta),
          data(data)
      {
        output.reserve(reserved);
      }
      unsigned pointId;
      basisFunction::SampledBasisFunctions<real> basisFunctions;
      kernels::LocalData data;
      std::vector<real> output;
    };

    class ReceiverCluster {
    public:
      ReceiverCluster()
        : m_nonZeroFlops(0), m_hardwareFlops(0),
          m_samplingInterval(1.0e99), m_syncPointInterval(0.0)
      {}

      ReceiverCluster(  GlobalData const*             global,
                        std::vector<unsigned> const&  quantities,
                        double                        samplingInterval,
                        double                        syncPointInterval )
        : m_quantities(quantities),
          m_samplingInterval(samplingInterval), m_syncPointInterval(syncPointInterval) {
        m_timeKernel.setGlobalData(global);
        m_timeKernel.flopsAder(m_nonZeroFlops, m_hardwareFlops);
      }

      void addReceiver( unsigned          meshId,
                        unsigned          pointId,
                        glm::dvec3 const& point,
                        MeshReader const& mesh,
                        seissol::initializers::Lut const& ltsLut,
                        seissol::initializers::LTS const& lts );

      //! Returns new receiver time
      double calcReceivers( double time,
                            double expansionPoint,
                            double timeStepWidth );

      std::vector<Receiver>::iterator begin() {
        return m_receivers.begin();
      }

      std::vector<Receiver>::iterator end() {
        return m_receivers.end();
      }

      size_t ncols() const {
        size_t ncols = m_quantities.size();
#ifdef MULTIPLE_SIMULATIONS
        ncols *= init::QAtPoint::Stop[0]-init::QAtPoint::Start[0];
#endif
        return 1 + ncols;
      }

    private:
      std::vector<Receiver>   m_receivers;
      seissol::kernels::Time  m_timeKernel;
      std::vector<unsigned>   m_quantities;
      unsigned                m_nonZeroFlops;
      unsigned                m_hardwareFlops;
      double                  m_samplingInterval;
      double                  m_syncPointInterval;
    };
  }
}

#endif
