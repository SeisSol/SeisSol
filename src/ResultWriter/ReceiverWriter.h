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

#ifndef RESULTWRITER_RECEIVERWRITER_H_
#define RESULTWRITER_RECEIVERWRITER_H_

#include <vector>
#include <glm/vec3.hpp>
#include <Geometry/MeshReader.h>
#include <Numerical_aux/BasisFunction.h>
#include <Initializer/tree/Lut.hpp>
#include <Initializer/LTS.h>
#include <Kernels/Time.h>

class LocalIntegrationData;
class GlobalData;
namespace seissol {
  namespace writer {
    struct Receiver {
      Receiver(std::string const& fileName, double xi, double eta, double zeta, real* dofs, LocalIntegrationData* local)
        : fileName(fileName),
          basisFunctions(CONVERGENCE_ORDER, xi, eta, zeta),
          dofs(dofs),
          local(local)
        {}
      std::string fileName;
      BasisFunction::SampledBasisFunctions<double> basisFunctions;
      real* dofs;
      LocalIntegrationData* local;
    };

    class ReceiverWriterCluster {
    public:
      ReceiverWriterCluster()
        : m_nonZeroFlops(0), m_hardwareFlops(0)
      {}

      ReceiverWriterCluster(  GlobalData const*             global,
                              std::vector<unsigned> const&  quantities,
                              std::string const&            fileNamePrefix )
        : m_quantities(quantities), m_fileNamePrefix(fileNamePrefix) {
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
      double writeReceivers(  double time,
                              double expansionPoint,
                              double timeStepWidth,
                              double samplingInterval );

    private:
      std::vector<Receiver>   m_receivers;
      seissol::kernels::Time  m_timeKernel;
      std::vector<unsigned>   m_quantities;
      std::string             m_fileNamePrefix;
      unsigned                m_nonZeroFlops;
      unsigned                m_hardwareFlops;
    };

    class ReceiverWriter {
    public:    
      void addPoints( std::vector<glm::dvec3> const&    points,
                      MeshReader const&                 mesh,
                      seissol::initializers::Lut const& ltsLut,
                      seissol::initializers::LTS const& lts,
                      GlobalData const*                 global,
                      std::string const&                fileNamePrefix );

      //! Returns new receiver time
      double writeReceivers(  unsigned cluster,
                              double time,
                              double expansionPoint,
                              double timeStepWidth,
                              double samplingInterval ) {
        if (cluster < m_writerClusters.size()) {
          return m_writerClusters[cluster].writeReceivers(time, expansionPoint, timeStepWidth, samplingInterval);
        }
        return 0.0;
      }

    private:
      std::vector<ReceiverWriterCluster> m_writerClusters;
    };
  }
}

#endif
