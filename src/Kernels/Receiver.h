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

#ifndef KERNELS_RECEIVER_H_
#define KERNELS_RECEIVER_H_

#include "Geometry/MeshReader.h"
#include "Initializer/LTS.h"
#include "Initializer/PointMapper.h"
#include "Initializer/tree/Lut.hpp"
#include "Kernels/Interface.hpp"
#include "Kernels/Time.h"
#include "Numerical_aux/BasisFunction.h"
#include "Numerical_aux/Transformation.h"
#include "Parallel/DataCollector.h"
#include "generated_code/init.h"
#include <Common/Executor.hpp>
#include <Eigen/Dense>
#include <optional>
#include <vector>

namespace seissol {
struct GlobalData;
class SeisSol;

namespace kernels {
struct Receiver {
  Receiver(unsigned pointId,
           Eigen::Vector3d position,
           const double* elementCoords[4],
           kernels::LocalData dataHost,
           kernels::LocalData dataDevice,
           size_t reserved);
  unsigned pointId;
  Eigen::Vector3d position;
  basisFunction::SampledBasisFunctions<real> basisFunctions;
  basisFunction::SampledBasisFunctionDerivatives<real> basisFunctionDerivatives;
  kernels::LocalData dataHost;
  kernels::LocalData dataDevice;
  std::vector<real> output;
};

struct DerivedReceiverQuantity {
  virtual std::vector<std::string> quantities() const = 0;
  virtual void compute(size_t sim,
                       std::vector<real>&,
                       seissol::init::QAtPoint::view::type&,
                       seissol::init::QDerivativeAtPoint::view::type&) = 0;
};

struct ReceiverRotation : public DerivedReceiverQuantity {
  std::vector<std::string> quantities() const override;
  void compute(size_t sim,
               std::vector<real>&,
               seissol::init::QAtPoint::view::type&,
               seissol::init::QDerivativeAtPoint::view::type&) override;
};

struct ReceiverStrain : public DerivedReceiverQuantity {
  std::vector<std::string> quantities() const override;
  void compute(size_t sim,
               std::vector<real>&,
               seissol::init::QAtPoint::view::type&,
               seissol::init::QDerivativeAtPoint::view::type&) override;
};

class ReceiverCluster {
  public:
  ReceiverCluster(seissol::SeisSol& seissolInstance);

  ReceiverCluster(const GlobalData* global,
                  const std::vector<unsigned>& quantities,
                  double samplingInterval,
                  double syncPointInterval,
                  const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
                  seissol::SeisSol& seissolInstance);

  void addReceiver(unsigned meshId,
                   unsigned pointId,
                   const Eigen::Vector3d& point,
                   const seissol::geometry::MeshReader& mesh,
                   const seissol::initializer::Lut& ltsLut,
                   seissol::initializer::LTS const& lts);

  //! Returns new receiver time
  double calcReceivers(
      double time, double expansionPoint, double timeStepWidth, Executor executor, void* stream);

  inline std::vector<Receiver>::iterator begin() { return m_receivers.begin(); }

  inline std::vector<Receiver>::iterator end() { return m_receivers.end(); }

  size_t ncols() const;

  void allocateData();
  void freeData();

  private:
  std::unique_ptr<seissol::parallel::DataCollector> deviceCollector{nullptr};
  std::vector<size_t> deviceIndices;
  std::vector<Receiver> m_receivers;
  seissol::kernels::Time m_timeKernel;
  std::vector<unsigned> m_quantities;
  unsigned m_nonZeroFlops;
  unsigned m_hardwareFlops;
  double m_samplingInterval;
  double m_syncPointInterval;
  std::vector<std::shared_ptr<DerivedReceiverQuantity>> derivedQuantities;
  seissol::SeisSol& seissolInstance;
};
} // namespace kernels
} // namespace seissol

#endif
