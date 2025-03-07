// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_RECEIVER_H_
#define SEISSOL_SRC_KERNELS_RECEIVER_H_

#include "Geometry/MeshReader.h"
#include "Initializer/PointMapper.h"
#include "Kernels/Interface.h"
#include "Kernels/Time.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Lut.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Transformation.h"
#include "Parallel/DataCollector.h"
#include "generated_code/init.h"
#include <Common/Executor.h>
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
  virtual ~DerivedReceiverQuantity() = default;
  [[nodiscard]] virtual std::vector<std::string> quantities() const = 0;
  virtual void compute(size_t sim,
                       std::vector<real>&,
                       seissol::init::QAtPoint::view::type&,
                       seissol::init::QDerivativeAtPoint::view::type&) = 0;
};

struct ReceiverRotation : public DerivedReceiverQuantity {
  ~ReceiverRotation() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override;
  void compute(size_t sim,
               std::vector<real>& /*output*/,
               seissol::init::QAtPoint::view::type& /*qAtPoint*/,
               seissol::init::QDerivativeAtPoint::view::type& /*qDerivativeAtPoint*/) override;
};

struct ReceiverStrain : public DerivedReceiverQuantity {
  ~ReceiverStrain() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override;
  void compute(size_t sim,
               std::vector<real>& /*output*/,
               seissol::init::QAtPoint::view::type& /*qAtPoint*/,
               seissol::init::QDerivativeAtPoint::view::type& /*qDerivativeAtPoint*/) override;
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

  std::vector<Receiver>::iterator begin() { return m_receivers.begin(); }

  std::vector<Receiver>::iterator end() { return m_receivers.end(); }

  [[nodiscard]] size_t ncols() const;

  void allocateData();
  void freeData();

  private:
  std::unique_ptr<seissol::parallel::DataCollector> deviceCollector{nullptr};
  std::vector<size_t> deviceIndices;
  std::vector<Receiver> m_receivers;
  seissol::kernels::Time m_timeKernel;
  std::vector<unsigned> m_quantities;
  unsigned m_nonZeroFlops{};
  unsigned m_hardwareFlops{};
  double m_samplingInterval;
  double m_syncPointInterval;
  std::vector<std::shared_ptr<DerivedReceiverQuantity>> derivedQuantities;
  seissol::SeisSol& seissolInstance;
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_RECEIVER_H_
