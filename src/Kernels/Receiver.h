// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_RECEIVER_H_
#define SEISSOL_SRC_KERNELS_RECEIVER_H_

#include "Common/Executor.h"
#include "GeneratedCode/init.h"
#include "Geometry/MeshReader.h"
#include "Initializer/PointMapper.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Transformation.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Runtime/Stream.h"

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
           LTS::Ref dataHost,
           LTS::Ref dataDevice,
           size_t reserved);
  unsigned pointId;
  Eigen::Vector3d position;
  basisFunction::SampledBasisFunctions<real> basisFunctions;
  basisFunction::SampledBasisFunctionDerivatives<real> basisFunctionDerivatives;
  LTS::Ref dataHost;
  LTS::Ref dataDevice;
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
  explicit ReceiverCluster(seissol::SeisSol& seissolInstance);

  ReceiverCluster(const CompoundGlobalData& global,
                  const std::vector<unsigned>& quantities,
                  double samplingInterval,
                  double syncPointInterval,
                  const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
                  seissol::SeisSol& seissolInstance);

  void addReceiver(unsigned meshId,
                   unsigned pointId,
                   const Eigen::Vector3d& point,
                   const seissol::geometry::MeshReader& mesh,
                   const LTS::Backmap& backmap);

  //! Returns new receiver time
  double calcReceivers(double time,
                       double expansionPoint,
                       double timeStepWidth,
                       Executor executor,
                       parallel::runtime::StreamRuntime& runtime);

  std::vector<Receiver>::iterator begin() { return receivers_.begin(); }

  std::vector<Receiver>::iterator end() { return receivers_.end(); }

  [[nodiscard]] size_t ncols() const;

  void allocateData();
  void freeData();

  private:
  std::optional<parallel::runtime::StreamRuntime> extraRuntime_;
  std::unique_ptr<seissol::parallel::DataCollector<real>> deviceCollector_{nullptr};
  std::vector<size_t> deviceIndices_;
  std::vector<Receiver> receivers_;
  seissol::kernels::Spacetime spacetimeKernel_;
  seissol::kernels::Time timeKernel_;
  std::vector<unsigned> quantities_;
  std::uint64_t nonZeroFlops_{};
  std::uint64_t hardwareFlops_{};
  double samplingInterval_;
  double syncPointInterval_;
  std::vector<std::shared_ptr<DerivedReceiverQuantity>> derivedQuantities_;
  seissol::SeisSol& seissolInstance_;
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_RECEIVER_H_
