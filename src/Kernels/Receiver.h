// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_RECEIVER_H_
#define SEISSOL_SRC_KERNELS_RECEIVER_H_

#include "GeneratedCode/init.h"
#include "Geometry/MeshReader.h"
#include "Initializer/PointMapper.h"
#include "Kernels/Interface.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Transformation.h"
#include "Parallel/DataCollector.h"
#include <Common/Executor.h>
#include <Eigen/Dense>
#include <Initializer/Typedefs.h>
#include <Memory/Tree/Backmap.h>
#include <Parallel/Runtime/Stream.h>
#include <optional>
#include <vector>

namespace seissol {
struct GlobalData;
class SeisSol;

namespace kernels {
template <typename Cfg>
struct Receiver {
  using real = Real<Cfg>;
  Receiver(std::size_t pointId,
           Eigen::Vector3d position,
           const double* elementCoords[4],
           LTS::Ref<Cfg> dataHost,
           LTS::Ref<Cfg> dataDevice,
           size_t reserved);
  std::size_t pointId;
  Eigen::Vector3d position;
  basisFunction::SampledBasisFunctions<real> basisFunctions;
  basisFunction::SampledBasisFunctionDerivatives<Cfg> basisFunctionDerivatives;
  LTS::Ref<Cfg> dataHost;
  LTS::Ref<Cfg> dataDevice;
  std::vector<real> output;
};

template <typename Cfg>
struct DerivedReceiverQuantity {
  virtual ~DerivedReceiverQuantity() = default;
  [[nodiscard]] virtual std::vector<std::string> quantities() const = 0;
  virtual void compute(size_t sim,
                       std::vector<Real<Cfg>>&,
                       typename seissol::init::QAtPoint<Cfg>::view::type&,
                       typename seissol::init::QDerivativeAtPoint<Cfg>::view::type&) = 0;
};

template <typename Cfg>
struct ReceiverRotation : public DerivedReceiverQuantity<Cfg> {
  ~ReceiverRotation() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override {
    return {"rot1", "rot2", "rot3"};
  }
  void compute(
      size_t sim,
      std::vector<Real<Cfg>>& output,
      typename seissol::init::QAtPoint<Cfg>::view::type& qAtPoint,
      typename seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint) override {
    output.push_back(seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 8, 1) -
                     seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 7, 2));
    output.push_back(seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 6, 2) -
                     seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 8, 0));
    output.push_back(seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 7, 0) -
                     seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 6, 1));
  }
};

template <typename Cfg>
struct ReceiverStrain : public DerivedReceiverQuantity<Cfg> {
  ~ReceiverStrain() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override {
    return {"epsxx", "epsxy", "epsxz", "epsyy", "epsyz", "epszz"};
  }
  void compute(
      size_t sim,
      std::vector<Real<Cfg>>& output,
      typename seissol::init::QAtPoint<Cfg>::view::type& qAtPoint,
      typename seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint) override {
    // actually 9 quantities; 3 removed due to symmetry

    output.push_back(seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 6, 0));
    output.push_back((seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 6, 1) +
                      seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 7, 0)) /
                     2);
    output.push_back((seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 6, 2) +
                      seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 8, 0)) /
                     2);
    output.push_back(seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 7, 1));
    output.push_back((seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 7, 2) +
                      seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 8, 1)) /
                     2);
    output.push_back(seissol::multisim::multisimWrap<Cfg>(qDerivativeAtPoint, sim, 8, 2));
  }
};

class ReceiverCluster {
  public:
  virtual ~ReceiverCluster() = default;

  virtual void addReceiver(std::size_t meshId,
                           std::size_t pointId,
                           const Eigen::Vector3d& point,
                           const seissol::geometry::MeshReader& mesh,
                           const LTS::Backmap& backmap) = 0;

  //! Returns new receiver time
  virtual double calcReceivers(double time,
                               double expansionPoint,
                               double timeStepWidth,
                               Executor executor,
                               parallel::runtime::StreamRuntime& runtime) = 0;

  [[nodiscard]] virtual size_t ncols() const = 0;

  [[nodiscard]] virtual std::vector<std::string> header() const = 0;

  virtual void allocateData() = 0;
  virtual void freeData() = 0;
};

template <typename Cfg>
class ReceiverClusterImpl : public ReceiverCluster {
  public:
  using real = Real<Cfg>;
  explicit ReceiverClusterImpl(seissol::SeisSol& seissolInstance);

  ReceiverClusterImpl(
      const GlobalData& global,
      const std::vector<std::size_t>& quantities,
      double samplingInterval,
      double syncPointInterval,
      const std::vector<std::shared_ptr<DerivedReceiverQuantity<Cfg>>>& derivedQuantities,
      seissol::SeisSol& seissolInstance);

  void addReceiver(std::size_t meshId,
                   std::size_t pointId,
                   const Eigen::Vector3d& point,
                   const seissol::geometry::MeshReader& mesh,
                   const LTS::Backmap& backmap) override;

  //! Returns new receiver time
  double calcReceivers(double time,
                       double expansionPoint,
                       double timeStepWidth,
                       Executor executor,
                       parallel::runtime::StreamRuntime& runtime) override;

  typename std::vector<Receiver<Cfg>>::iterator begin() { return m_receivers.begin(); }

  typename std::vector<Receiver<Cfg>>::iterator end() { return m_receivers.end(); }

  [[nodiscard]] size_t ncols() const override;

  [[nodiscard]] std::vector<std::string> header() const override;

  void allocateData() override;
  void freeData() override;

  private:
  std::optional<parallel::runtime::StreamRuntime> extraRuntime;
  std::unique_ptr<seissol::parallel::DataCollector<real>> deviceCollector{nullptr};
  std::vector<size_t> deviceIndices;
  std::vector<Receiver<Cfg>> m_receivers;
  seissol::kernels::Spacetime<Cfg> spacetimeKernel;
  seissol::kernels::Time<Cfg> timeKernel;
  std::vector<std::size_t> m_quantities;
  std::uint64_t m_nonZeroFlops{};
  std::uint64_t m_hardwareFlops{};
  double m_samplingInterval;
  double m_syncPointInterval;
  std::vector<std::shared_ptr<DerivedReceiverQuantity<Cfg>>> derivedQuantities;
  seissol::SeisSol& seissolInstance;
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_RECEIVER_H_
