// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_BENCH_SYNTHETICDATA_H_
#define SEISSOL_SRC_SIMULATION_BENCH_SYNTHETICDATA_H_

#include "Core/Common/Real.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Setup/Factory.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Simulation/Bench/LoopBench.h"
#include "Simulation/Bench/Proxy.h"

#include <cstddef>
#include <memory>

namespace seissol::proxy {

/// the size of the time steps of the loops: a step of a cell of some hundred meters, as the
/// friction laws take it to be one
constexpr double Timestep = 1e-3;

constexpr double GravitationalAcceleration = 9.81;

/// the relaxation time of the plasticity: long against the time step, so that the stress of a cell
/// that yields approaches the yield surface only slowly, and the cell yields in every step (for a
/// million steps, the stress beyond the surface still goes down to a third)
constexpr double RelaxationTime = 1e6 * Timestep;

/// a global time step from time 0
constexpr solver::StepParams Step{.time = 0,
                                  .timeStepSize = Timestep,
                                  .subTimeStart = 0,
                                  .neighborTimeStepSize = Timestep,
                                  .resetBuffers = true};

/**
 * The data the loops of the proxy run on, in the configuration `Cfg`: an interior layer of cells
 * with random neighbors and arbitrary values. Either all faces of the cells are regular, or all of
 * them lie on the fault, with four times as many fault faces as cells, whose adjacent cells are
 * random as well. The cells keep the time integrals of their steps in their buffers, which their
 * neighbors read (GTS).
 *
 * With plasticity, the cells either stay within the yield surface, or have neither cohesion nor
 * friction and yield in every step (see RelaxationTime). The faults have the friction law of the
 * setup, with the parameters of a benchmark of it, the same at all points (see
 * docs/miniapp.rst), and the stress sources of a run by default; what the friction law derives
 * from these, its initializer derives (see BaseDRInitializer::deriveInitialState()).
 */
template <typename Cfg>
class SyntheticData {
  public:
  using real = Real<Cfg>;

  explicit SyntheticData(const BenchSetup& setup);

  SyntheticData(const SyntheticData&) = delete;
  SyntheticData(SyntheticData&&) = delete;
  auto operator=(const SyntheticData&) = delete;
  auto operator=(SyntheticData&&) = delete;
  ~SyntheticData() = default;

  [[nodiscard]] LTS::Layer& cells() { return ltsStorage_.layer(layerId_); }
  [[nodiscard]] const LTS::Layer& cells() const { return ltsStorage_.layer(layerId_); }

  /// the fault faces; only with dynamic rupture
  [[nodiscard]] DynamicRupture::Layer& faces() { return drStorage_.layer(layerId_); }
  [[nodiscard]] const DynamicRupture::Layer& faces() const { return drStorage_.layer(layerId_); }

  [[nodiscard]] bool hasFaces() const { return dynamicRupture_; }
  [[nodiscard]] CellPlasticity plasticity() const { return plasticity_; }

  [[nodiscard]] GlobalDataStore& globalData() { return globalData_; }

  /// the friction solvers of the faults on the host and on the device; only with dynamic rupture
  [[nodiscard]] const dr::friction_law::FrictionSolverFactory& frictionLaw() const {
    return frictionLaw_;
  }
  [[nodiscard]] const dr::friction_law::FrictionSolverFactory& frictionLawDevice() const {
    return frictionLawDevice_;
  }

  private:
  void initGlobalData();
  void initCells(FaceType faceType);
  void initPlasticity();
  void initFaults(FrictionLaw frictionLaw);
  void initFaultState(dr::factory::AbstractFactory& factory, FrictionLaw frictionLaw);
  void initDevice();

  std::size_t cellCount_;
  bool dynamicRupture_;
  CellPlasticity plasticity_;
  initializer::LayerIdentifier layerId_;
  //! the parameters of the faults, which their storage, their friction law and its initializer
  //! take alike
  std::shared_ptr<initializer::parameters::DRParameters> drParameters_;
  std::unique_ptr<DynamicRupture> faultVariables_;
  dr::friction_law::FrictionSolverFactory frictionLaw_;
  dr::friction_law::FrictionSolverFactory frictionLawDevice_;

  memory::ManagedAllocator allocator_;
  GlobalDataStore globalData_;
  LTS::Storage ltsStorage_;
  DynamicRupture::Storage drStorage_;
  //! the derivatives of random cells, which the fault faces read
  real* faceDerivatives_{nullptr};
  real* faceDerivativesOnDevice_{nullptr};
};

} // namespace seissol::proxy

#endif // SEISSOL_SRC_SIMULATION_BENCH_SYNTHETICDATA_H_
