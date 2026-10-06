// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Alignment.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Plasticity.h"
#include "Model/CommonDatastructures.h"
#include "Model/Plasticity.h"
#include "Solver/MultipleSimulations.h"
#include "TestConfigs.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::unit_test {

// ---------------------------------------------------------------------------
// computeRelaxTime: pure constexpr math
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("Plasticity computeRelaxTime" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using Plasticity = seissol::kernels::Plasticity<Cfg>;
  SUBCASE("tV = 0 → factor = 1 (instantaneous relaxation)") {
    CHECK(Plasticity::computeRelaxTime(0.0, 1.0) == doctest::Approx(1.0));
    CHECK(Plasticity::computeRelaxTime(0.0, 0.5) == doctest::Approx(1.0));
    CHECK(Plasticity::computeRelaxTime(0.0, 0.001) == doctest::Approx(1.0));
  }

  SUBCASE("tV > 0 → exponential decay factor") {
    // computeRelaxTime(tV, dt) = -expm1(-dt/tV) = 1 - exp(-dt/tV)
    const double tV = 0.1;
    const double dt = 0.05;
    const double expected = 1.0 - std::exp(-dt / tV);
    CHECK(Plasticity::computeRelaxTime(tV, dt) == doctest::Approx(expected));
  }

  SUBCASE("Small dt/tV ratio (nearly no relaxation)") {
    const double tV = 100.0;
    const double dt = 0.001;
    double result = Plasticity::computeRelaxTime(tV, dt);
    // For small dt/tV: result ≈ dt/tV
    CHECK(result == doctest::Approx(dt / tV).epsilon(1e-4));
    CHECK(result > 0.0);
    CHECK(result < 1.0);
  }

  SUBCASE("Large dt/tV ratio (full relaxation)") {
    const double tV = 0.001;
    const double dt = 10.0;
    double result = Plasticity::computeRelaxTime(tV, dt);
    // For large dt/tV: result ≈ 1
    CHECK(result == doctest::Approx(1.0).epsilon(1e-6));
  }

  SUBCASE("dt = tV → known value") {
    const double tV = 1.0;
    const double dt = 1.0;
    // 1 - exp(-1) ≈ 0.6321
    CHECK(Plasticity::computeRelaxTime(tV, dt) == doctest::Approx(1.0 - std::exp(-1.0)));
  }

  SUBCASE("Result is always in (0, 1] for tV > 0") {
    for (const double tV : {0.01, 0.1, 1.0, 10.0, 100.0}) {
      for (const double dt : {0.001, 0.01, 0.1, 1.0, 10.0}) {
        double result = Plasticity::computeRelaxTime(tV, dt);
        CHECK(result > 0.0);
        CHECK(result <= 1.0);
      }
    }
  }

  SUBCASE("Monotone in dt for fixed tV") {
    const double tV = 1.0;
    double prev = 0.0;
    for (const double dt : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) {
      double result = Plasticity::computeRelaxTime(tV, dt);
      CHECK(result > prev);
      prev = result;
    }
  }
}

// ---------------------------------------------------------------------------
// flopsPlasticity: generated code constants
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("Plasticity metrics" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using Plasticity = seissol::kernels::Plasticity<Cfg>;
  const auto [metricsCheck, metricsYield] = Plasticity::metrics();

  SUBCASE("Check flops are positive") {
    CHECK(metricsCheck.nonzeroFlop > 0);
    CHECK(metricsCheck.hardwareFlop > 0);
  }

  SUBCASE("Yield flops are positive") {
    CHECK(metricsYield.nonzeroFlop > 0);
    CHECK(metricsYield.hardwareFlop > 0);
  }

  SUBCASE("Hardware flops >= nonzero flops") {
    CHECK(metricsCheck.hardwareFlop >= metricsCheck.nonzeroFlop);
    CHECK(metricsYield.hardwareFlop >= metricsYield.nonzeroFlop);
  }
}

// ---------------------------------------------------------------------------
// computePlasticity: which cells get the plastic correction
// ---------------------------------------------------------------------------

// One cell under uniform pure shear: the modal DOFs are zero, and the initial loading is
// sigma_xy = ShearStress at every node. The mean stress then vanishes and tau = ShearStress at
// every node, and with zero bulk friction the yield stress of a node is its cohesion. Each node is
// given a cohesion of either half or twice ShearStress, so it is chosen exactly which nodes yield.
template <typename Cfg>
class ShearedPlasticityCell {
  public:
  using Plasticity = seissol::kernels::Plasticity<Cfg>;
  using real = Real<Cfg>;

  static constexpr std::size_t NumNodes = model::PlasticityData<Cfg>::PointCount;
  static constexpr std::size_t ComponentXY = 3;

  static constexpr double ShearStress = 1.0e6;
  static constexpr double Mu = 3.0e10;
  static constexpr double RelaxationTime = 0.1;
  static constexpr double TimeStep = 1.0e-3;

  // At a yielding node, the yield factor is (taulim / tau - 1) r = -r / 2, hence (Wollherr et al.,
  // eq. 10) d/dt strain_xy = -yield s_xy / (2 mu tV r) = ShearStress / (4 mu tV), while the other
  // components stay zero. eta grows by dt sqrt(0.5 d/dt strain_ij d/dt strain_ij).
  static constexpr double YieldingNodeStrainXY = TimeStep * ShearStress / (4 * Mu * RelaxationTime);
  static inline const double YieldingNodeEta = YieldingNodeStrainXY * std::sqrt(0.5);

  // Runs the kernel once on the zero DOFs and plastic strains; node `i` yields iff `yields(i)`.
  template <typename YieldsT>
  std::size_t run(const YieldsT& yields) {
    std::vector<model::Plasticity> parameters(NumNodes);
    for (std::size_t node = 0; node < NumNodes; ++node) {
      parameters[node].bulkFriction = 0;
      parameters[node].plastCo = yields(node) ? ShearStress / 2 : ShearStress * 2;
      parameters[node].sXY = ShearStress;
    }
    std::array<const model::Plasticity*, Cfg::NumSimulations> perSimulation{};
    perSimulation.fill(parameters.data());

    model::ElasticMaterial material;
    material.mu = Mu;
    material.lambda = Mu;
    const model::PlasticityData<Cfg> plasticityData(perSimulation, &material, true);

    dofs_.fill(0);
    pstrain_.fill(0);
    const GlobalData<Cfg> global = seissol::Pool<Cfg>::host();
    return Plasticity::computePlasticity(
        static_cast<real>(Plasticity::computeRelaxTime(RelaxationTime, TimeStep)),
        static_cast<real>(TimeStep),
        static_cast<real>(RelaxationTime),
        &global,
        &plasticityData,
        dofs_.data(),
        pstrain_.data());
  }

  [[nodiscard]] bool dofsUnchanged() const {
    return std::all_of(dofs_.begin(), dofs_.end(), [](real value) { return value == 0; });
  }

  [[nodiscard]] bool plasticStrainUnchanged() const {
    return std::all_of(pstrain_.begin(), pstrain_.end(), [](real value) { return value == 0; });
  }

  // pstrain_ holds the plastic strain (in the layout of QStressNodal), followed by eta
  [[nodiscard]] double
      plasticStrain(std::size_t sim, std::size_t node, std::size_t component) const {
    auto view = init::QStressNodal<Cfg>::view::create(pstrain_.data());
    return multisim::simtensor<Cfg>(view, static_cast<int>(sim))(node, component);
  }

  [[nodiscard]] double eta(std::size_t sim, std::size_t node) const {
    auto view =
        init::QEtaNodal<Cfg>::view::create(pstrain_.data() + tensor::QStressNodal<Cfg>::size());
    return multisim::simtensor<Cfg>(view, static_cast<int>(sim))(node);
  }

  private:
  // the kernel reads and writes the six stress quantities of the DOFs, via the tensor QStress
  static constexpr std::size_t DofsSize =
      std::max(tensor::Q<Cfg>::size(), tensor::QStress<Cfg>::size());

  alignas(Alignment) std::array<real, DofsSize> dofs_{};
  alignas(Alignment)
      std::array<real,
                 tensor::QStressNodal<Cfg>::size() + tensor::QEtaNodal<Cfg>::size()> pstrain_{};
};

TEST_CASE_TEMPLATE("Plasticity computePlasticity corrects every cell with a yielding node" *
                       doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using Cell = ShearedPlasticityCell<Cfg>;
  constexpr auto NumNodes = Cell::NumNodes;
  Cell cell;

  SUBCASE("No node yields: the cell stays unchanged") {
    CHECK(cell.run([](std::size_t /*node*/) { return false; }) == 0);
    CHECK(cell.dofsUnchanged());
    CHECK(cell.plasticStrainUnchanged());
  }

  SUBCASE("Exactly one node yields, at every position in turn") {
    // The yield check reduces over all nodes in SIMD chunks; a single yielding node has to trigger
    // the correction wherever it is, from the first to the last real node.
    for (std::size_t yieldingNode = 0; yieldingNode < NumNodes; ++yieldingNode) {
      CAPTURE(yieldingNode);
      CHECK(cell.run([&](std::size_t node) { return node == yieldingNode; }) == 1);
      CHECK_FALSE(cell.dofsUnchanged());
      for (std::size_t sim = 0; sim < Cfg::NumSimulations; ++sim) {
        for (std::size_t node = 0; node < NumNodes; ++node) {
          if (node == yieldingNode) {
            CHECK(cell.plasticStrain(sim, node, Cell::ComponentXY) / Cell::YieldingNodeStrainXY ==
                  doctest::Approx(1.0));
            CHECK(cell.eta(sim, node) / Cell::YieldingNodeEta == doctest::Approx(1.0));
          } else {
            CHECK(cell.plasticStrain(sim, node, Cell::ComponentXY) == 0);
            CHECK(cell.eta(sim, node) == 0);
          }
        }
      }
    }
  }

  SUBCASE("Every node yields") {
    CHECK(cell.run([](std::size_t /*node*/) { return true; }) == 1);
    CHECK_FALSE(cell.dofsUnchanged());
    for (std::size_t sim = 0; sim < Cfg::NumSimulations; ++sim) {
      for (std::size_t node = 0; node < NumNodes; ++node) {
        CAPTURE(node);
        CHECK(cell.eta(sim, node) / Cell::YieldingNodeEta == doctest::Approx(1.0));
      }
    }
  }
}

} // namespace seissol::unit_test
