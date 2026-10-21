// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Expr/Backend.h"
#include "Expr/SderivFrontend.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/GhostStateBatch.h"
#include "Model/ScriptParameters.h"
#include "Model/StateNames.h"
#include "Physics/GhostState.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace seissol::unit_test::ghost_state_batch {

using Cfg = Config0;
using RealT = Real<Cfg>;

/// The kernel of a batch of three faces (of the cells 2, 0 and 3 of a layer of four), against the
/// integral of the host evaluated face by face: the same values where the program gives a
/// component, and the tensor untouched where it does not.
template <bool Aligned>
void checkAgainstHost(const std::string& source, bool readsInner) {
  using Batch = kernels::GhostStateBatch<Cfg, Aligned>;
  constexpr std::size_t Faces = 3;
  constexpr std::size_t Nodes = Batch::NodeCount;
  constexpr std::size_t Size = Batch::TargetSize;
  constexpr std::size_t TimePoints = Cfg::ConvergenceOrder;
  constexpr std::size_t ParameterCount = Batch::ParameterCount;

  auto components = model::stateNames<model::MaterialOf<Cfg>>();
  components.resize(Batch::Columns);
  const auto parameters = model::scriptParameterNames<model::MaterialOf<Cfg>>();
  expr::SderivOptions options;
  if (readsInner) {
    options.inputs.insert(components.begin(), components.end());
  }
  const auto program = expr::compileSderivModule(source, options);
  const physics::GhostStateIntegral integral(program,
                                             "test",
                                             components,
                                             readsInner ? components : std::vector<std::string>{},
                                             parameters,
                                             physics::GhostStateIntegral::Missing::Zero,
                                             TimePoints);
  REQUIRE(integral.innerCount() == (readsInner ? components.size() : 0));

  const std::vector<std::uint32_t> cells = {2, 0, 3};
  std::vector<RealT> nodes(Faces * Nodes * 3);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    nodes[i] = static_cast<RealT>(0.125 * static_cast<double>(i % 17) - 1.0);
  }
  std::vector<RealT> materials(4 * ParameterCount);
  for (std::size_t i = 0; i < materials.size(); ++i) {
    materials[i] = static_cast<RealT>(1.0 + 0.5 * static_cast<double>(i));
  }
  std::vector<RealT> inner(Faces * TimePoints * Size);
  for (std::size_t i = 0; i < inner.size(); ++i) {
    inner[i] = static_cast<RealT>(0.01 * static_cast<double>(i % 101) - 0.3);
  }
  std::vector<RealT> ghost(Faces * Size, -7.0);

  typename Batch::Data data;
  data.faces = Faces;
  data.nodes = nodes.data();
  data.parameters = materials.data();
  data.cells = cells.data();
  data.inner = readsInner ? inner.data() : nullptr;
  data.innerStride = TimePoints * Size;
  data.ghost = ghost.data();
  expr::BackendOptions backend;
  backend.preferred = expr::BackendKind::Interpreter;
  Batch batch(integral, data, backend);

  const double start = 0.5;
  const double width = 0.25;
  batch.evaluate(start, width, 0, 1.0, nullptr);

  std::vector<bool> given(components.size(), false);
  for (std::size_t j = 0; j < integral.program().outputs().size(); ++j) {
    given[integral.outputComponent(j)] = true;
  }
  for (std::size_t e = 0; e < Faces; ++e) {
    std::vector<std::array<double, 3>> points(Nodes);
    for (std::size_t p = 0; p < Nodes; ++p) {
      for (std::size_t d = 0; d < 3; ++d) {
        points[p][d] = nodes[(e * Nodes + p) * 3 + d];
      }
    }
    std::vector<double> material(ParameterCount);
    for (std::size_t k = 0; k < ParameterCount; ++k) {
      material[k] = materials[cells[e] * ParameterCount + k];
    }
    std::vector<double> innerOfFace(TimePoints * integral.innerCount() * Nodes);
    for (std::size_t t = 0; t < TimePoints; ++t) {
      const auto view = Batch::TargetInit::view::create(const_cast<RealT*>(inner.data()) +
                                                        (e * TimePoints + t) * Size);
      const auto slice = multisim::simtensor<Cfg>(view, 0);
      for (std::size_t k = 0; k < integral.innerCount(); ++k) {
        for (std::size_t p = 0; p < Nodes; ++p) {
          innerOfFace[(t * integral.innerCount() + k) * Nodes + p] = slice(p, k);
        }
      }
    }
    std::vector<double> expected(components.size() * Nodes);
    integral.evaluate(start,
                      width,
                      1,
                      points.data(),
                      Nodes,
                      material.data(),
                      innerOfFace.data(),
                      expected.data());

    const auto view = Batch::TargetInit::view::create(ghost.data() + e * Size);
    const auto slice = multisim::simtensor<Cfg>(view, 0);
    for (std::size_t k = 0; k < components.size(); ++k) {
      for (std::size_t p = 0; p < Nodes; ++p) {
        CAPTURE(e);
        CAPTURE(k);
        CAPTURE(p);
        if (given[k]) {
          CHECK(slice(p, k) ==
                doctest::Approx(static_cast<RealT>(expected[k * Nodes + p])).epsilon(1e-6));
        } else {
          CHECK(slice(p, k) == -7.0);
        }
      }
    }
  }
}

TEST_CASE("GhostStateBatch: the faces of a batch in one kernel, as the host does them face by "
          "face") {
  const auto components = model::stateNames<model::MaterialOf<Cfg>>();
  const auto& first = components.front();
  const auto& last = components[kernels::GhostStateBatch<Cfg, false>::Columns - 1];

  SUBCASE("the analytical boundary: the nodes, the time, sim and the material") {
    checkAgainstHost<false>("out def " + last + " = x * t + rho + sim\nout def " + first +
                                " = y * y - z * t * t\n",
                            false);
  }
  SUBCASE("the nonlinear Dirichlet boundary: the state inside at the time points") {
    checkAgainstHost<false>("out def " + first + " = -2.0 * " + first + " + x * t\nout def " +
                                last + " = " + last + " * " + first + "\n",
                            true);
  }
  SUBCASE("the same in the face-aligned basis") {
    const auto& lastAligned = components[kernels::GhostStateBatch<Cfg, true>::Columns - 1];
    checkAgainstHost<true>("out def " + first + " = -" + first + " * rho\nout def " + lastAligned +
                               " = " + lastAligned + " + z\n",
                           true);
  }
}

} // namespace seissol::unit_test::ghost_state_batch
