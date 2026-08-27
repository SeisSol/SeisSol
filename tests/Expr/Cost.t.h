// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Cost.h"
#include "Expr/Lower.h"
#include "Expr/SderivFrontend.h"
namespace seissol::expr::test {
TEST_SUITE("ExprCost") {
  TEST_CASE("costs are counted per point, exactly") {
    // Exact rather than estimated: after lowering there is no branching left,
    // so the instruction list IS the execution.
    const Program simple = compileSderivModule("out def u = x + 1.0\n");
    const auto one = cost(lower(simple), simple.computeType());
    CHECK(one.run.additive == 1);
    CHECK(one.run.transcendentals == 0);
    CHECK(one.run.stores == 1);

    const Program wave =
        compileSderivModule("def phi = x + y\nout def u = sin(phi)\nout def v = cos(phi)\n");
    const auto two = cost(lower(wave), wave.computeType());
    CHECK(two.run.transcendentals == 2);
    // A transcendental is weighted, so it dominates a handful of adds -- which
    // is the whole reason the weighted figure exists.
    CHECK(two.run.weighted() > two.run.operations());
  }

  TEST_CASE("a Select costs both arms") {
    // Interp.h evaluates both, so the count has to as well.
    const Program program = compileSderivModule("out def u = select(lt(x,y), sqrt(x), sqrt(y))\n");
    const auto costs = cost(lower(program), program.computeType());
    CHECK(costs.run.transcendentals == 2);
  }

  TEST_CASE("hoisting moves cost out of the run stage") {
    const char* source = "def vp = sqrt((lambda + 2.0*mu)/rho)\n"
                         "out def u = rho*vp*sin(x - vp*t)\n";
    const Program program = compileSderivModule(source);

    const auto plain = cost(lower(program), program.computeType());
    LowerOptions options;
    options.invariantInputs = {"rho", "mu", "lambda"};
    options.hoistThreshold = 1;
    const auto hoisted = cost(lower(program, options), program.computeType());

    CHECK(hoisted.precompute.operations() > 0);
    CHECK(hoisted.run.weighted() < plain.run.weighted());
  }

  TEST_CASE("a lookup is counted, not costed") {
    // Its cost is a stencil gather and a tensor-product reduction, not an
    // instruction; a number that pretended otherwise would be worse than none.
    const Program program =
        compileSderivModule("grid m = \"asagi\", \"model.nc\", \"linear\", \"rho\"\n"
                            "out def u = m_rho(x, y, z)\n");
    const auto costs = cost(lower(program), program.computeType());
    CHECK(costs.run.lookups == 1);
    CHECK(costs.run.weighted() == 0);
  }
}
} // namespace seissol::expr::test
