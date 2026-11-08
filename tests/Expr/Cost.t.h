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
    const auto one = cost(simple, lower(simple), simple.computeType());
    CHECK(one.run.additive == 1);
    CHECK(one.run.transcendentals == 0);
    CHECK(one.run.stores == 1);

    const Program wave =
        compileSderivModule("def phi = x + y\nout def u = sin(phi)\nout def v = cos(phi)\n");
    const auto two = cost(wave, lower(wave), wave.computeType());
    CHECK(two.run.transcendentals == 2);
    // A transcendental is weighted, so it dominates a handful of adds -- which
    // is the whole reason the weighted figure exists.
    CHECK(two.run.weighted() > two.run.operations());
  }

  TEST_CASE("the floating-point operations count every operation once, a multiply-add twice") {
    // as SeisSol counts the operations of its generated kernels: x * y, + and sin; the loads,
    // the stores and the literals are none
    const Program program = compileSderivModule("out def u = x * y + sin(x)\n");
    CHECK(flopsPerPoint(program) == 3);
    // x is loaded in binary64 (no operation), * 2.0, the conversion back, + y
    const Program converted = compileSderivModule("out def u = real(f64(x) * 2.0) + y\n");
    CHECK(flopsPerPoint(converted) == 3);
    StageCost contraction;
    contraction.fma = 20;
    contraction.additive = 1;
    CHECK(contraction.flops() == 41);
  }

  TEST_CASE("a Select costs both arms") {
    // Interp.h evaluates both, so the count has to as well.
    const Program program = compileSderivModule("out def u = select(lt(x,y), sqrt(x), sqrt(y))\n");
    const auto costs = cost(program, lower(program), program.computeType());
    CHECK(costs.run.transcendentals == 2);
  }

  TEST_CASE("hoisting moves cost out of the run stage") {
    const char* source = "def vp = sqrt((lambda + 2.0*mu)/rho)\n"
                         "out def u = rho*vp*sin(x - vp*t)\n";
    const Program program = compileSderivModule(source);

    const auto plain = cost(program, lower(program), program.computeType());
    LowerOptions options;
    options.invariantInputs = {"rho", "mu", "lambda"};
    options.hoistThreshold = 1;
    const auto hoisted = cost(program, lower(program, options), program.computeType());

    CHECK(hoisted.precompute.operations() > 0);
    CHECK(hoisted.run.weighted() < plain.run.weighted());
  }

  TEST_CASE("a lookup is counted as what it computes") {
    // written out: the weights and the tensor product as operations, and the elements of the
    // stencil -- 2^3 for a linear one in three dimensions -- as the random accesses they are
    const Program program =
        compileSderivModule("grid m = \"asagi\", \"model.nc\", \"linear\", \"rho\"\n"
                            "out def u = m_rho(x, y, z)\n");
    const auto costs = cost(program, lower(program), program.computeType());
    CHECK(costs.run.gridLoads == 8);
    // a sum of eight products, and per axis two weights from the coordinate
    CHECK(costs.run.multiplicative >= 8 + 3);
    CHECK(costs.run.additive >= 8 + 3 * 3);
  }
}
} // namespace seissol::expr::test
