// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Codegen.h"
#include "Expr/Lower.h"
#include "Expr/RtcCpu.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::Direction;

/// Deliberately not "nice" values. Zero, both signed zeros, the subnormal end
/// and the overflow end are where an emitted expression and an interpreted one
/// stop agreeing -- a ladder of 0.5, 1.0, 2.0 would pass for almost any bug.
const std::vector<double> Ladder = {
    -1e6, -1e3, -7.5, -1.0, -1e-3, -0.0, 0.0, 1e-3, 0.5, 1.0, 2.0, 3.7, 1e3, 1e6};

/// Evaluate `source` on both backends and return true when every output is
/// bitwise equal. `compiled` reports whether the RTC backend was actually used;
/// on a machine with no usable C++ compiler it falls back, and a test that
/// silently passed on the interpreter twice would be worthless.
bool agreesBitwise(const std::string& source, bool& compiled) {
  const Program program = compileSderivModule(source);
  const std::size_t numPoints = Ladder.size();

  std::vector<double> x(numPoints);
  std::vector<double> y(numPoints);
  std::vector<double> z(numPoints);
  std::vector<double> t(numPoints);
  for (std::size_t i = 0; i < numPoints; ++i) {
    x[i] = Ladder[i];
    y[i] = Ladder[(i + 3) % numPoints];
    z[i] = Ladder[(i + 7) % numPoints];
    t[i] = Ladder[(i + 5) % numPoints];
  }

  const std::size_t outputs = program.outputs().size();
  std::vector<double> interpreted(outputs * numPoints, -1.0);
  std::vector<double> emitted(outputs * numPoints, -1.0);

  const auto evaluate = [&](std::vector<double>& out, BackendKind preferred) {
    DataTable table(numPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindViewConst<double>("z", Direction::In, z.data());
    table.bindViewConst<double>("t", Direction::In, t.data());
    for (std::size_t i = 0; i < outputs; ++i) {
      table.bindView<double>(program.outputs()[i].name, Direction::Out, &out[i * numPoints]);
    }
    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    BackendOptions options;
    options.preferred = preferred;
    const auto kernel = makeKernel(program, binding, store, options);
    kernel->precompute(table);
    kernel->run(table);
    return kernel->kind();
  };

  evaluate(interpreted, BackendKind::Interpreter);
  compiled = evaluate(emitted, BackendKind::RtcCpu) == BackendKind::RtcCpu;

  // memcmp rather than ==, so NaN counts as equal to NaN and -0.0 does NOT
  // count as equal to 0.0. Both distinctions matter: a sign-of-zero difference
  // is a real divergence that == would hide.
  return std::memcmp(interpreted.data(), emitted.data(), interpreted.size() * sizeof(double)) == 0;
}

} // namespace

TEST_SUITE("ExprRtcCpu") {

  // These need no compiler and so run everywhere.

  TEST_CASE("the emitted source turns transients into locals, not slots") {
    const Program program = compileSderivModule("out def u = sqrt(x*x + y*y)\n");
    const std::string source = emitCpuSource(program, lower(program));

    CHECK(source.find("seissol_expr_run") != std::string::npos);
    CHECK(source.find("for (unsigned long l = 0") != std::string::npos);
    // The whole reason this is faster than the interpreter: the values live in
    // locals the vectoriser can see through, not in a scratch array.
    CHECK(source.find("double s0") != std::string::npos);
    CHECK(source.find("std::sqrt") != std::string::npos);
  }

  TEST_CASE("substitution is token aware") {
    // std::exp contains an `x`. A plain string replace of the operand
    // placeholder would corrupt the function name into something that still
    // compiles for some operand names and computes nonsense.
    const Program program = compileSderivModule("out def u = exp(x)\n");
    const std::string source = emitCpuSource(program, lower(program));
    CHECK(source.find("std::exp(s") != std::string::npos);
  }

  TEST_CASE("constants keep every bit") {
    // %.17g, not std::to_string: the latter truncates at six decimals, which
    // would quietly emit a different program than the one that was lowered.
    const Program program = compileSderivModule("out def u = x * 0.1234567890123456789\n");
    const std::string source = emitCpuSource(program, lower(program));
    CHECK(source.find("0.12345678901234568") != std::string::npos);
  }

  TEST_CASE("the first call of a compiled kernel is compared with the interpreter") {
    // x * y + z, at x = 1 + 2^-27, y = 1 - 2^-27, z = -1: rounded on its own, the product is 1 and
    // the sum 0; contracted into a fused multiply-add, the sum is -2^-54
    const Program program = compileSderivModule("out def u = x * y + z\n");
    std::vector<double> x = {1.0 + std::ldexp(1.0, -27), 2.0};
    std::vector<double> y = {1.0 - std::ldexp(1.0, -27), 3.0};
    std::vector<double> z = {-1.0, 1.0};
    std::vector<double> u(2, 0.0);
    reader::scripting::DataTable table(2);
    table.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
    table.bindViewConst<double>("y", reader::scripting::Direction::In, y.data());
    table.bindViewConst<double>("z", reader::scripting::Direction::In, z.data());
    table.bindView<double>("u", reader::scripting::Direction::Out, u.data());
    auto binding = Binding::bind(program, table);
    df::GridStore store;
    BackendOptions options;
    options.quiet = true;
    options.preferred = BackendKind::RtcCpu;

    // with the flags of SeisSol, the same bits, which the check finds at its first call
    {
      const auto kernel = makeKernel(program, binding, store, options);
      if (kernel->kind() != BackendKind::RtcCpu) {
        WARN_MESSAGE(false, "no usable C++ compiler; the compiled CPU backend was not run");
        return;
      }
      kernel->precompute(table);
      kernel->run(table);
      const auto comparison = compareWithInterpreter(program, binding, store, &table, nullptr);
      REQUIRE(comparison.has_value());
      CHECK(comparison->points.size() == 2);
      CHECK(comparison->deviating == 0);
      CHECK(u[0] == 0.0);
      REQUIRE(comparison->inputs.size() == 6);
      CHECK(comparison->interpreter[1] == 7.0);
    }

    // a kernel compiled with contraction deviates, and the check says where
    setenv("SEISSOL_EXPR_CXXFLAGS", "-O2 -ffp-contract=fast -march=native", 1);
    setenv("SEISSOL_EXPR_CHECK", "warn", 1);
    const auto contracted = makeKernel(program, binding, store, options);
    contracted->precompute(table);
    contracted->run(table);
    unsetenv("SEISSOL_EXPR_CXXFLAGS");
    unsetenv("SEISSOL_EXPR_CHECK");
    const auto comparison = compareWithInterpreter(program, binding, store, &table, nullptr);
    REQUIRE(comparison.has_value());
    if (u[0] == 0.0) {
      WARN_MESSAGE(false, "the compiler did not contract; the check had nothing to find");
      return;
    }
    CHECK(comparison->deviating == 1);
    CHECK(comparison->worst.find("`u` at point 0") != std::string::npos);

    // against the perturbed reference, as a GPU is checked, the contraction is within the spread
    // of products moved by an ulp; a value beyond it is not
    const auto tolerant = compareWithInterpreter(program, binding, store, &table, nullptr, true);
    REQUIRE(tolerant.has_value());
    CHECK(tolerant->perturbed);
    CHECK(tolerant->deviating == 0);
    CHECK(tolerant->tolerance[0] >= std::ldexp(1.0, -52));
    CHECK(tolerant->worstRatio <= 1.0);
    CHECK(tolerant->unstable.empty());
    u[0] = 1e-10;
    const auto beyond = compareWithInterpreter(program, binding, store, &table, nullptr, true);
    REQUIRE(beyond.has_value());
    CHECK(beyond->deviating == 1);
    CHECK(beyond->worstRatio > 1.0);
  }

  TEST_CASE("a program with a grid lookup compiles and samples as the interpreter does") {
    // written out (expandLookups()): the compiled kernel reads the grid through its descriptor, in
    // both shapes, at points inside the grid and far beyond it
    const Program program = compileSderivModule(
        "grid m = \"asagi\", \"ramp.nc\", \"linear\", \"a\", \"b\"\n"
        "out def u = m_b(x, y, z) * 0.5 + m_a(z, x, y)\nout def v = m_a(x, y, z)\n");
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    for (std::size_t i = 0; i < Ladder.size(); ++i) {
      x.push_back(Ladder[i]);
      y.push_back(3.25 + 0.37 * static_cast<double>(i));
      z.push_back(7.0 - 0.61 * static_cast<double>(i));
    }
    const std::size_t numPoints = x.size();
    const auto evaluate = [&](BackendKind backend, KernelShape shape, std::vector<double>& out) {
      out.assign(2 * numPoints, -1.0);
      DataTable table(numPoints);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindViewConst<double>("y", Direction::In, y.data());
      table.bindViewConst<double>("z", Direction::In, z.data());
      table.bindView<double>("u", Direction::Out, out.data());
      table.bindView<double>("v", Direction::Out, out.data() + numPoints);
      Binding binding = Binding::bind(program, table);
      df::GridStore store;
      const std::size_t id = store.intern(program.grids()[0]);
      store.injectForTesting(id, std::make_unique<LinearRampGrid>());
      BackendOptions options;
      options.preferred = backend;
      options.shape = shape;
      options.quiet = true;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);
      kernel->run(table);
      return kernel->kind();
    };
    std::vector<double> interpreted;
    evaluate(BackendKind::Interpreter, KernelShape::Auto, interpreted);
    for (std::size_t p = 0; p < numPoints; ++p) {
      // a ramp, which linear interpolation reproduces where the point is inside
      const auto inside = [](double c) { return c >= 0.0 && c <= 7.0; };
      if (inside(x[p]) && inside(y[p]) && inside(z[p])) {
        CHECK(interpreted[numPoints + p] ==
              doctest::Approx(LinearRampGrid::value_(x[p], y[p], z[p])).epsilon(1e-12));
      }
    }
    for (const auto shape : {KernelShape::Point, KernelShape::Tile}) {
      CAPTURE(static_cast<int>(shape));
      std::vector<double> compiled;
      if (evaluate(BackendKind::RtcCpu, shape, compiled) != BackendKind::RtcCpu) {
        WARN_MESSAGE(false, "no usable C++ compiler; the RTC backend fell back");
        return;
      }
      CHECK(std::memcmp(interpreted.data(), compiled.data(), interpreted.size() * sizeof(double)) ==
            0);
    }
  }

  // These need a working C++ compiler. Where there is none the backend falls
  // back and the case reports that rather than passing vacuously.

  TEST_CASE("the compiled kernel is bitwise identical to the interpreter") {
    bool compiled = false;
    const std::vector<std::pair<const char*, const char*>> programs = {
        {"arithmetic and powers", "out def u = (x*y - z)/(t*t + 1.0) + x**3.0\n"},
        {"roots, exp and logs",
         "def a = abs(x) + 1.0\nout def u = sqrt(a) + exp(0.0-a) + log(a) + log2(a) + log10(a)\n"},
        {"trigonometry", "out def u = sin(x)+cos(y)+tan(z)+atan(t)+atan2(x,y)\n"},
        {"hyperbolics", "out def u = sinh(x/1e3)+cosh(y/1e3)+tanh(z)+erf(t)\n"},
        {"rounding and sign", "out def u = floor(x)+ceil(y)+round(z)+sign(t)+abs(x)\n"},
        {"min, max, mod", "out def u = min(x,y)+max(z,t)+mod(x, 3.0)\n"},
        {"comparisons and select",
         "out def u = select(land(lt(x,y), ge(z,t)), x, y) + select(lnot(eq(x,y)), 1.0, 2.0)\n"},
        {"several outputs", "def a = x*x+y*y\nout def u = a\nout def v = a-z\nout def w = a*t\n"},
        {"a shared subexpression", "def a = sqrt(abs(x))\nout def u = a+a*a\nout def v = a*a*a\n"},
    };

    for (const auto& program : programs) {
      // (a structured binding cannot be captured where OpenMP is enabled)
      const char* const label = program.first;
      const char* const source = program.second;
      CAPTURE(label);
      const bool same = agreesBitwise(source, compiled);
      if (!compiled) {
        WARN_MESSAGE(compiled, "no usable C++ compiler; the RTC backend fell back");
        break;
      }
      CHECK(same);
    }
  }

  TEST_CASE("run(KernelArgs) and run(table) agree on the compiled backend too") {
    const Program program = compileSderivModule(
        "def a = x*x + y*y\nout def u = a + z\nout def v = a - z\nout def w = a*z\n");
    constexpr std::size_t Nodes = 9;

    std::vector<double> face(3 * Nodes);
    std::vector<double> out(3 * Nodes, -1.0);
    for (std::size_t i = 0; i < face.size(); ++i) {
      face[i] = 0.1 * static_cast<double>(i + 1);
    }

    DataTable table(Nodes);
    table.bindViewConst<double>("x", Direction::In, face.data(), 3, 0);
    table.bindViewConst<double>("y", Direction::In, face.data(), 3, 1);
    table.bindViewConst<double>("z", Direction::In, face.data(), 3, 2);
    table.bindView<double>("u", Direction::Out, out.data(), 3, 0);
    table.bindView<double>("v", Direction::Out, out.data(), 3, 1);
    table.bindView<double>("w", Direction::Out, out.data(), 3, 2);

    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    BackendOptions options;
    options.preferred = BackendKind::RtcCpu;
    const auto kernel = makeKernel(program, binding, store, options);
    kernel->precompute(table);
    kernel->run(table);
    const std::vector<double> reference = out;

    KernelArgs args{};
    const void* inputs[3] = {face.data(), face.data(), face.data()};
    void* outputs[3] = {out.data(), out.data(), out.data()};
    args.inputs = inputs;
    args.inputCount = 3;
    args.outputs = outputs;
    args.outputCount = 3;

    // One point per call: the shape EasiBoundary::query needs, and it must not
    // be a second numerical path.
    std::fill(out.begin(), out.end(), -1.0);
    for (std::size_t n = 0; n < Nodes; ++n) {
      args.first = n;
      args.count = 1;
      kernel->run(args);
    }
    for (std::size_t i = 0; i < out.size(); ++i) {
      CHECK(out[i] == reference[i]);
    }
  }

  TEST_CASE("two kernels for one program share the compiled artifact") {
    const Program program = compileSderivModule("out def u = x * 17.0 + 1.0\n");
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> out(2, 0.0);

    DataTable table(2);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindView<double>("u", Direction::Out, out.data());

    df::GridStore store;
    BackendOptions options;
    options.preferred = BackendKind::RtcCpu;

    Binding first = Binding::bind(program, table);
    const auto a = makeKernel(program, first, store, options);
    const std::size_t afterFirst = rtcCpuCacheSize();

    Binding second = Binding::bind(program, table);
    const auto b = makeKernel(program, second, store, options);
    // Same program, same options, same flags: compiling it twice would be pure
    // waste, and the fingerprint exists to notice.
    CHECK(rtcCpuCacheSize() == afterFirst);
    CHECK(a->kind() == b->kind());
  }

  TEST_CASE("a deep program takes the shape of a tile, with the bits of the interpreter") {
    // Newton iterations, unrolled: a chain of dependent operations, one of them through the maths
    // library, which keeps the loop over the points from being vectorized
    std::string source = "def v0 = 1.0 + x * x\ndef w0 = 0.5\n";
    constexpr int Steps = 12;
    for (int i = 1; i <= Steps; ++i) {
      const auto previous = std::to_string(i - 1);
      const auto current = std::to_string(i);
      source += "def v" + current + " = v" + previous + " - (v" + previous + " * v" + previous +
                " * v" + previous + " - 2.0 - y) / (3.0 * v" + previous + " * v" + previous +
                " + 1.0)\n";
      source += "def w" + current + " = w" + previous + " - (w" + previous + " + sinh(w" +
                previous + ") - x) / (1.0 + cosh(w" + previous + "))\n";
    }
    source +=
        "out def u = v" + std::to_string(Steps) + "\nout def w = w" + std::to_string(Steps) + "\n";
    const Program program = compileSderivModule(source);
    const auto lowered = lower(program);
    REQUIRE(codegen::criticalPath(lowered.run(), lowered.operands()) >= TileDepth);
    CHECK(emitCpuSource(program, lowered).find("#pragma omp simd") != std::string::npos);
    const Program shallow = compileSderivModule("out def u = x * y + 1.0\n");
    CHECK(emitCpuSource(shallow, lower(shallow)).find("#pragma omp simd") == std::string::npos);

    // more points than a block, and not a multiple of one
    constexpr std::size_t NumPoints = 1000;
    std::vector<double> x(NumPoints);
    std::vector<double> y(NumPoints);
    for (std::size_t i = 0; i < NumPoints; ++i) {
      x[i] = Ladder[i % Ladder.size()] * 1e-3 + 0.25 * static_cast<double>(i % 7);
      y[i] = 0.1 * static_cast<double>(i % 13) - 0.5;
    }
    const auto evaluate = [&](BackendKind preferred, KernelShape shape, BackendKind& kind) {
      std::vector<double> out(2 * NumPoints, -1.0);
      DataTable table(NumPoints);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindViewConst<double>("y", Direction::In, y.data());
      table.bindView<double>("u", Direction::Out, out.data());
      table.bindView<double>("w", Direction::Out, out.data() + NumPoints);
      Binding binding = Binding::bind(program, table);
      df::GridStore store;
      BackendOptions options;
      options.preferred = preferred;
      options.shape = shape;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);
      kernel->run(table);
      kind = kernel->kind();
      return out;
    };
    BackendKind kind = BackendKind::Interpreter;
    const auto interpreted = evaluate(BackendKind::Interpreter, KernelShape::Auto, kind);
    const auto tiled = evaluate(BackendKind::RtcCpu, KernelShape::Tile, kind);
    if (kind != BackendKind::RtcCpu) {
      WARN_MESSAGE(false, "no usable host compiler; the shape of a tile was not compiled");
      return;
    }
    const auto pointwise = evaluate(BackendKind::RtcCpu, KernelShape::Point, kind);
    REQUIRE(kind == BackendKind::RtcCpu);
    for (std::size_t i = 0; i < interpreted.size(); ++i) {
      CAPTURE(i);
      CHECK(std::memcmp(&tiled[i], &interpreted[i], sizeof(double)) == 0);
      CHECK(std::memcmp(&pointwise[i], &interpreted[i], sizeof(double)) == 0);
    }
  }

} // TEST_SUITE

} // namespace seissol::expr::test
