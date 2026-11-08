// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcCpu.h"
#include "Expr/SderivFrontend.h"
#include "Expr/Types.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::Direction;

/// The message of what compiling `source` throws, or empty.
std::string errorOf(const std::string& source) {
  try {
    static_cast<void>(compileSderivModule(source));
  } catch (const std::exception& error) {
    return error.what();
  }
  return {};
}

/// The outputs of `program` (in its compute type `type`) at the points `x` and `t` (double
/// columns), into double columns, after `calls` runs, on the backend `preferred`; whether that
/// backend ran.
bool evaluate(const std::string& source,
              ComputeType type,
              const std::vector<double>& x,
              const std::vector<double>& t,
              std::size_t calls,
              BackendKind preferred,
              std::vector<double>& out) {
  Program program = compileSderivModule(source);
  program.setComputeType(type);
  const std::size_t points = x.size();
  out.assign(program.outputs().size() * points, -1.0);
  DataTable table(points);
  table.bindViewConst<double>("x", Direction::In, x.data());
  table.bindViewConst<double>("t", Direction::In, t.data());
  for (std::size_t i = 0; i < program.outputs().size(); ++i) {
    table.bindView<double>(program.outputs()[i].name, Direction::Out, &out[i * points]);
  }
  Binding binding = Binding::bind(program, table);
  df::GridStore store;
  BackendOptions options;
  options.preferred = preferred;
  options.quiet = true;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  for (std::size_t call = 0; call < calls; ++call) {
    kernel->run(table);
  }
  return kernel->kind() == preferred;
}

} // namespace

TEST_SUITE("ExprTypes") {

  TEST_CASE("conversions and typed states parse, and type the program") {
    const Program program = compileSderivModule(
        "state f64 s = 0.0\nout def s = s + f64(x) * 2.0\nout def u = real(s) * x\n");
    REQUIRE(program.state().size() == 1);
    CHECK(program.state()[0].type == ValueType::F64);

    Program single = program;
    single.setComputeType(ComputeType::F32);
    CHECK(dataType(single) == ComputeType::F64);
    const LoweredProgram lowered = lower(single);
    CHECK(lowered.dataType() == ComputeType::F64);
    CHECK(lowered.computeType() == ComputeType::F32);
    CHECK(lowered.mixesTypes());
    bool narrow = false;
    bool wide = false;
    for (const Instruction& inst : lowered.run().code) {
      narrow = narrow || inst.type == ComputeType::F32;
      wide = wide || inst.type == ComputeType::F64;
    }
    CHECK(narrow);
    CHECK(wide);
  }

  TEST_CASE("a program without conversions moves its data in its compute type") {
    Program program = compileSderivModule("state s = 0.0\nout def s = s + x\n");
    program.setComputeType(ComputeType::F32);
    CHECK(dataType(program) == ComputeType::F32);
    const LoweredProgram lowered = lower(program);
    CHECK_FALSE(lowered.mixesTypes());
    // the same source as before there were types
    const std::string source = emitCpuSource(program, lowered);
    CHECK(source.find("(double)") == std::string::npos);
    CHECK(source.find("double s") == std::string::npos);
  }

  TEST_CASE("typing a program without conversions keeps the order of its values") {
    // so that the lowering emits the same instructions as before there were types: a device kernel
    // that contracts into fused multiply-adds would round otherwise where its source changed
    const Program program = compileSderivModule("def a = x * x + t\ndef b = sin(a) * 2.0 + a\nout "
                                                "def u = b * b - a\nout def v = a + 1.0\n");
    const Program result = typed(program);
    const auto sequence = [](const Program& p) {
      std::vector<char> reachable(p.arena().size(), 0);
      std::vector<NodeId> stack(p.roots().begin(), p.roots().end());
      while (!stack.empty()) {
        const NodeId id = stack.back();
        stack.pop_back();
        if (reachable[static_cast<std::size_t>(id)] == 0) {
          reachable[static_cast<std::size_t>(id)] = 1;
          for (const NodeId kid : p.arena().children(id)) {
            stack.push_back(kid);
          }
        }
      }
      std::vector<std::string> nodes;
      for (std::size_t id = 0; id < reachable.size(); ++id) {
        if (reachable[id] != 0) {
          const Node& node = p.arena()[static_cast<NodeId>(id)];
          nodes.push_back(std::string(name(node.kind)) + ":" +
                          (node.kind == Kind::PW ? name(node.fn) : "") + ":" +
                          std::to_string(node.value));
        }
      }
      return nodes;
    };
    CHECK(sequence(program) == sequence(result));
    CHECK(emitCpuSource(program, lower(program)) == emitCpuSource(result, lower(result)));
  }

  TEST_CASE("operands of different types are an error, literals take the type they meet") {
    CHECK(errorOf("out def u = f64(x) * x\n").find("are f64 and real") != std::string::npos);
    CHECK(errorOf("out def u = f64(x) * 2.0 + 1.0\n").empty());
    CHECK(errorOf("out def u = real(f64(x) * 2.0) + x\n").empty());
    // also a condition: no implicit conversion
    CHECK(errorOf("out def u = select(lt(f64(x), 1.0), x, 2.0)\n").find("are f64 and real") !=
          std::string::npos);
    CHECK(errorOf("out def u = select(real(lt(f64(x), 1.0)), x, 2.0)\n").empty());
    CHECK(errorOf("state f64 s = 0.0\nout def s = x\n").find("the state `s` is f64") !=
          std::string::npos);
    CHECK(errorOf("state f16 s = 0.0\nout def s = s\n").find("no type of a state") !=
          std::string::npos);
    CHECK(errorOf("def f64 = x\nout def u = f64\n").find("function") != std::string::npos);
    CHECK(errorOf("out def u = f64(x, x)\n").find("takes 1 argument") != std::string::npos);
  }

  TEST_CASE("binary64 where the script asks for it, in a program of binary32") {
    // t late in a run, where binary32 resolves 1.9e-6: the difference survives in binary64 only
    const std::vector<double> x = {20.0, 20.0, 1.0};
    const std::vector<double> t = {20.0000001, 20.25, 1.0000001};
    for (const auto backend : {BackendKind::Interpreter, BackendKind::RtcCpu}) {
      CAPTURE(static_cast<int>(backend));
      std::vector<double> out;
      const bool ran = evaluate("out def d = f64(t) - f64(x)\nout def e = t - x\n",
                                ComputeType::F32,
                                x,
                                t,
                                1,
                                backend,
                                out);
      if (!ran) {
        WARN_MESSAGE(ran, "the backend fell back");
        continue;
      }
      for (std::size_t p = 0; p < x.size(); ++p) {
        CAPTURE(p);
        CHECK(out[p] == t[p] - x[p]);
        CHECK(out[x.size() + p] ==
              static_cast<double>(static_cast<float>(t[p]) - static_cast<float>(x[p])));
      }
    }
  }

  TEST_CASE("a state in binary64 adds up in binary64, in a program of binary32") {
    const std::vector<double> x = {0.1, 0.25};
    const std::vector<double> t = {0.0, 0.0};
    for (const auto backend : {BackendKind::Interpreter, BackendKind::RtcCpu}) {
      CAPTURE(static_cast<int>(backend));
      std::vector<double> wide;
      std::vector<double> narrow;
      const bool ran =
          evaluate("state f64 s = 0.0\nout def s = s + f64(x)\n",
                   ComputeType::F32,
                   x,
                   t,
                   10,
                   backend,
                   wide) &&
          evaluate(
              "state s = 0.0\nout def s = s + x\n", ComputeType::F32, x, t, 10, backend, narrow);
      if (!ran) {
        WARN_MESSAGE(ran, "the backend fell back");
        continue;
      }
      for (std::size_t p = 0; p < x.size(); ++p) {
        CAPTURE(p);
        double expectedWide = 0.0;
        float expectedNarrow = 0.0F;
        for (int call = 0; call < 10; ++call) {
          expectedWide += x[p];
          expectedNarrow += static_cast<float>(x[p]);
        }
        CHECK(wide[p] == expectedWide);
        CHECK(narrow[p] == static_cast<double>(expectedNarrow));
      }
    }
  }

  TEST_CASE("the compiled kernel computes each type as the interpreter does, bit for bit") {
    const std::vector<std::string> sources = {
        // binary64 in a program of binary32, also through the maths library
        "out def u = real(sin(f64(t)) * f64(x) + exp(f64(x) / 1e3))\nout def v = sin(t) * x\n",
        // binary32 in a program of binary64
        "out def u = real(f32(x) * f32(t) - sqrt(abs(f32(x))))\nout def v = x * t\n",
        // a deep chain through the maths library, in the shape of a tile
        "def a0 = asinh(f64(x))\ndef a1 = asinh(a0 + a0 * a0)\ndef a2 = asinh(a1 + a1 * a0)\n"
        "def a3 = asinh(a2 + a2 * a1)\ndef a4 = asinh(a3 + a3 * a2)\ndef a5 = asinh(a4 + a4 * a3)\n"
        "def a6 = asinh(a5 + a5 * a4)\ndef a7 = asinh(a6 + a6 * a5)\ndef a8 = asinh(a7 + a7 * a6)\n"
        "def a9 = asinh(a8 + a8 * a7)\ndef a10 = asinh(a9 + a9 * a8)\nout def u = real(a10) + x\n",
    };
    const std::vector<double> x = {-1e3, -7.5, -1e-3, -0.0, 0.0, 0.5, 1.0, 3.7, 1e3};
    const std::vector<double> t = {20.0000001, 0.5, 1e-3, 1.0, 2.5, -3.0, 7.0, 0.25, 1e6};
    for (const auto& source : sources) {
      for (const auto type : {ComputeType::F32, ComputeType::F64}) {
        CAPTURE(source);
        CAPTURE(static_cast<int>(type));
        std::vector<double> interpreted;
        std::vector<double> compiled;
        evaluate(source, type, x, t, 1, BackendKind::Interpreter, interpreted);
        if (!evaluate(source, type, x, t, 1, BackendKind::RtcCpu, compiled)) {
          WARN_MESSAGE(false, "no usable C++ compiler; the RTC backend fell back");
          return;
        }
        CHECK(std::memcmp(
                  interpreted.data(), compiled.data(), interpreted.size() * sizeof(double)) == 0);
      }
    }
  }
}

} // namespace seissol::expr::test
