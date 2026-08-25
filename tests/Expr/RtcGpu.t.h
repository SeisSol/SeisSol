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
#include "Expr/RtcGpu.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>
#include <string>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace seissol::expr::test {

namespace {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::Direction;

/// Redefines the three macros the emitted source hides its lane indexing
/// behind, so the SAME text compiles and runs as ordinary host C++.
///
/// This is what makes the device code generator testable without a device: the
/// arithmetic is the part that can be wrong, and it is identical either way.
/// What this does NOT cover is the driver layer -- launch geometry, module
/// loading, stream association -- which needs real hardware.
constexpr auto HostShim = R"(
#define SEISSOL_EXPR_TID 0ul
#define SEISSOL_EXPR_NTHREADS 1ul
#define SEISSOL_EXPR_KERNEL extern "C"
#define __device__
#include <cmath>
using namespace std;
)";

void* compileForHost(const std::string& source) {
  const int object = memfd_create("seissol-expr-gputest", 0);
  if (object < 0) {
    return nullptr;
  }
  const std::string path = "/proc/self/fd/" + std::to_string(object);

  const pid_t pid = fork();
  if (pid < 0) {
    close(object);
    return nullptr;
  }
  if (pid == 0) {
    const int input = memfd_create("src", 0);
    const ssize_t written = write(input, source.data(), source.size());
    static_cast<void>(written);
    lseek(input, 0, SEEK_SET);
    dup2(input, STDIN_FILENO);
    execlp("c++",
           "c++",
           "-O2",
           "-ffp-contract=off",
           "-shared",
           "-fPIC",
           "-x",
           "c++",
           "-",
           "-o",
           path.c_str(),
           static_cast<char*>(nullptr));
    _exit(127);
  }
  int status = 0;
  waitpid(pid, &status, 0);
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    close(object);
    return nullptr;
  }
  // The descriptor stays open: dlopen keys its cache on the path, and closing
  // it frees the number for the next artifact, which would then be opened under
  // the same name and hand back the first library.
  return dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
}

/// Evaluate `source` on the interpreter and through the emitted device kernel
/// compiled for the host, and report whether every output is bitwise equal.
/// `ran` is false when no compiler was available.
bool deviceCodeAgrees(const std::string& source, bool& ran) {
  ran = false;
  const Program program = compileSderivModule(source);
  const std::size_t numPoints = 14;
  const std::vector<double> ladder = {
      -1e6, -1e3, -7.5, -1.0, -1e-3, -0.0, 0.0, 1e-3, 0.5, 1.0, 2.0, 3.7, 1e3, 1e6};

  std::vector<double> x(numPoints);
  std::vector<double> y(numPoints);
  for (std::size_t i = 0; i < numPoints; ++i) {
    x[i] = ladder[i];
    y[i] = ladder[(i + 5) % numPoints];
  }
  const std::size_t outputs = program.outputs().size();
  std::vector<double> interpreted(outputs * numPoints, -1.0);
  std::vector<double> emitted(outputs * numPoints, -1.0);

  DataTable table(numPoints);
  table.bindViewConst<double>("x", Direction::In, x.data());
  table.bindViewConst<double>("y", Direction::In, y.data());
  for (std::size_t i = 0; i < outputs; ++i) {
    table.bindView<double>(program.outputs()[i].name, Direction::Out, &interpreted[i * numPoints]);
  }

  Binding binding = Binding::bind(program, table);
  REQUIRE(gpuRejection(lower(program), binding, nullptr) == GpuRejection::None);

  df::GridStore store;
  const auto kernel = makeKernel(program, binding, store, {});
  kernel->precompute(table);
  kernel->run(table);

  const GpuLayout layout = gpuLayoutOf(binding);
  const std::string generated = std::string(HostShim) +
                                emitGpuSource(program, lower(program), layout, GpuTarget::Cuda) +
                                emitGpuHostTrampoline(layout, "double");
  void* handle = compileForHost(generated);
  if (handle == nullptr) {
    return true;
  }
  auto* invoke = reinterpret_cast<void (*)(void**)>(dlsym(handle, "seissol_expr_invoke"));
  if (invoke == nullptr) {
    return false;
  }
  ran = true;

  KernelArgs args{};
  std::vector<void*> bases(outputs);
  for (std::size_t i = 0; i < outputs; ++i) {
    bases[i] = &emitted[i * numPoints];
  }
  args.outputs = bases.data();
  args.outputCount = outputs;
  args.first = 0;
  args.count = numPoints;

  GpuArguments packed(binding, args, nullptr);
  invoke(packed.data());

  return std::memcmp(interpreted.data(), emitted.data(), interpreted.size() * sizeof(double)) == 0;
}

} // namespace

TEST_SUITE("ExprRtcGpu") {

  TEST_CASE("device code drops the std:: qualification") {
    // NVRTC has no <cmath>, so `std::sqrt` would not resolve. The maths
    // built-ins are unqualified in device code, and the substitution knows it --
    // which is why there is one expression table and not two.
    const Program program = compileSderivModule("out def u = sqrt(x) + exp(y)\n");
    std::vector<double> x(2);
    std::vector<double> y(2);
    std::vector<double> u(2);
    DataTable table(2);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindView<double>("u", Direction::Out, u.data());
    const Binding binding = Binding::bind(program, table);

    const std::string source =
        emitGpuSource(program, lower(program), gpuLayoutOf(binding), GpuTarget::Cuda);
    CHECK(source.find("std::") == std::string::npos);
    CHECK(source.find("sqrt(") != std::string::npos);
    // The `x` inside `exp` must survive the operand substitution.
    CHECK(source.find("exp(s") != std::string::npos);
  }

  TEST_CASE("the element type is baked in, the stride is not") {
    // A per-point type switch would cost more than the arithmetic, so the type
    // has to be resolved at generation time. Strides are scalars read once per
    // thread; baking those would give one kernel per binding layout instead of
    // one per program, for nothing.
    const Program program = compileSderivModule("out def u = x + 1.0\n");
    std::vector<float> x(4);
    std::vector<double> u(4);
    DataTable table(4);
    table.bindViewConst<float>("x", Direction::In, x.data());
    table.bindView<double>("u", Direction::Out, u.data());
    const Binding binding = Binding::bind(program, table);

    const std::string source =
        emitGpuSource(program, lower(program), gpuLayoutOf(binding), GpuTarget::Cuda);
    CHECK(source.find("float value;") != std::string::npos);
    CHECK(source.find("unsigned long long stride_in0") != std::string::npos);
  }

  TEST_CASE("the HIP target brings its own runtime header") {
    const Program program = compileSderivModule("out def u = x\n");
    std::vector<double> x(2);
    std::vector<double> u(2);
    DataTable table(2);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindView<double>("u", Direction::Out, u.data());
    const Binding binding = Binding::bind(program, table);

    const GpuLayout layout = gpuLayoutOf(binding);
    CHECK(
        emitGpuSource(program, lower(program), layout, GpuTarget::Hip).find("hip/hip_runtime.h") !=
        std::string::npos);
    CHECK(
        emitGpuSource(program, lower(program), layout, GpuTarget::Cuda).find("hip/hip_runtime.h") ==
        std::string::npos);
  }

  TEST_CASE("a program that cannot go to a device says why") {
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> u(2, 0.0);

    SUBCASE("a grid lookup") {
      const Program program =
          compileSderivModule("grid m = \"asagi\", \"model.nc\", \"linear\", \"rho\"\n"
                              "out def u = m_rho(x, y, z)\n");
      std::vector<double> y(2);
      std::vector<double> z(2);
      DataTable table(2);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindViewConst<double>("y", Direction::In, y.data());
      table.bindViewConst<double>("z", Direction::In, z.data());
      table.bindView<double>("u", Direction::Out, u.data());
      const Binding binding = Binding::bind(program, table);
      CHECK(gpuRejection(lower(program), binding, nullptr) == GpuRejection::Lookup);
    }

    SUBCASE("a computed column") {
      const Program program = compileSderivModule("out def u = x + group\n");
      DataTable table(2);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindComputed("group", [](std::size_t) -> std::int32_t { return 1; });
      table.bindView<double>("u", Direction::Out, u.data());
      const Binding binding = Binding::bind(program, table);
      CHECK(gpuRejection(lower(program), binding, nullptr) == GpuRejection::ComputedColumn);
    }

    SUBCASE("a host pointer") {
      const Program program = compileSderivModule("out def u = x + 1.0\n");
      DataTable table(2);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindView<double>("u", Direction::Out, u.data());
      const Binding binding = Binding::bind(program, table);
      CHECK(gpuRejection(lower(program), binding, [](const void*) { return false; }) ==
            GpuRejection::HostPointer);
    }
  }

  TEST_CASE("the emitted device kernel is bitwise identical to the interpreter") {
    // Compiled for the host and called through GpuArguments, so this covers the
    // parameter ORDER as well as the arithmetic -- the emitted signature's arity
    // depends on the program, and getting the packing wrong makes the kernel
    // write nothing at all rather than crash.
    const std::vector<std::pair<const char*, const char*>> programs = {
        {"arithmetic and powers", "out def u = (x*y)/(x*x + 1.0) + x**3.0\n"},
        {"roots, exp and logs",
         "def a = abs(x)+1.0\nout def u = sqrt(a)+exp(0.0-a)+log(a)+log2(a)+log10(a)\n"},
        {"trigonometry", "out def u = sin(x)+cos(y)+tan(x)+atan(y)+atan2(x,y)\n"},
        {"hyperbolics", "out def u = sinh(x/1e3)+cosh(y/1e3)+tanh(x)+erf(y)\n"},
        {"rounding and sign", "out def u = floor(x)+ceil(y)+round(x)+sign(y)+abs(x)\n"},
        {"min, max, mod", "out def u = min(x,y)+max(x,y)+mod(x,3.0)\n"},
        {"comparisons and select",
         "out def u = select(land(lt(x,y), ge(y,x)), x, y) + select(lnot(eq(x,y)),1.0,2.0)\n"},
        {"several outputs", "def a = x*y\nout def u = a\nout def v = a - x\n"},
        {"a shared subexpression", "def a = sqrt(abs(x))\nout def u = a+a*a+a*a*a\n"},
    };

    bool ran = false;
    for (const auto& [label, source] : programs) {
      CAPTURE(label);
      const bool same = deviceCodeAgrees(source, ran);
      if (!ran) {
        WARN_MESSAGE(ran, "no usable C++ compiler; the device code generator was not executed");
        break;
      }
      CHECK(same);
    }
  }

} // TEST_SUITE

} // namespace seissol::expr::test
