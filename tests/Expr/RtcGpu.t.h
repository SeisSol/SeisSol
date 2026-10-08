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
#include "Expr/Rewrite.h"
#include "Expr/RtcGpu.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "TestHelper.h"

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

/// Whether NVRTC, if this machine has it, compiles `source` as the driver does (set `ran`
/// accordingly; `log` gets its diagnostics). Loaded at run time, so the test needs no toolkit to
/// build and is skipped where there is none -- where there is one, the CUDA dialect of the
/// emitted source is checked by the real compiler rather than by the host shim alone.
bool nvrtcAccepts(const std::string& source, bool& ran, std::string& log) {
  ran = false;
  void* library = dlopen("libnvrtc.so", RTLD_NOW | RTLD_LOCAL);
  if (library == nullptr) {
    library = dlopen("libnvrtc.so.12", RTLD_NOW | RTLD_LOCAL);
  }
  if (library == nullptr) {
    return true;
  }
  using Create =
      int (*)(void**, const char*, const char*, int, const char* const*, const char* const*);
  using Compile = int (*)(void*, int, const char* const*);
  using LogSize = int (*)(void*, std::size_t*);
  using Log = int (*)(void*, char*);
  using Destroy = int (*)(void**);
  auto* create = reinterpret_cast<Create>(dlsym(library, "nvrtcCreateProgram"));
  auto* compile = reinterpret_cast<Compile>(dlsym(library, "nvrtcCompileProgram"));
  auto* logSize = reinterpret_cast<LogSize>(dlsym(library, "nvrtcGetProgramLogSize"));
  auto* getLog = reinterpret_cast<Log>(dlsym(library, "nvrtcGetProgramLog"));
  auto* destroy = reinterpret_cast<Destroy>(dlsym(library, "nvrtcDestroyProgram"));
  if (create == nullptr || compile == nullptr || logSize == nullptr || getLog == nullptr ||
      destroy == nullptr) {
    return true;
  }
  void* program = nullptr;
  if (create(&program, source.c_str(), "seissol_expr.cu", 0, nullptr, nullptr) != 0) {
    return true;
  }
  ran = true;
  const std::vector<const char*> options = {
      "--gpu-architecture=compute_70", "--std=c++17", "-default-device"};
  const bool compiled = compile(program, static_cast<int>(options.size()), options.data()) == 0;
  std::size_t size = 0;
  logSize(program, &size);
  log.assign(size, '\0');
  getLog(program, log.data());
  destroy(&program);
  return compiled;
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
  REQUIRE(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);

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
    // A typed load, not a byte copy: gpuRejection() has established that the
    // address is element-aligned, so the cast is safe and a frontend that
    // lowers __builtin_memcpy badly cannot hurt us.
    CHECK(source.find("(const float*)bytes") != std::string::npos);
    // SeissolU64, because the dialects spell a 64-bit unsigned differently and
    // `unsigned long` is 32-bit on a Windows host.
    CHECK(source.find("SeissolU64 stride_in0") != std::string::npos);
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
      CHECK(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::Lookup);
    }

    SUBCASE("a computed column") {
      const Program program = compileSderivModule("out def u = x + group\n");
      DataTable table(2);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindComputed("group", [](std::size_t) -> std::int32_t { return 1; });
      table.bindView<double>("u", Direction::Out, u.data());
      const Binding binding = Binding::bind(program, table);
      CHECK(gpuRejection(program, lower(program), binding, nullptr) ==
            GpuRejection::ComputedColumn);
    }

    SUBCASE("a host pointer") {
      const Program program = compileSderivModule("out def u = x + 1.0\n");
      DataTable table(2);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindView<double>("u", Direction::Out, u.data());
      const Binding binding = Binding::bind(program, table);
      CHECK(gpuRejection(program, lower(program), binding, [](const void*) { return false; }) ==
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
    for (const auto& program : programs) {
      // (a structured binding cannot be captured where OpenMP is enabled)
      const char* const label = program.first;
      const char* const source = program.second;
      CAPTURE(label);
      const bool same = deviceCodeAgrees(source, ran);
      if (!ran) {
        WARN_MESSAGE(ran, "no usable C++ compiler; the device code generator was not executed");
        break;
      }
      CHECK(same);
    }
  }

  TEST_CASE("the emitted device kernel reads and writes a state the table keeps in place") {
    // Compiled for the host and called through GpuArguments: the state descriptors in the
    // argument block, a gathered subset of cells, and a non-zero initial value, which a state the
    // kernel kept itself could not have.
    const Program program = compileSderivModule("state peak = -1.0\n"
                                                "state total = 0.5\n"
                                                "out def peak = max(peak, x)\n"
                                                "out def ratio = total / (1.0 + abs(peak))\n"
                                                "def total = total + x * y\n");
    constexpr std::size_t Rows = 2;
    constexpr std::size_t NumPoints = 3 * Rows;
    constexpr std::size_t CellValues = 2 * Rows;
    const std::vector<std::uint32_t> cellIndex = {2, 0, 3};
    std::vector<double> x(NumPoints);
    std::vector<double> y(NumPoints);
    const auto keptStates = [&]() {
      std::vector<double> storage(4 * CellValues);
      for (std::size_t cell = 0; cell < 4; ++cell) {
        for (std::size_t row = 0; row < Rows; ++row) {
          storage[cell * CellValues + row] = program.state()[0].initial;
          storage[cell * CellValues + Rows + row] = program.state()[1].initial;
        }
      }
      return storage;
    };
    std::vector<double> interpretedStates = keptStates();
    std::vector<double> emittedStates = keptStates();

    std::vector<double> interpreted(2 * NumPoints, -1.0);
    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindView<double>("peak", Direction::Out, interpreted.data());
    table.bindView<double>("ratio", Direction::Out, interpreted.data() + NumPoints);
    table.bindState<double>(
        "peak", interpretedStates.data(), Rows, CellValues, 1, cellIndex.data());
    table.bindState<double>(
        "total", interpretedStates.data() + Rows, Rows, CellValues, 1, cellIndex.data());
    Binding binding = Binding::bind(program, table);
    REQUIRE(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);
    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);

    const GpuLayout layout = gpuLayoutOf(binding);
    CHECK(layout.states == 2);
    const std::string generated = std::string(HostShim) +
                                  emitGpuSource(program, lower(program), layout, GpuTarget::Cuda) +
                                  emitGpuHostTrampoline(layout, "double");
    void* handle = compileForHost(generated);
    if (handle == nullptr) {
      WARN_MESSAGE(false, "no usable C++ compiler; the device code generator was not executed");
      return;
    }
    auto* invoke = reinterpret_cast<void (*)(void**)>(dlsym(handle, "seissol_expr_invoke"));
    REQUIRE(invoke != nullptr);

    std::vector<double> emitted(2 * NumPoints, -1.0);
    for (std::size_t call = 0; call < 3; ++call) {
      CAPTURE(call);
      for (std::size_t p = 0; p < NumPoints; ++p) {
        x[p] = std::cos(0.9 * static_cast<double>(call + p)) * static_cast<double>(p);
        y[p] = 0.3 * static_cast<double>(call) - 0.1 * static_cast<double>(p);
      }
      kernel->run(table);

      KernelArgs args{};
      std::vector<void*> outputs = {emitted.data(), emitted.data() + NumPoints};
      std::vector<void*> states = {emittedStates.data(), emittedStates.data() + Rows};
      args.outputs = outputs.data();
      args.outputCount = outputs.size();
      args.states = states.data();
      args.stateCount = states.size();
      args.first = 0;
      args.count = NumPoints;
      GpuArguments packed(binding, args, nullptr);
      invoke(packed.data());
      // five fields per state in the argument block
      CHECK(packed.fieldCount() == 3 * (2 + 2) + 5 * 2 + 4);

      CHECK(unit_test::bitwiseEqual(interpreted.data(), emitted.data(), interpreted.size()));
      CHECK(unit_test::bitwiseEqual(
          interpretedStates.data(), emittedStates.data(), interpretedStates.size()));
    }

    // a state with a non-zero initial value that the kernel had to keep itself
    DataTable own(NumPoints);
    own.bindViewConst<double>("x", Direction::In, x.data());
    own.bindViewConst<double>("y", Direction::In, y.data());
    own.bindView<double>("peak", Direction::Out, interpreted.data());
    own.bindView<double>("ratio", Direction::Out, interpreted.data() + NumPoints);
    const Binding ownBinding = Binding::bind(program, own);
    CHECK(gpuRejection(program, lower(program), ownBinding, nullptr) ==
          GpuRejection::StatefulProgram);
  }

  TEST_CASE("the emitted device kernel contracts bit for bit as the interpreter does") {
    // A contraction compiled for the host and called through GpuArguments: covers the block and
    // matrix descriptors in the argument block as well as the loop.
    Program program =
        compileSderivModule("out def u = 2.0 * v - w * x\nout def r = sqrt(abs(v))\n");
    constexpr std::size_t Rows = 2;
    constexpr std::size_t Cols = 3;
    constexpr std::size_t Ld = 4;
    constexpr std::size_t NumCells = 3;
    constexpr std::size_t NumPoints = Rows * NumCells;
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto v = program.internBlock("v", Cols);
    const auto w = program.internBlock("w", Cols);
    substituteByContraction(program, matrix, {{"v", v}, {"w", w}});

    std::vector<double> m(Ld * Cols);
    std::vector<float> dofs(std::size_t{4} * 2 *
                            Cols); // f32 coefficients, [cell][quantity][coefficient]
    std::vector<double> x(NumPoints);
    for (std::size_t i = 0; i < m.size(); ++i) {
      m[i] = 0.25 * static_cast<double>(i) - 1.3;
    }
    for (std::size_t i = 0; i < dofs.size(); ++i) {
      dofs[i] = static_cast<float>(1.7 * static_cast<double>(i % 5) - 2.1);
    }
    for (std::size_t i = 0; i < x.size(); ++i) {
      x[i] = 0.5 + static_cast<double>(i);
    }
    const std::vector<std::uint32_t> cellIndex = {3, 0, 2};

    std::vector<double> interpreted(2 * NumPoints, -1.0);
    std::vector<double> emitted(2 * NumPoints, -1.0);
    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindBlock<float>("v", dofs.data(), Cols, 2 * Cols, 1, cellIndex.data());
    table.bindBlock<float>("w", dofs.data() + Cols, Cols, 2 * Cols, 1, cellIndex.data());
    table.bindMatrix<double>("proj", m.data(), Rows, Cols, Ld);
    table.bindView<double>("u", Direction::Out, interpreted.data());
    table.bindView<double>("r", Direction::Out, interpreted.data() + NumPoints);

    Binding binding = Binding::bind(program, table);
    REQUIRE(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);
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
      WARN_MESSAGE(false, "no usable C++ compiler; the device code generator was not executed");
      return;
    }
    auto* invoke = reinterpret_cast<void (*)(void**)>(dlsym(handle, "seissol_expr_invoke"));
    REQUIRE(invoke != nullptr);

    KernelArgs args{};
    std::vector<void*> bases = {emitted.data(), emitted.data() + NumPoints};
    args.outputs = bases.data();
    args.outputCount = 2;
    args.first = 0;
    args.count = NumPoints;
    GpuArguments packed(binding, args, nullptr);
    invoke(packed.data());

    CHECK(std::memcmp(interpreted.data(), emitted.data(), interpreted.size() * sizeof(double)) ==
          0);

    // one pointer per matrix and four fields per block in the argument block
    CHECK(packed.fieldCount() == 3 * (1 + 2) + 1 + 4 * 2 + 4);
  }

  TEST_CASE("a column per cell reaches the device kernel with its divisor and index") {
    const Program program = compileSderivModule("out def u = x * j\n");
    constexpr std::size_t PointsPerCell = 2;
    constexpr std::size_t NumPoints = 3 * PointsPerCell;
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const std::vector<float> j = {10.0F, 20.0F, 30.0F, 40.0F};
    const std::vector<std::uint32_t> cellIndex = {3, 1, 0};
    std::vector<double> interpreted(NumPoints, -1.0);
    std::vector<double> emitted(NumPoints, -1.0);

    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindCellView<float>("j", j.data(), PointsPerCell, 1, 0, cellIndex.data());
    table.bindView<double>("u", Direction::Out, interpreted.data());
    Binding binding = Binding::bind(program, table);
    REQUIRE(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);
    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);
    for (std::size_t p = 0; p < NumPoints; ++p) {
      CHECK(interpreted[p] == x[p] * j[cellIndex[p / PointsPerCell]]);
    }

    const GpuLayout layout = gpuLayoutOf(binding);
    const std::string generated = std::string(HostShim) +
                                  emitGpuSource(program, lower(program), layout, GpuTarget::Cuda) +
                                  emitGpuHostTrampoline(layout, "double");
    void* handle = compileForHost(generated);
    if (handle == nullptr) {
      WARN_MESSAGE(false, "no usable C++ compiler; the device code generator was not executed");
      return;
    }
    auto* invoke = reinterpret_cast<void (*)(void**)>(dlsym(handle, "seissol_expr_invoke"));
    REQUIRE(invoke != nullptr);
    KernelArgs args{};
    void* base = emitted.data();
    args.outputs = &base;
    args.outputCount = 1;
    args.first = 0;
    args.count = NumPoints;
    GpuArguments packed(binding, args, nullptr);
    invoke(packed.data());
    CHECK(unit_test::bitwiseEqual(interpreted.data(), emitted.data(), NumPoints));
    // a column per cell carries its divisor and index along
    CHECK(packed.fieldCount() == 3 * 3 + 2 + 4);
  }

  TEST_CASE("a uniform input reaches the device kernel by value, per call") {
    // `t` is bound as a constant and moved per call to a host variable: the device kernel takes
    // its value in the argument block, so neither the binding nor the call needs device memory.
    for (const auto computeType : {ComputeType::F64, ComputeType::F32}) {
      Program program = compileSderivModule("out def u = x * t + t\n");
      program.setComputeType(computeType);
      const bool f32 = computeType == ComputeType::F32;
      constexpr std::size_t NumPoints = 5;
      std::vector<double> x = {1.0, -2.0, 0.5, 3.0, 1e3};
      std::vector<double> interpreted(NumPoints, -1.0);
      std::vector<double> emitted(NumPoints, -1.0);

      DataTable table(NumPoints);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindConstant<double>("t", 0.0);
      table.bindView<double>("u", Direction::Out, interpreted.data());
      Binding binding = Binding::bind(program, table);
      // the constant lives on the host, which a uniform input may
      const auto nothingIsOnTheDevice = [](const void* /*pointer*/) { return false; };
      CHECK(gpuRejection(program, lower(program), binding, +nothingIsOnTheDevice) ==
            GpuRejection::HostPointer);
      const GpuLayout layout = gpuLayoutOf(binding);
      REQUIRE(layout.inputUniform.size() == 2);
      CHECK(!layout.inputUniform[0]);
      CHECK(layout.inputUniform[1]);

      const std::string source = emitGpuSource(program, lower(program), layout, GpuTarget::Cuda);
      CHECK(source.find(std::string(f32 ? "float" : "double") + " value_in1;") !=
            std::string::npos);
      const std::string generated =
          std::string(HostShim) + source + emitGpuHostTrampoline(layout, f32 ? "float" : "double");
      void* handle = compileForHost(generated);
      if (handle == nullptr) {
        WARN_MESSAGE(false, "no usable C++ compiler; the device code generator was not executed");
        return;
      }
      auto* invoke = reinterpret_cast<void (*)(void**)>(dlsym(handle, "seissol_expr_invoke"));
      REQUIRE(invoke != nullptr);

      df::GridStore store;
      const auto kernel = makeKernel(program, binding, store, {});
      for (const double time : {0.25, -3.0}) {
        // the time is the base of the call, for the interpreter and the device kernel alike
        std::vector<const void*> inputs = {nullptr, &time};
        KernelArgs args{};
        args.inputs = inputs.data();
        args.inputCount = inputs.size();
        args.first = 0;
        args.count = NumPoints;
        void* interpretedBase = interpreted.data();
        args.outputs = &interpretedBase;
        args.outputCount = 1;
        kernel->run(args);

        void* emittedBase = emitted.data();
        args.outputs = &emittedBase;
        GpuArguments packed(binding, args, nullptr);
        invoke(packed.data());
        CHECK(unit_test::bitwiseEqual(interpreted.data(), emitted.data(), NumPoints));
        for (std::size_t p = 0; p < NumPoints; ++p) {
          const double expected = x[p] * time + time;
          if (f32) {
            CHECK(emitted[p] == doctest::Approx(expected).epsilon(1e-6));
          } else {
            CHECK(emitted[p] == expected);
          }
        }
        // one field for the uniform input, eight bytes wide whatever the compute type
        CHECK(packed.fieldCount() == 3 + 1 + 3 + 4);
        CHECK(packed.fieldSize(3) == (f32 ? sizeof(float) : sizeof(double)));
        CHECK(static_cast<const char*>(packed.fieldData(4)) -
                  static_cast<const char*>(packed.fieldData(3)) ==
              8);
      }
    }
  }

  TEST_CASE("NVRTC compiles the emitted CUDA source") {
    // Columns per point and per cell, a uniform input, a state, a contraction, in both compute
    // types: every form the argument block and the accessors take.
    for (const auto computeType : {ComputeType::F64, ComputeType::F32}) {
      Program program = compileSderivModule(
          "state s = 0.0\nout def s = s + v * t\nout def u = 2.0 * v - j * x\n");
      program.setComputeType(computeType);
      const auto matrix = program.internMatrix("proj", MatrixShape{2, 3, 4});
      const auto v = program.internBlock("v", 3);
      substituteByContraction(program, matrix, {{"v", v}});

      constexpr std::size_t NumPoints = 4;
      std::vector<double> x(NumPoints, 1.0);
      std::vector<float> j(2, 2.0F);
      std::vector<double> m(12, 0.5);
      std::vector<double> dofs(6, 1.0);
      std::vector<double> state(NumPoints, 0.0);
      std::vector<float> stateF32(NumPoints, 0.0F);
      std::vector<double> s(NumPoints);
      std::vector<double> u(NumPoints);
      DataTable table(NumPoints);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindCellView<float>("j", j.data(), 2);
      table.bindConstant<double>("t", 0.5);
      table.bindBlock<double>("v", dofs.data(), 3, 3);
      table.bindMatrix<double>("proj", m.data(), 2, 3, 4);
      if (computeType == ComputeType::F32) {
        table.bindState<float>("s", stateF32.data(), 2, 2, 1);
      } else {
        table.bindState<double>("s", state.data(), 2, 2, 1);
      }
      table.bindView<double>("s", Direction::Out, s.data());
      table.bindView<double>("u", Direction::Out, u.data());
      const Binding binding = Binding::bind(program, table);
      const std::string source =
          emitGpuSource(program, lower(program), gpuLayoutOf(binding), GpuTarget::Cuda);
      bool ran = false;
      std::string log;
      const bool accepted = nvrtcAccepts(source, ran, log);
      if (!ran) {
        WARN_MESSAGE(false, "no NVRTC on this machine; the CUDA dialect was not compiled");
        return;
      }
      INFO(log);
      INFO(source);
      CHECK(accepted);
    }
  }

  TEST_CASE("the kernel splits into a point function and a wrapper") {
    // The split is what lets a kernel that already owns a loop over the same
    // points -- a batched contraction, say -- call the expression from inside,
    // with the interior state still in registers. Without it, fusing a
    // projection would mean teaching this generator what a GEMM is.
    const Program program = compileSderivModule("out def u = x * 2.0\n");
    std::vector<double> x(4);
    std::vector<double> u(4);
    DataTable table(4);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindView<double>("u", Direction::Out, u.data());
    const Binding binding = Binding::bind(program, table);

    const std::string source =
        emitGpuSource(program, lower(program), gpuLayoutOf(binding), GpuTarget::Cuda);
    CHECK(source.find("__device__ inline void seissol_expr_run_point") != std::string::npos);
    CHECK(source.find("SEISSOL_EXPR_KERNEL void seissol_expr_run(") != std::string::npos);
    // The wrapper must be a loop over the point function and nothing else, or
    // the two would be two implementations of the same expression.
    CHECK(source.find("seissol_expr_run_point(&a, a.first + l)") != std::string::npos);
  }

  TEST_CASE("the OpenCL dialect differs where it has to and nowhere else") {
    // sign, the comparisons and the constant cast to the compute type, which OpenCL C can only
    // spell as a C cast
    const Program program =
        compileSderivModule("out def u = sqrt(x) + select(lt(y, 0.5), sign(y), y)\n");
    std::vector<double> x(2);
    std::vector<double> y(2);
    std::vector<double> u(2);
    DataTable table(2);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindView<double>("u", Direction::Out, u.data());
    const Binding binding = Binding::bind(program, table);
    const GpuLayout layout = gpuLayoutOf(binding);

    const std::string opencl = emitGpuSource(program, lower(program), layout, GpuTarget::OpenCl);
    const std::string cuda = emitGpuSource(program, lower(program), layout, GpuTarget::Cuda);
    // no functional casts and no C++ casts
    CHECK(opencl.find("double(") == std::string::npos);
    CHECK(opencl.find("static_cast") == std::string::npos);

    // Address spaces are required in OpenCL C and inferred everywhere else.
    CHECK(opencl.find("__global const void* in0") != std::string::npos);
    // Not "__global": the CUDA keyword __global__ contains it. The OpenCL
    // address space is what must be absent.
    CHECK(cuda.find("__global const void*") == std::string::npos);
    // f64 is an extension there, and needed even for an f32 program that
    // carries one f64 column.
    CHECK(opencl.find("cl_khr_fp64") != std::string::npos);
    // No `unsigned long` anywhere: 32-bit on Windows, and not an OpenCL C type.
    CHECK(opencl.find("unsigned long") == std::string::npos);
    CHECK(opencl.find("typedef ulong SeissolU64") != std::string::npos);
    // A kernel-argument struct may not portably hold pointers, so this target
    // spells the parameters out and gathers them into a local struct.
    // The kernel qualifier goes through the macro, so both the definition and
    // the use have to be there.
    CHECK(opencl.find("#define SEISSOL_EXPR_KERNEL __kernel") != std::string::npos);
    CHECK(opencl.find("SEISSOL_EXPR_KERNEL void seissol_expr_run(") != std::string::npos);
    CHECK(opencl.find("SeissolExprArgs a;") != std::string::npos);
    // The point function is identical in shape on both, which is what keeps
    // the fusion seam target-independent.
    CHECK(opencl.find("seissol_expr_run_point(&a, a.first + l)") != std::string::npos);
    CHECK(cuda.find("seissol_expr_run_point(&a, a.first + l)") != std::string::npos);
  }

  TEST_CASE("both argument views come from one packing") {
    const Program program = compileSderivModule("out def u = x + y\n");
    std::vector<double> x = {1.0};
    std::vector<double> y = {2.0};
    std::vector<double> u = {0.0};
    DataTable table(1);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindView<double>("u", Direction::Out, u.data());
    const Binding binding = Binding::bind(program, table);

    KernelArgs args{};
    args.first = 0;
    args.count = 1;
    GpuArguments packed(binding, args, nullptr);

    // Two inputs and one output at three fields each, plus persistent,
    // numPoints, first and count.
    CHECK(packed.fieldCount() == 3 * 3 + 4);
    // The struct view is one entry pointing at the whole image; the flat view
    // is one entry per field of that same image. They cannot disagree because
    // there is one packing behind both.
    CHECK(packed.size() == 1);
    CHECK(packed.fieldData(0) == *packed.data());
  }

} // TEST_SUITE

} // namespace seissol::expr::test
