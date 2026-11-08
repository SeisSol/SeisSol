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
#include "Expr/Rewrite.h"
#include "Expr/RtcGpu.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "Reader/Scripting/DataTable.h"
#include "TestHelper.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>
#include <limits>
#include <memory>
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
#define SEISSOL_EXPR_DADD(x, y) ((x) + (y))
#define SEISSOL_EXPR_DSUB(x, y) ((x) - (y))
#define SEISSOL_EXPR_DMUL(x, y) ((x) * (y))
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
      "--gpu-architecture=compute_75", "--std=c++17", "-default-device"};
  const bool compiled = compile(program, static_cast<int>(options.size()), options.data()) == 0;
  std::size_t size = 0;
  logSize(program, &size);
  log.assign(size, '\0');
  getLog(program, log.data());
  destroy(&program);
  return compiled;
}

/// A grid over an array of its own, sampled by the arithmetic of Interpolation.h, which it exposes
/// as the array a device kernel copies (Grid::arrayView()). Values that no scheme reproduces
/// exactly, so that every weight and every rounding shows in the samples.
class ArrayGrid final : public df::Grid {
  public:
  static constexpr unsigned Components = 2;

  ArrayGrid(const std::vector<unsigned>& num, df::ElementType type, df::Interpolation interpolation)
      : interpolation_(interpolation) {
    view_.type = type;
    view_.dimensions = static_cast<unsigned>(num.size());
    view_.components = Components;
    view_.componentStride = 1;
    std::size_t stride = Components;
    for (std::size_t d = 0; d < num.size(); ++d) {
      view_.min[d] = -1.0 + 0.3 * static_cast<double>(d);
      delta_[d] = 0.5 + 0.25 * static_cast<double>(d);
      view_.deltaInv[d] = 1.0 / delta_[d];
      view_.num[d] = num[d];
      view_.stride[d] = stride;
      stride *= num[d];
    }
    values_.resize(stride);
    for (std::size_t i = 0; i < stride; ++i) {
      values_[i] = std::sin(0.37 * static_cast<double>(i) + 0.11) * 10.0 +
                   0.013 * static_cast<double>(i % 7);
    }
    if (type == df::ElementType::Float) {
      floats_.assign(values_.begin(), values_.end());
      view_.data = floats_.data();
    } else {
      view_.data = values_.data();
    }
  }

  [[nodiscard]] std::size_t dimensions() const override { return view_.dimensions; }
  [[nodiscard]] std::size_t components() const override { return Components; }

  void geometry(df::GridGeometry& out) const override {
    out.dimensions = view_.dimensions;
    for (unsigned d = 0; d < view_.dimensions; ++d) {
      out.min[d] = view_.min[d];
      out.delta[d] = delta_[d];
      out.num[d] = view_.num[d];
    }
  }

  void sample(const double* x, double* out) const override {
    df::sample(view_, interpolation_, x, out);
  }
  void sampleBatch(const double* const* x, std::size_t count, double* const* out) const override {
    df::sampleBatch(view_, interpolation_, x, count, out, scratch_);
  }
  void sampleBatch(const double* const* x, std::size_t count, float* const* out) const override {
    df::sampleBatch(view_, interpolation_, x, count, out, scratch_);
  }
  [[nodiscard]] const df::ArrayView* arrayView() const override { return &view_; }

  private:
  df::Interpolation interpolation_;
  df::ArrayView view_;
  double delta_[df::MaxGridDimensions]{};
  std::vector<double> values_;
  std::vector<float> floats_;
  mutable df::SampleScratch scratch_;
};

/// The program that samples both components of a grid `desc` of `dimensions` axes at coordinates
/// computed from the inputs x and y, in the compute type `type`, and gives the samples and a value
/// computed from one; and its grid.
Program lookupProgram(const df::GridDesc& desc, unsigned dimensions, ComputeType type) {
  Program program;
  program.setComputeType(type);
  const GridId grid = program.internGrid(desc);
  program.addInput("x", reader::scripting::DataType::F64);
  program.addInput("y", reader::scripting::DataType::F64);
  auto& arena = program.arena();
  const NodeId x = arena.field("x");
  const NodeId y = arena.field("y");
  // the coordinates, each its own expression: x, y, and mixes of both beyond
  std::vector<NodeId> coordinates = {x, y};
  coordinates.push_back(arena.pw(Fn::Add, arena.pw(Fn::Mul, x, arena.konst(0.5)), y));
  coordinates.push_back(arena.pw(Fn::Sub, y, arena.pw(Fn::Mul, x, arena.konst(0.25))));
  coordinates.resize(dimensions);
  const NodeId first = arena.lookup(grid, 0, coordinates);
  const NodeId second = arena.lookup(grid, 1, coordinates);
  program.addOutput("u", reader::scripting::DataType::F64, first);
  program.addOutput("v", reader::scripting::DataType::F64, second);
  program.addOutput("w", reader::scripting::DataType::F64, arena.pw(Fn::Mul, first, x));
  return program;
}

/// The outputs of `program`, which reads x and y and samples the grid `desc`, an ArrayGrid, on
/// `backend` in the shape `shape`; whether that backend ran.
bool evaluateLookups(const Program& program,
                     const df::GridDesc& desc,
                     const std::vector<unsigned>& num,
                     df::ElementType element,
                     const std::vector<double>& x,
                     const std::vector<double>& y,
                     BackendKind backend,
                     KernelShape shape,
                     std::vector<double>& out) {
  const std::size_t numPoints = x.size();
  const std::size_t outputs = program.outputs().size();
  out.assign(outputs * numPoints, -1.0);
  DataTable table(numPoints);
  table.bindViewConst<double>("x", Direction::In, x.data());
  table.bindViewConst<double>("y", Direction::In, y.data());
  for (std::size_t i = 0; i < outputs; ++i) {
    table.bindView<double>(program.outputs()[i].name, Direction::Out, &out[i * numPoints]);
  }
  Binding binding = Binding::bind(program, table);
  df::GridStore store;
  const std::size_t id = store.intern(desc);
  store.injectForTesting(id, std::make_unique<ArrayGrid>(num, element, desc.interpolation));
  BackendOptions options;
  options.preferred = backend;
  options.shape = shape;
  options.quiet = true;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);
  return kernel->kind() == backend;
}

/// The samples of a grid by the interpreter, by the compiled CPU kernel and by the emitted device
/// kernel compiled for the host, at points inside, on the edges and outside of the grid (and NaN
/// and infinities): whether they are the same bits. `ran` is false when no compiler was available;
/// `outside` is what the device kernel left in GridDescriptor::outside. The backends of the host
/// evaluate the program with a grid that clamps: one that fails would stop the test at a point
/// beyond it, and the samples are the same.
bool lookupsAgree(df::Interpolation interpolation,
                  df::BoundaryMode boundary,
                  const std::vector<unsigned>& num,
                  df::ElementType element,
                  ComputeType type,
                  const std::vector<double>& x,
                  const std::vector<double>& y,
                  KernelShape shape,
                  bool& ran,
                  std::uint64_t& outside,
                  bool withCpu = true) {
  ran = false;
  df::GridDesc desc;
  desc.path = "grid.nc";
  desc.variable = "z";
  desc.interpolation = interpolation;
  desc.boundary = boundary;
  df::GridDesc clamping = desc;
  clamping.boundary = df::BoundaryMode::Clamp;
  const auto dimensions = static_cast<unsigned>(num.size());
  const Program program = lookupProgram(desc, dimensions, type);
  const Program clampingProgram = lookupProgram(clamping, dimensions, type);

  std::vector<double> interpreted;
  evaluateLookups(clampingProgram,
                  clamping,
                  num,
                  element,
                  x,
                  y,
                  BackendKind::Interpreter,
                  KernelShape::Auto,
                  interpreted);
  bool same = true;
  std::vector<double> compiled;
  if (withCpu &&
      evaluateLookups(
          clampingProgram, clamping, num, element, x, y, BackendKind::RtcCpu, shape, compiled)) {
    same = unit_test::bitwiseEqual(interpreted.data(), compiled.data(), interpreted.size());
  }

  const std::size_t numPoints = x.size();
  const std::size_t outputs = program.outputs().size();
  std::vector<double> emitted(outputs * numPoints, -1.0);
  DataTable table(numPoints);
  table.bindViewConst<double>("x", Direction::In, x.data());
  table.bindViewConst<double>("y", Direction::In, y.data());
  for (std::size_t i = 0; i < outputs; ++i) {
    table.bindView<double>(program.outputs()[i].name, Direction::Out, &emitted[i * numPoints]);
  }
  Binding binding = Binding::bind(program, table);
  REQUIRE(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);

  GpuLayout layout = gpuLayoutOf(binding);
  layout.grids = {element};
  const LoweredProgram lowered = lower(program);
  const std::string generated =
      std::string(HostShim) + emitGpuSource(program, lowered, layout, GpuTarget::Cuda) +
      emitGpuHostTrampoline(layout, codegen::typeName(lowered.dataType()));
  void* handle = compileForHost(generated);
  if (handle == nullptr) {
    return true;
  }
  auto* invoke = reinterpret_cast<void (*)(void**)>(dlsym(handle, "seissol_expr_invoke"));
  if (invoke == nullptr) {
    return false;
  }
  ran = true;

  const ArrayGrid grid(num, element, interpolation);
  const df::ArrayView* view = grid.arrayView();
  codegen::GridDescriptor descriptor = codegen::gridDescriptor(*view, view->data);
  KernelArgs args{};
  std::vector<void*> bases(outputs);
  for (std::size_t i = 0; i < outputs; ++i) {
    bases[i] = &emitted[i * numPoints];
  }
  args.outputs = bases.data();
  args.outputCount = outputs;
  args.first = 0;
  args.count = numPoints;
  GpuArguments packed(binding, args, nullptr, nullptr, {&descriptor});
  invoke(packed.data());
  outside = descriptor.outside;
  return same && unit_test::bitwiseEqual(interpreted.data(), emitted.data(), interpreted.size());
}

/// The samples of a grid by the instructions of the written-out lookups, in the interpreter,
/// against those of datafield::sampleBatch() at the same points -- inputs of the program, so that
/// nothing is computed before the lookup: whether they are the same bits. In a program of binary32,
/// the coordinates arrive in binary32, and the samples leave in it.
bool samplesAgree(df::Interpolation interpolation,
                  const std::vector<unsigned>& num,
                  df::ElementType element,
                  ComputeType type,
                  const std::vector<double>& x,
                  const std::vector<double>& y) {
  df::GridDesc desc;
  desc.path = "grid.nc";
  desc.interpolation = interpolation;
  const std::size_t dimensions = num.size();
  const std::size_t numPoints = x.size();
  // the coordinates: x, y and mixes of both
  std::vector<std::vector<double>> columns(dimensions, std::vector<double>(numPoints));
  for (std::size_t p = 0; p < numPoints; ++p) {
    const double mixes[4] = {x[p], y[p], 0.5 * x[p] + y[p], y[p] - 0.25 * x[p]};
    for (std::size_t d = 0; d < dimensions; ++d) {
      columns[d][p] = mixes[d];
    }
  }
  Program program;
  program.setComputeType(type);
  const GridId gridId = program.internGrid(desc);
  std::vector<NodeId> coordinates;
  const std::vector<std::string> names = {"c0", "c1", "c2", "c3"};
  for (std::size_t d = 0; d < dimensions; ++d) {
    program.addInput(names[d], reader::scripting::DataType::F64);
    coordinates.push_back(program.arena().field(names[d]));
  }
  program.addOutput(
      "u", reader::scripting::DataType::F64, program.arena().lookup(gridId, 0, coordinates));
  program.addOutput(
      "v", reader::scripting::DataType::F64, program.arena().lookup(gridId, 1, coordinates));

  std::vector<double> sampled(2 * numPoints, -1.0);
  DataTable table(numPoints);
  for (std::size_t d = 0; d < dimensions; ++d) {
    table.bindViewConst<double>(names[d], Direction::In, columns[d].data());
  }
  table.bindView<double>("u", Direction::Out, sampled.data());
  table.bindView<double>("v", Direction::Out, sampled.data() + numPoints);
  Binding binding = Binding::bind(program, table);
  df::GridStore store;
  const std::size_t id = store.intern(desc);
  store.injectForTesting(id, std::make_unique<ArrayGrid>(num, element, interpolation));
  BackendOptions options;
  options.quiet = true;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);

  // the reference, at the coordinates as the program reads them
  const bool narrow = type == ComputeType::F32;
  std::vector<std::vector<double>> read = columns;
  std::vector<const double*> axes;
  for (auto& column : read) {
    if (narrow) {
      for (double& value : column) {
        value = static_cast<double>(static_cast<float>(value));
      }
    }
    axes.push_back(column.data());
  }
  const ArrayGrid& grid = dynamic_cast<const ArrayGrid&>(store.get(id));
  std::vector<double> expected(2 * numPoints);
  if (narrow) {
    std::vector<float> samples(2 * numPoints);
    float* const outs[2] = {samples.data(), samples.data() + numPoints};
    grid.sampleBatch(axes.data(), numPoints, outs);
    for (std::size_t i = 0; i < samples.size(); ++i) {
      expected[i] = static_cast<double>(samples[i]);
    }
  } else {
    double* const outs[2] = {expected.data(), expected.data() + numPoints};
    grid.sampleBatch(axes.data(), numPoints, outs);
  }
  return unit_test::bitwiseEqual(expected.data(), sampled.data(), expected.size());
}

/// Points around a grid of the geometry of ArrayGrid: inside, on its nodes and edges, beyond them,
/// and NaN and the infinities.
void lookupPoints(std::vector<double>& x, std::vector<double>& y) {
  constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
  constexpr double Inf = std::numeric_limits<double>::infinity();
  x = {-1.0, -0.75, -0.3, 0.0,  0.123, 0.5, 0.77, 1.0, 1.25, 1.5,  2.0,
       2.2,  2.999, 3.0,  -5.0, 10.0,  0.4, 1.9,  Nan, Inf,  -Inf, 0.6};
  y.resize(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    y[i] = -0.7 + 0.29 * static_cast<double>(i);
  }
  y[3] = Nan;
  y[5] = Inf;
  y[7] = -0.7;
}

/// Points inside a grid of the geometry of ArrayGrid in two dimensions.
void insidePoints(std::vector<double>& x, std::vector<double>& y) {
  x = {-0.9, -0.2, 0.0, 0.5, 1.1};
  y = {-0.6, 0.0, 0.3, 0.9, 1.4};
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

  TEST_CASE("a written-out lookup samples as Interpolation.cpp does, bit for bit") {
    // the instructions of expandLookups() against datafield::sampleBatch(), the reference
    std::vector<double> x;
    std::vector<double> y;
    lookupPoints(x, y);
    const std::vector<std::vector<unsigned>> shapes = {
        {7}, {7, 5}, {7, 3, 6}, {6, 5, 2, 4}, {1, 9}, {2, 3}};
    for (const auto interpolation :
         {df::Interpolation::Nearest, df::Interpolation::Linear, df::Interpolation::Cubic}) {
      for (const auto& shape : shapes) {
        for (const auto element : {df::ElementType::Double, df::ElementType::Float}) {
          for (const auto type : {ComputeType::F64, ComputeType::F32}) {
            CAPTURE(static_cast<int>(interpolation));
            CAPTURE(shape.size());
            CAPTURE(shape[0]);
            CAPTURE(static_cast<int>(element));
            CAPTURE(static_cast<int>(type));
            CHECK(samplesAgree(interpolation, shape, element, type, x, y));
          }
        }
      }
    }
  }

  TEST_CASE("the emitted device kernel samples a grid bit for bit as the interpreter does") {
    // Every scheme in one to four dimensions, of binary32 and binary64 elements, in programs of
    // either compute type: the locate branches (the edges, the fallback of the cubic stencil near
    // them and on axes shorter than it, the clamp of NaN and the infinities), the gather and the
    // order of the sums of the tensor product. Called through GpuArguments, which covers the
    // descriptor of the grid in the argument block.
    std::vector<double> x;
    std::vector<double> y;
    lookupPoints(x, y);
    const std::vector<std::vector<unsigned>> shapes = {
        {7}, {7, 5}, {7, 3, 6}, {6, 5, 2, 4}, {1, 9}, {2, 3}};
    for (const auto interpolation :
         {df::Interpolation::Nearest, df::Interpolation::Linear, df::Interpolation::Cubic}) {
      for (const auto& shape : shapes) {
        for (const auto element : {df::ElementType::Double, df::ElementType::Float}) {
          for (const auto type : {ComputeType::F64, ComputeType::F32}) {
            CAPTURE(static_cast<int>(interpolation));
            CAPTURE(shape.size());
            CAPTURE(shape[0]);
            CAPTURE(static_cast<int>(element));
            CAPTURE(static_cast<int>(type));
            if (interpolation == df::Interpolation::Cubic && shape.size() >= 3) {
              // thousands of instructions per lookup, which take a compiler a while: the
              // instructions are those of fewer dimensions, and the interpreter checks them
              // against datafield::sampleBatch() (above)
              continue;
            }
            bool ran = false;
            std::uint64_t outside = 0;
            // the compiled CPU kernel for a quarter of them, which the time of its compiler buys,
            // alternately in both of its shapes
            const bool withCpu = element == df::ElementType::Double && type == ComputeType::F64;
            const KernelShape kernelShape =
                shape.size() % 2 == 0 ? KernelShape::Point : KernelShape::Tile;
            const bool same = lookupsAgree(interpolation,
                                           df::BoundaryMode::Clamp,
                                           shape,
                                           element,
                                           type,
                                           x,
                                           y,
                                           kernelShape,
                                           ran,
                                           outside,
                                           withCpu);
            if (!ran) {
              WARN_MESSAGE(ran, "no usable C++ compiler; the device code generator was not run");
              return;
            }
            CHECK(same);
            // clamping is no error for a grid that clamps
            CHECK(outside == 0);
          }
        }
      }
    }
  }

  TEST_CASE("the emitted device kernel marks a lookup outside a grid that fails there") {
    for (const bool beyond : {false, true}) {
      CAPTURE(beyond);
      std::vector<double> x;
      std::vector<double> y;
      insidePoints(x, y);
      if (beyond) {
        x.push_back(4.0);
        y.push_back(0.0);
      }
      bool ran = false;
      std::uint64_t outside = 0;
      const bool same = lookupsAgree(df::Interpolation::Linear,
                                     df::BoundaryMode::Fail,
                                     {7, 5},
                                     df::ElementType::Double,
                                     ComputeType::F64,
                                     x,
                                     y,
                                     KernelShape::Point,
                                     ran,
                                     outside);
      if (!ran) {
        WARN_MESSAGE(ran, "no usable C++ compiler; the device code generator was not run");
        return;
      }
      CHECK(same);
      CHECK(outside == (beyond ? 1U : 0U));
    }
  }

  TEST_CASE("the device source of a lookup is valid for every target") {
    df::GridDesc desc;
    desc.path = "grid.nc";
    desc.interpolation = df::Interpolation::Cubic;
    const Program program = lookupProgram(desc, 3, ComputeType::F32);
    std::vector<double> x(2);
    std::vector<double> y(2);
    std::vector<double> out(2);
    DataTable table(2);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    for (const auto& output : program.outputs()) {
      table.bindView<double>(output.name, Direction::Out, out.data());
    }
    const Binding binding = Binding::bind(program, table);
    GpuLayout layout = gpuLayoutOf(binding);
    layout.grids = {df::ElementType::Float};
    const LoweredProgram lowered = lower(program);

    // each product and sum rounded on its own, on every target
    const std::string cuda = emitGpuSource(program, lowered, layout, GpuTarget::Cuda);
    CHECK(cuda.find("__dmul_rn") != std::string::npos);
    CHECK(cuda.find("SeissolExprGrid* grid0;") != std::string::npos);
    CHECK(cuda.find("(const float*)(g->data)") != std::string::npos);
    const std::string opencl = emitGpuSource(program, lowered, layout, GpuTarget::OpenCl);
    CHECK(opencl.find("#pragma OPENCL FP_CONTRACT OFF") != std::string::npos);
    CHECK(opencl.find("__global SeissolExprGrid* grid0,") != std::string::npos);
    CHECK(opencl.find("llround") == std::string::npos);

    // the element type of the array is part of the kernel
    GpuLayout wide = layout;
    wide.grids = {df::ElementType::Double};
    CHECK(wide.fingerprint() != layout.fingerprint());

    // a program without grids keeps its source
    const Program plain = compileSderivModule("out def u = x * y\n");
    DataTable plainTable(2);
    plainTable.bindViewConst<double>("x", Direction::In, x.data());
    plainTable.bindViewConst<double>("y", Direction::In, y.data());
    plainTable.bindView<double>("u", Direction::Out, out.data());
    const Binding plainBinding = Binding::bind(plain, plainTable);
    CHECK(emitGpuSource(plain, lower(plain), gpuLayoutOf(plainBinding), GpuTarget::Cuda)
              .find("SeissolExprGrid") == std::string::npos);

    // where the toolkit is there, its compiler takes the source as the driver gives it
    bool ran = false;
    std::string log;
    CHECK_MESSAGE(nvrtcAccepts(cuda, ran, log), log);
  }

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

    SUBCASE("a grid lookup is none: the kernel samples the grid itself") {
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
      CHECK(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);
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
        {"binary32 in a program of binary64",
         "def a = f32(x) * f32(y) - sqrt(abs(f32(x))) + sin(f32(y)) + exp(f32(x) / 1e3)\n"
         "out def u = real(a) + x\nout def v = f32(y) * 2.0\n"},
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
    // a column per cell carries its divisor, index and point stride along
    CHECK(packed.fieldCount() == 3 * 3 + 3 + 4);
  }

  TEST_CASE("a column per point of a cell reaches the device kernel, to read and to write") {
    // the layout of a tensor per face: three points per cell, the cells eight values apart, the
    // points two apart; read through a cell index, written in the order of the cells
    const Program program = compileSderivModule("out def u = x * j + 1.0\n");
    constexpr std::size_t PointsPerCell = 3;
    constexpr std::size_t Cells = 3;
    constexpr std::size_t NumPoints = Cells * PointsPerCell;
    constexpr std::size_t CellStride = 8;
    constexpr std::size_t PointStride = 2;
    std::vector<float> x(4 * CellStride);
    for (std::size_t i = 0; i < x.size(); ++i) {
      x[i] = 0.25F * static_cast<float>(i);
    }
    const std::vector<double> j = {1.0, -2.0, 3.0, -4.0};
    const std::vector<std::uint32_t> cellIndex = {3, 0, 2};
    std::vector<double> interpreted(Cells * CellStride, -1.0);
    std::vector<double> emitted(Cells * CellStride, -1.0);

    DataTable table(NumPoints);
    table.bindCellPointViewConst<float>(
        "x", x.data(), PointsPerCell, CellStride, PointStride, 1, cellIndex.data());
    table.bindCellView<double>("j", j.data(), PointsPerCell, 1, 0, cellIndex.data());
    table.bindCellPointView<double>(
        "u", Direction::Out, interpreted.data(), PointsPerCell, CellStride, PointStride, 3);
    Binding binding = Binding::bind(program, table);
    REQUIRE(gpuRejection(program, lower(program), binding, nullptr) == GpuRejection::None);
    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);
    for (std::size_t p = 0; p < NumPoints; ++p) {
      const std::size_t cell = p / PointsPerCell;
      const std::size_t within = (p % PointsPerCell) * PointStride;
      const double xp = x[cellIndex[cell] * CellStride + within + 1];
      CHECK(interpreted[cell * CellStride + within + 3] == xp * j[cellIndex[cell]] + 1.0);
    }

    const GpuLayout layout = gpuLayoutOf(binding);
    REQUIRE(layout.inputPerCell.size() == 2);
    CHECK(layout.inputPerCell[0]);
    CHECK(layout.inputPerCell[1]);
    REQUIRE(layout.outputPerCell.size() == 1);
    CHECK(layout.outputPerCell[0]);
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
    // the same elements written, the same values; the others untouched
    CHECK(unit_test::bitwiseEqual(interpreted.data(), emitted.data(), interpreted.size()));
    // every column carries a divisor, an index and a point stride along
    CHECK(packed.fieldCount() == 3 * 3 + 3 * 3 + 4);
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
    // Columns per point, per cell and per point of a cell, a uniform input, a state, a
    // contraction, in both compute types: every form the argument block and the accessors take.
    for (const auto computeType : {ComputeType::F64, ComputeType::F32}) {
      Program program =
          compileSderivModule("state s = 0.0\nout def s = s + v * t\nout def u = 2.0 * v - "
                              "j * x\nout def w = y * t\n");
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
      std::vector<double> y(2 * NumPoints, 1.0);
      std::vector<double> w(2 * NumPoints);
      const std::vector<std::uint32_t> cellIndex = {1, 0};
      DataTable table(NumPoints);
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindCellView<float>("j", j.data(), 2);
      table.bindCellPointViewConst<double>("y", y.data(), 2, 4, 2, 1, cellIndex.data());
      table.bindCellPointView<double>("w", Direction::Out, w.data(), 2, 4, 2);
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
