// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/RtcCpu.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Codegen.h"
#include "Expr/Interp.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <fcntl.h>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <sys/mman.h>
#include <sys/wait.h>
#include <tuple>
#include <unistd.h>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

using reader::scripting::DataTable;

// The CPU kernel works on a gathered tile, so its addressing is tile-relative.
// The arithmetic in between is emitted by codegen::emitStageBody and is the
// same text the GPU backends get.
std::string cpuLoadInput(std::int32_t index) {
  return "inputTile[" + std::to_string(index) + "ul * count + l]";
}
std::string cpuLoadPersistent(std::int32_t slot) {
  return "persistent[" + std::to_string(slot) + "ul * numPoints + first + l]";
}
std::string cpuStoreOutput(std::int32_t index) {
  return "outputTile[" + std::to_string(index) + "ul * count + l]";
}
std::string cpuStorePersistent(std::int32_t slot) {
  return "persistent[" + std::to_string(slot) + "ul * numPoints + first + l]";
}

void emitStage(std::ostringstream& out,
               const char* name,
               const StageCode& stage,
               const std::vector<std::int32_t>& operands,
               const std::string& computeType) {
  out << "extern \"C\" void " << name << "(const " << computeType << "* __restrict inputTile,\n"
      << "                                 " << computeType << "* __restrict outputTile,\n"
      << "                                 " << computeType << "* __restrict persistent,\n"
      << "                                 unsigned long numPoints,\n"
      << "                                 unsigned long first,\n"
      << "                                 unsigned long count) {\n";
  if (stage.code.empty() && stage.outputs.empty() && stage.persistent.empty()) {
    out << "  (void)inputTile; (void)outputTile; (void)persistent;\n"
        << "  (void)numPoints; (void)first; (void)count;\n}\n\n";
    return;
  }

  codegen::StageAddressing addressing;
  addressing.loadInput = cpuLoadInput;
  addressing.loadPersistent = cpuLoadPersistent;
  addressing.storeOutput = cpuStoreOutput;
  addressing.storePersistent = cpuStorePersistent;

  out << "  for (unsigned long l = 0; l < count; ++l) {\n";
  codegen::emitStageBody(
      out, stage, operands, computeType, codegen::MathStyle::Namespaced, addressing, "    ");
  out << "  }\n}\n\n";
}

// --- compilation ------------------------------------------------------------

struct Artifact {
  void* handle{nullptr};
  void* precompute{nullptr};
  void* run{nullptr};
  /// Held open for the artifact's whole life, and NOT as a leak.
  ///
  /// dlopen keys its "already loaded" cache on the path it was given, and the
  /// path here is /proc/self/fd/N. Closing the descriptor frees N for the next
  /// memfd, so the second kernel would be opened under the same path as the
  /// first -- and dlopen hands back the FIRST library. That is silent: the
  /// second program then computes the first one's expression. Keeping the
  /// descriptor keeps the path unique, which is the property dlopen is relying
  /// on. Artifacts live in the cache for the process's lifetime anyway.
  int fd{-1};
};

const char* compilerCommand() {
  if (const char* env = std::getenv("SEISSOL_EXPR_CXX"); env != nullptr) {
    return env;
  }
  return "c++";
}

// Compile `source` and return the loaded artifact, or an empty one.
//
// The artifact lives in an anonymous memfd and is dlopen'd through
// /proc/self/fd, so it never has a filesystem name. That removes three problems
// at once: no collision between ranks on a node writing the same hashed path,
// no O_EXCL-and-rename dance, and no shared filesystem being hammered by every
// rank at once. The compiler's own intermediates still go to TMPDIR; they are
// per-process unique and short-lived, and TMPDIR=/dev/shm keeps even those in
// RAM.
Artifact compileAndLoad(const std::string& source, const std::string& flags) {
  Artifact artifact;

  const int object = memfd_create("seissol-expr-kernel", 0);
  if (object < 0) {
    logWarning() << "expr: memfd_create failed; not using the compiled CPU backend.";
    return artifact;
  }
  std::string objectPath = "/proc/self/fd/" + std::to_string(object);

  const int input = memfd_create("seissol-expr-source", 0);
  if (input < 0) {
    close(object);
    return artifact;
  }
  if (write(input, source.data(), source.size()) != static_cast<ssize_t>(source.size())) {
    close(object);
    close(input);
    return artifact;
  }

  std::vector<std::string> argv{compilerCommand()};
  std::istringstream split(flags);
  for (std::string token; split >> token;) {
    argv.push_back(token);
  }
  for (const char* fixed : {"-shared", "-fPIC", "-x", "c++", "-", "-o"}) {
    argv.emplace_back(fixed);
  }
  argv.push_back(objectPath);

  std::vector<char*> raw;
  raw.reserve(argv.size() + 1);
  for (auto& arg : argv) {
    raw.push_back(arg.data());
  }
  raw.push_back(nullptr);

  const pid_t pid = fork();
  if (pid < 0) {
    close(object);
    close(input);
    return artifact;
  }
  if (pid == 0) {
    lseek(input, 0, SEEK_SET);
    dup2(input, STDIN_FILENO);
    execvp(raw[0], raw.data());
    _exit(127);
  }

  int status = 0;
  waitpid(pid, &status, 0);
  close(input);

  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    logWarning() << "expr: the CPU kernel did not compile (" << compilerCommand() << " exited with"
                 << WEXITSTATUS(status) << "); falling back to the interpreter.";
    close(object);
    return artifact;
  }

  artifact.handle = dlopen(objectPath.c_str(), RTLD_NOW | RTLD_LOCAL);
  artifact.fd = object; // see the note on Artifact::fd -- must outlive the dlopen
  if (artifact.handle == nullptr) {
    close(object);
    artifact.fd = -1;
    // Hardened kernels can make a memfd non-executable (vm.memfd_noexec), and
    // /proc may not be mounted. Both are configuration, not programming errors,
    // so this is a warning and the interpreter takes over.
    logWarning() << "expr: could not load the compiled CPU kernel (" << dlerror()
                 << "); falling back to the interpreter.";
    return artifact;
  }
  artifact.precompute = dlsym(artifact.handle, "seissol_expr_precompute");
  artifact.run = dlsym(artifact.handle, "seissol_expr_run");
  if (artifact.run == nullptr) {
    dlclose(artifact.handle);
    close(artifact.fd);
    artifact.handle = nullptr;
    artifact.fd = -1;
  }
  return artifact;
}

// --- cache ------------------------------------------------------------------

// Keyed on the program AND on everything that shapes the emitted code but is
// not part of the program: the lowering options, the compute type and the
// compiler flags. Program::fingerprint() alone would hand an -march=skylake
// artifact to a run configured for something else.
struct CacheKey {
  std::uint64_t program{0};
  std::uint64_t lowering{0};
  ComputeType type{ComputeType::F64};
  std::string flags;

  bool operator<(const CacheKey& other) const {
    return std::tie(program, lowering, type, flags) <
           std::tie(other.program, other.lowering, other.type, other.flags);
  }
};

std::mutex& cacheMutex() {
  static std::mutex mutex;
  return mutex;
}

std::map<CacheKey, Artifact>& cache() {
  static std::map<CacheKey, Artifact> artifacts;
  return artifacts;
}

// --- kernel -----------------------------------------------------------------

template <typename T>
using StageFn = void (*)(const T*, T*, T*, unsigned long, unsigned long, unsigned long);

template <typename T>
class RtcCpuKernel final : public Kernel {
  public:
  RtcCpuKernel(Binding& binding,
               LoweredProgram lowered,
               const Artifact& artifact,
               std::size_t tileSize)
      : binding_(&binding), lowered_(std::move(lowered)),
        precompute_(reinterpret_cast<StageFn<T>>(artifact.precompute)),
        run_(reinterpret_cast<StageFn<T>>(artifact.run)), tileSize_(tileSize),
        inputTile_(static_cast<std::size_t>(std::max<std::size_t>(1, binding.inputs().size())) *
                   tileSize),
        outputTile_(static_cast<std::size_t>(std::max<std::size_t>(1, binding.outputs().size())) *
                    tileSize),
        needsPrecompute_(lowered_.hasPrecompute()) {}

  void precompute(const DataTable& table) override {
    if (!needsPrecompute_ || precompute_ == nullptr) {
      return;
    }
    const BoundIo io(*binding_, &table, nullptr);
    sweep(precompute_, io, 0, binding_->numPoints());
    precomputed_ = true;
  }

  void run(const DataTable& table) override {
    guard();
    const BoundIo io(*binding_, &table, nullptr);
    // The dense path may thread; the element-wise one below may not.
    sweep(run_, io, 0, binding_->numPoints());
  }

  void run(const KernelArgs& args) override {
    if (!binding_->addressable()) {
      logError() << "expr: this program has a computed column and cannot be evaluated from raw "
                    "bases; call run(table).";
      return;
    }
    guard();
    const BoundIo io(*binding_, nullptr, &args);
    sweep(run_, io, args.first, args.first + args.count);
  }

  [[nodiscard]] BackendKind kind() const override { return BackendKind::RtcCpu; }

  private:
  /// Gather/scatter over either a bound table or per-call bases. One type
  /// rather than two TileIo subclasses because the compiled kernel is not
  /// virtual — there is nothing here for a vtable to buy.
  struct BoundIo {
    BoundIo(const Binding& binding, const DataTable* table, const KernelArgs* args)
        : binding(&binding), table(table), args(args) {}

    // An addressable binding reads through the StridedView even on the dense
    // path. Not an optimisation to be tuned away: measured, the accessor route
    // costs about 30 ns/point in std::function calls, which is MORE than the
    // whole compiled kernel, so leaving run(table) on it would have hidden the
    // entire benefit of compiling. The two routes read the same bytes by
    // construction -- the view is where the accessor's closure got its pointer.
    void gather(std::size_t first, std::size_t count, T* dst) const {
      if (args != nullptr) {
        binding->gatherFrom(args->inputs, args->inputCount, first, count, dst);
      } else if (binding->addressable()) {
        binding->gatherFrom(nullptr, 0, first, count, dst);
      } else {
        binding->gather(*table, first, count, dst);
      }
    }
    void scatter(std::size_t first, std::size_t count, const T* src) const {
      if (args != nullptr) {
        binding->scatterTo(args->outputs, args->outputCount, first, count, src);
      } else if (binding->addressable()) {
        binding->scatterTo(nullptr, 0, first, count, src);
      } else {
        binding->scatter(*table, first, count, src);
      }
    }

    const Binding* binding;
    const DataTable* table;
    const KernelArgs* args;
  };

  void guard() const {
    if (needsPrecompute_ && !precomputed_) {
      logError() << "expr: the kernel has a precompute stage that was never run; call "
                    "Kernel::precompute() from prepare().";
    }
  }

  void sweep(StageFn<T> fn, const BoundIo& io, std::size_t begin, std::size_t end) {
    if (fn == nullptr) {
      return;
    }
    for (std::size_t first = begin; first < end; first += tileSize_) {
      const std::size_t count = std::min(tileSize_, end - first);
      io.gather(first, count, inputTile_.data());
      fn(inputTile_.data(), outputTile_.data(), persistent(), binding_->numPoints(), first, count);
      io.scatter(first, count, outputTile_.data());
    }
  }

  T* persistent();

  Binding* binding_;
  LoweredProgram lowered_;
  StageFn<T> precompute_{nullptr};
  StageFn<T> run_{nullptr};
  std::size_t tileSize_{0};
  std::vector<T> inputTile_;
  std::vector<T> outputTile_;
  bool needsPrecompute_{false};
  bool precomputed_{false};
};

template <>
double* RtcCpuKernel<double>::persistent() {
  return binding_->persistentF64();
}
template <>
float* RtcCpuKernel<float>::persistent() {
  return binding_->persistentF32();
}

std::string defaultFlags(const BackendOptions& options) {
  if (const char* env = std::getenv("SEISSOL_EXPR_CXXFLAGS"); env != nullptr) {
    return env;
  }
  // -ffp-contract=off is not optional and not a tuning knob. Without it GCC
  // fuses a*b+c into an FMA, which rounds once where the interpreter rounds
  // twice; measured, that alone moves ordinary arithmetic by ~1e-6 and breaks
  // the acceptance criterion. It costs nothing measurable.
  std::string flags = "-O3 -ffp-contract=off -fno-math-errno";
  flags += options.arch.empty() || options.arch == "native" ? " -march=native"
                                                            : " -march=" + options.arch;
  return flags;
}

} // namespace

bool cpuCompilable(const LoweredProgram& lowered) { return !codegen::containsLookup(lowered); }

std::string emitCpuSource(const Program& program, const LoweredProgram& lowered) {
  const std::string computeType = program.computeType() == ComputeType::F32 ? "float" : "double";

  std::ostringstream out;
  out << "// Generated by seissol::expr for program fingerprint 0x" << std::hex
      << program.fingerprint() << std::dec << ".\n"
      << "// The arithmetic below is stringified from SEISSOL_EXPR_PW_LIST, the same\n"
      << "// table the interpreter evaluates, so the two cannot disagree.\n"
      << "#include <cmath>\n\n";
  emitStage(out, "seissol_expr_precompute", lowered.precompute(), lowered.operands(), computeType);
  emitStage(out, "seissol_expr_run", lowered.run(), lowered.operands(), computeType);
  return out.str();
}

std::size_t rtcCpuCacheSize() {
  const std::lock_guard<std::mutex> lock(cacheMutex());
  return cache().size();
}

std::unique_ptr<Kernel> makeRtcCpuKernel(const Program& program,
                                         Binding& binding,
                                         LoweredProgram lowered,
                                         const BackendOptions& options) {
  if (!cpuCompilable(lowered)) {
    logWarning() << "expr: this program samples a data grid, which the compiled CPU backend does "
                    "not support yet; using the interpreter.";
    return nullptr;
  }

  const std::string flags = defaultFlags(options);
  const CacheKey key{
      program.fingerprint(), options.lowering.fingerprint(), program.computeType(), flags};

  Artifact artifact;
  {
    // Held across the compile on purpose. This runs from prepare(), i.e. at
    // init, where serialising a handful of 40 ms compiles costs nothing and
    // two threads compiling the same program would be pure waste. If this path
    // is ever reached during the timestep loop, that is the bug to fix rather
    // than the lock.
    const std::lock_guard<std::mutex> lock(cacheMutex());
    auto found = cache().find(key);
    if (found == cache().end()) {
      const std::string source = emitCpuSource(program, lowered);
      found = cache().emplace(key, compileAndLoad(source, flags)).first;
    }
    artifact = found->second;
  }
  if (artifact.handle == nullptr) {
    return nullptr;
  }

  binding.allocatePersistent(program, lowered.persistentSlotCount());
  const std::size_t tileSize =
      options.tileSize != 0
          ? options.tileSize
          : chooseTileSize(lowered.peakSlotCount(), program.computeType(), DefaultTileBudgetBytes);

  logInfo() << "expr: compiled CPU kernel --" << lowered.summary().c_str() << "-- with" << flags;

  if (program.computeType() == ComputeType::F32) {
    return std::make_unique<RtcCpuKernel<float>>(binding, std::move(lowered), artifact, tileSize);
  }
  return std::make_unique<RtcCpuKernel<double>>(binding, std::move(lowered), artifact, tileSize);
}

} // namespace seissol::expr
