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
#include "Expr/Cost.h"
#include "Expr/Expand.h"
#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
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
#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <system_error>
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
std::string cpuLoadState(std::int32_t state) {
  return "stateTile[" + std::to_string(state) + "ul * count + l]";
}
std::string cpuLoadPersistent(std::int32_t slot) {
  return "persistent[" + std::to_string(slot) + "ul * numPoints + first + l]";
}
std::string cpuStoreOutput(std::int32_t index, const std::string& value) {
  return "outputTile[" + std::to_string(index) + "ul * count + l] = " + value;
}
std::string cpuStoreState(std::int32_t state, const std::string& value) {
  return "stateTile[" + std::to_string(state) + "ul * count + l] = " + value;
}
std::string cpuStorePersistent(std::int32_t slot, const std::string& value) {
  return "persistent[" + std::to_string(slot) + "ul * numPoints + first + l] = " + value;
}
std::string cpuMatrix(std::int32_t matrix) { return "matrices[" + std::to_string(matrix) + "]"; }
std::string cpuBlockBase(std::int32_t block) {
  return "blocks[" + std::to_string(block) + "].base";
}
std::string cpuCellStride(std::int32_t block) {
  return "blocks[" + std::to_string(block) + "].cellStride";
}
std::string cpuModeStride(std::int32_t block) {
  return "blocks[" + std::to_string(block) + "].modeStride";
}
std::string cpuCellIndex(std::int32_t block) {
  return "blocks[" + std::to_string(block) + "].cellIndex";
}
std::string cpuGridParam(GridId grid, GridField field, std::int32_t axis) {
  const std::string g = "grids[" + std::to_string(grid) + "].";
  const std::string a = "[" + std::to_string(axis) + "]";
  switch (field) {
  case GridField::Min:
    return g + "min" + a;
  case GridField::DeltaInv:
    return g + "deltaInv" + a;
  case GridField::Num:
    return g + "num" + a;
  case GridField::Stride:
    return g + "stride" + a;
  case GridField::ComponentStride:
    return g + "componentStride";
  case GridField::Before:
    return g + "before" + a;
  case GridField::After:
    return g + "after" + a;
  }
  return "0";
}
std::string cpuGridLoad(GridId grid, std::int32_t component, const std::string& offset) {
  const std::string g = std::to_string(grid);
  return "seissol_expr_load" + g + "(&grids[" + g + "], " + offset + ", " +
         std::to_string(component) + ")";
}
std::string cpuStoreCheck(GridId grid, CheckKind kind, const std::string& value) {
  return "if (" + value + ") { grids[" + std::to_string(grid) + "]." +
         (kind == CheckKind::Window ? "stale" : "outside") + " = 1; }";
}

/// The element types of the arrays of the grids, by GridId: baked into the source, hence part of
/// the cache key.
std::uint64_t gridFingerprint(const std::vector<reader::datafield::ElementType>& elements) {
  std::uint64_t hash = 0xcbf29ce484222325ULL;
  for (const auto type : elements) {
    hash ^= static_cast<std::uint64_t>(type) + 1;
    hash *= 0x100000001b3ULL;
  }
  return hash;
}

/// The host-side image of the emitted SeissolExprBlock; the two have to agree field by field.
struct BlockDescriptor {
  const void* base;
  unsigned long cellStride;
  unsigned long modeStride;
  const std::uint32_t* cellIndex;
};

/// The element types the matrices and blocks are bound with: baked into the source, hence part of
/// the cache key.
struct ContractLayout {
  std::vector<reader::scripting::DataType> matrices;
  std::vector<reader::scripting::DataType> blocks;

  [[nodiscard]] std::uint64_t fingerprint() const {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    const auto mix = [&hash](std::uint64_t value) {
      hash ^= value;
      hash *= 0x100000001b3ULL;
    };
    mix(matrices.size());
    for (const auto type : matrices) {
      mix(static_cast<std::uint64_t>(type));
    }
    mix(blocks.size());
    for (const auto type : blocks) {
      mix(static_cast<std::uint64_t>(type));
    }
    return hash;
  }
};

ContractLayout contractLayoutOf(const Binding& binding) {
  ContractLayout layout;
  for (const auto& matrix : binding.matrices()) {
    layout.matrices.push_back(matrix.type);
  }
  for (const auto& block : binding.blocks()) {
    layout.blocks.push_back(block.type);
  }
  return layout;
}

void emitStage(std::ostringstream& out,
               const char* name,
               const StageCode& stage,
               const std::vector<std::int32_t>& operands,
               std::int32_t stateSlots,
               ComputeType dataType,
               const codegen::ContractAddressing& contract,
               KernelShape shape) {
  // the tiles in the type the data moves in; the instructions compute in theirs
  const std::string data = codegen::typeName(dataType);
  out << "extern \"C\" void " << name << "(const " << data << "* __restrict inputTile,\n"
      << "                                 " << data << "* __restrict outputTile,\n"
      << "                                 " << data << "* __restrict stateTile,\n"
      << "                                 " << data << "* __restrict persistent,\n"
      << "                                 unsigned long numPoints,\n"
      << "                                 unsigned long first,\n"
      << "                                 unsigned long count,\n"
      << "                                 const void* const* matrices,\n"
      << "                                 const SeissolExprBlock* blocks,\n"
      << "                                 const unsigned long* pointIndex,\n"
      << "                                 SeissolExprGrid* grids) {\n"
      << "  (void)matrices; (void)blocks; (void)pointIndex; (void)stateTile; (void)grids;\n";
  if (stage.code.empty() && stage.outputs.empty() && stage.persistent.empty() &&
      stage.checks.empty()) {
    out << "  (void)inputTile; (void)outputTile; (void)persistent;\n"
        << "  (void)numPoints; (void)first; (void)count;\n}\n\n";
    return;
  }

  // The states come in a tile of their own, gathered and scattered around the call like the
  // inputs and outputs; the hoisted values are read where they lie.
  codegen::StageAddressing addressing;
  addressing.loadInput = cpuLoadInput;
  addressing.stateSlots = stateSlots;
  addressing.loadState = cpuLoadState;
  addressing.loadPersistent = cpuLoadPersistent;
  addressing.storeOutput = cpuStoreOutput;
  addressing.storeState = cpuStoreState;
  addressing.storePersistent = cpuStorePersistent;
  addressing.contract = &contract;
  addressing.gridParam = cpuGridParam;
  addressing.gridLoad = cpuGridLoad;
  addressing.storeCheck = cpuStoreCheck;

  if (tiledStage(stage, operands, shape)) {
    // four cache lines of a transient per block: 64 lanes in binary32, 32 where one is in binary64
    const bool narrow = std::all_of(stage.slotTypes.begin(),
                                    stage.slotTypes.end(),
                                    [](ComputeType type) { return type == ComputeType::F32; });
    const std::size_t block = narrow ? 64 : 32;
    codegen::emitStageBodyTiled(
        out, stage, operands, dataType, codegen::MathStyle::Namespaced, addressing, block, "  ");
    out << "}\n\n";
    return;
  }
  out << "  for (unsigned long l = 0; l < count; ++l) {\n";
  codegen::emitStageBody(
      out, stage, operands, dataType, codegen::MathStyle::Namespaced, addressing, "    ");
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

/// Loads the shared object in the memfd `object` (see Artifact::fd); takes the descriptor.
Artifact loadArtifact(int object) {
  Artifact artifact;
  const std::string objectPath = "/proc/self/fd/" + std::to_string(object);
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

// Compile `source` into a memfd and return its descriptor, or -1.
//
// The artifact lives in an anonymous memfd and is dlopen'd through
// /proc/self/fd, so it never has a filesystem name. That removes three problems
// at once: no collision between ranks on a node writing the same hashed path,
// no O_EXCL-and-rename dance, and no shared filesystem being hammered by every
// rank at once. The compiler's own intermediates still go to TMPDIR; they are
// per-process unique and short-lived, and TMPDIR=/dev/shm keeps even those in
// RAM. (The cache across runs below keeps a copy under a name, written once.)
int compile(const std::string& source, const std::string& flags) {
  const int object = memfd_create("seissol-expr-kernel", 0);
  if (object < 0) {
    logWarning() << "expr: memfd_create failed; not using the compiled CPU backend.";
    return -1;
  }
  const std::string objectPath = "/proc/self/fd/" + std::to_string(object);

  const int input = memfd_create("seissol-expr-source", 0);
  if (input < 0) {
    close(object);
    return -1;
  }
  if (write(input, source.data(), source.size()) != static_cast<ssize_t>(source.size())) {
    close(object);
    close(input);
    return -1;
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
    return -1;
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
    return -1;
  }
  return object;
}

// --- cache across runs ------------------------------------------------------
//
// The host compiler takes 0.1 to 0.3 s per program and process, a run with a few programs a second
// or more, in every rank. The compiled kernels are therefore kept in a directory across runs:
// SEISSOL_EXPR_CACHE names it ("off" for none), by default seissol/expr in $XDG_CACHE_HOME or
// ~/.cache. An entry is the shared object (<name>.so) and the text it was compiled from, with
// everything else that shapes it (<name>.key): the flags, the compiler (its path, size and time of
// modification), and for -march=native the processor. The name is a hash of that text, and a hit
// compares the text in full, so a collision costs a compile and never a wrong kernel. Entries are
// written to a temporary file and renamed, so that ranks that compile the same program at once
// leave a whole entry behind; a reader that finds half of one compiles.

/// The directory of the cache; empty without one.
const std::string& cacheDirectory() {
  static const std::string Directory = []() -> std::string {
    std::string directory;
    if (const char* env = std::getenv("SEISSOL_EXPR_CACHE"); env != nullptr) {
      const std::string value(env);
      if (value.empty() || value == "off" || value == "0" || value == "none") {
        return {};
      }
      directory = value;
    } else if (const char* xdg = std::getenv("XDG_CACHE_HOME"); xdg != nullptr && *xdg != '\0') {
      directory = std::string(xdg) + "/seissol/expr";
    } else if (const char* home = std::getenv("HOME"); home != nullptr && *home != '\0') {
      directory = std::string(home) + "/.cache/seissol/expr";
    } else {
      return {};
    }
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error || !std::filesystem::is_directory(directory, error)) {
      logWarning() << "expr: no cache of compiled kernels in" << directory << "(" << error.message()
                   << "); every run compiles them.";
      return {};
    }
    return directory;
  }();
  return Directory;
}

/// The compiler as the cache tells it apart: where it is, its size, and when it changed.
const std::string& compilerIdentity() {
  static const std::string Identity = [] {
    std::string command = compilerCommand();
    std::filesystem::path path;
    if (command.find('/') != std::string::npos) {
      path = command;
    } else if (const char* env = std::getenv("PATH"); env != nullptr) {
      std::istringstream directories(env);
      for (std::string directory; std::getline(directories, directory, ':');) {
        const auto candidate = std::filesystem::path(directory) / command;
        std::error_code error;
        if (std::filesystem::is_regular_file(candidate, error)) {
          path = candidate;
          break;
        }
      }
    }
    std::error_code error;
    const auto resolved = std::filesystem::canonical(path, error);
    if (error) {
      return command;
    }
    struct stat info{};
    if (stat(resolved.c_str(), &info) != 0) {
      return resolved.string();
    }
    return resolved.string() + " " + std::to_string(info.st_size) + " " +
           std::to_string(info.st_mtime);
  }();
  return Identity;
}

/// The processor, for code compiled with -march=native: vendor, model and flags.
const std::string& processorIdentity() {
  static const std::string Identity = [] {
    std::string identity;
    std::FILE* cpuinfo = std::fopen("/proc/cpuinfo", "r");
    if (cpuinfo == nullptr) {
      return identity;
    }
    std::array<char, 8192> line{};
    while (std::fgets(line.data(), static_cast<int>(line.size()), cpuinfo) != nullptr) {
      const std::string text(line.data());
      if (text.starts_with("vendor_id") || text.starts_with("model name") ||
          text.starts_with("flags")) {
        identity += text;
      }
      if (text == "\n" && !identity.empty()) {
        // the first processor is enough
        break;
      }
    }
    std::fclose(cpuinfo);
    return identity;
  }();
  return Identity;
}

/// What a kernel is compiled from, with all that shapes it besides.
std::string cacheKey(const std::string& source, const std::string& flags) {
  std::string key =
      "seissol-expr-cache 1\ncompiler " + compilerIdentity() + "\nflags " + flags + "\n";
  if (flags.find("-march=native") != std::string::npos) {
    key += "processor\n" + processorIdentity();
  }
  key += "source\n" + source;
  return key;
}

std::string cacheName(const std::string& key) {
  // two independent 64-bit hashes; a hit compares the key in full anyway
  std::uint64_t first = 0xcbf29ce484222325ULL;
  std::uint64_t second = 0x84222325cbf29ce4ULL;
  for (const unsigned char c : key) {
    first = (first ^ c) * 0x100000001b3ULL;
    second = (second ^ (c + 0x9eU)) * 0x100000001b3ULL;
  }
  std::array<char, 40> name{};
  std::snprintf(name.data(),
                name.size(),
                "%016llx%016llx",
                static_cast<unsigned long long>(first),
                static_cast<unsigned long long>(second));
  return name.data();
}

std::string readFile(const std::string& path) {
  std::string content;
  std::FILE* file = std::fopen(path.c_str(), "rb");
  if (file == nullptr) {
    return content;
  }
  std::array<char, 65536> buffer{};
  std::size_t read = 0;
  while ((read = std::fread(buffer.data(), 1, buffer.size(), file)) > 0) {
    content.append(buffer.data(), read);
  }
  std::fclose(file);
  return content;
}

/// Writes `content` to `path` through a temporary file and a rename.
void writeFileAtomically(const std::string& path, const std::string& content) {
  const std::string temporary =
      path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(std::rand());
  std::FILE* file = std::fopen(temporary.c_str(), "wb");
  if (file == nullptr) {
    return;
  }
  const bool written = std::fwrite(content.data(), 1, content.size(), file) == content.size();
  if (std::fclose(file) != 0 || !written) {
    std::remove(temporary.c_str());
    return;
  }
  if (std::rename(temporary.c_str(), path.c_str()) != 0) {
    std::remove(temporary.c_str());
  }
}

/// The contents of the memfd `object`.
std::string contentsOf(int object) {
  std::string content;
  std::array<char, 65536> buffer{};
  ssize_t read = 0;
  off_t offset = 0;
  while ((read = pread(object, buffer.data(), buffer.size(), offset)) > 0) {
    content.append(buffer.data(), static_cast<std::size_t>(read));
    offset += read;
  }
  return content;
}

/// The kernel of the entry of `key`, in a memfd; -1 if there is none.
int fromCache(const std::string& key) {
  const auto& directory = cacheDirectory();
  if (directory.empty()) {
    return -1;
  }
  const std::string base = directory + "/" + cacheName(key);
  if (readFile(base + ".key") != key) {
    return -1;
  }
  const std::string object = readFile(base + ".so");
  if (object.empty()) {
    return -1;
  }
  const int memfd = memfd_create("seissol-expr-kernel", 0);
  if (memfd < 0) {
    return -1;
  }
  if (write(memfd, object.data(), object.size()) != static_cast<ssize_t>(object.size())) {
    close(memfd);
    return -1;
  }
  return memfd;
}

void toCache(const std::string& key, int object) {
  const auto& directory = cacheDirectory();
  if (directory.empty()) {
    return;
  }
  const std::string base = directory + "/" + cacheName(key);
  // the key first: an entry without its object is a miss, as is one with another key
  writeFileAtomically(base + ".key", key);
  writeFileAtomically(base + ".so", contentsOf(object));
}

// Compile `source` and return the loaded artifact, or an empty one; from the cache across runs if
// it has the kernel, which `cached` tells.
Artifact compileAndLoad(const std::string& source, const std::string& flags, bool& cached) {
  const std::string key = cacheKey(source, flags);
  int object = fromCache(key);
  cached = object >= 0;
  if (object < 0) {
    object = compile(source, flags);
    if (object < 0) {
      return {};
    }
    toCache(key, object);
  }
  auto artifact = loadArtifact(object);
  if (artifact.handle == nullptr && cached) {
    // a damaged entry: compile instead
    cached = false;
    object = compile(source, flags);
    if (object >= 0) {
      toCache(key, object);
      artifact = loadArtifact(object);
    }
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
  std::uint64_t contraction{0};
  std::uint64_t grids{0};
  ComputeType type{ComputeType::F64};
  std::string flags;
  KernelShape shape{KernelShape::Auto};

  bool operator<(const CacheKey& other) const {
    return std::tie(program, lowering, contraction, grids, type, flags, shape) <
           std::tie(other.program,
                    other.lowering,
                    other.contraction,
                    other.grids,
                    other.type,
                    other.flags,
                    other.shape);
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
using StageFn = void (*)(const T*,
                         T*,
                         T*,
                         T*,
                         unsigned long,
                         unsigned long,
                         unsigned long,
                         const void* const*,
                         const BlockDescriptor*,
                         const unsigned long*,
                         codegen::GridDescriptor*);

template <typename T>
class RtcCpuKernel final : public Kernel {
  public:
  RtcCpuKernel(Binding& binding,
               LoweredProgram lowered,
               const Artifact& artifact,
               std::size_t tileSize,
               std::vector<const reader::datafield::Grid*> grids,
               std::vector<std::string> gridPaths)
      : binding_(&binding), lowered_(std::move(lowered)), grids_(std::move(grids)),
        gridPaths_(std::move(gridPaths)), descriptors_(grids_.size()),
        precompute_(reinterpret_cast<StageFn<T>>(artifact.precompute)),
        run_(reinterpret_cast<StageFn<T>>(artifact.run)), tileSize_(tileSize),
        inputTile_(static_cast<std::size_t>(std::max<std::size_t>(1, binding.inputs().size())) *
                   tileSize),
        outputTile_(static_cast<std::size_t>(std::max<std::size_t>(1, binding.outputs().size())) *
                    tileSize),
        stateTile_(static_cast<std::size_t>(std::max<std::size_t>(1, binding.states().size())) *
                   tileSize),
        needsPrecompute_(lowered_.hasPrecompute()) {}

  void precompute(const DataTable& table) override {
    if (!needsPrecompute_ || precompute_ == nullptr) {
      return;
    }
    const BoundIo io(*binding_, &table, nullptr);
    sweep(precompute_, io, nullptr, 0, binding_->numPoints());
    precomputed_ = true;
  }

  void run(const DataTable& table) override {
    guard();
    const BoundIo io(*binding_, &table, nullptr);
    // The dense path may thread; the element-wise one below may not.
    sweep(run_, io, nullptr, 0, binding_->numPoints());
  }

  void run(const KernelArgs& args) override {
    if (!binding_->hostAddressable()) {
      logError() << "expr: this program has a computed column and cannot be evaluated from raw "
                    "bases; call run(table).";
      return;
    }
    guard();
    const BoundIo io(*binding_, nullptr, &args);
    sweep(run_, io, &args, args.first, args.first + args.count);
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
      } else if (binding->hostAddressable()) {
        binding->gatherFrom(nullptr, 0, first, count, dst);
      } else {
        binding->gather(*table, first, count, dst);
      }
    }
    void scatter(std::size_t first, std::size_t count, const T* src) const {
      if (args != nullptr) {
        binding->scatterTo(args->outputs, args->outputCount, first, count, src);
      } else if (binding->hostAddressable()) {
        binding->scatterTo(nullptr, 0, first, count, src);
      } else {
        binding->scatter(*table, first, count, src);
      }
    }
    void gatherState(std::size_t first, std::size_t count, T* dst) const {
      if (args != nullptr) {
        binding->gatherState(args->states, args->stateCount, first, count, dst);
      } else {
        binding->gatherState(nullptr, 0, first, count, dst);
      }
    }
    void scatterState(std::size_t first, std::size_t count, const T* src) const {
      if (args != nullptr) {
        binding->scatterState(args->states, args->stateCount, first, count, src);
      } else {
        binding->scatterState(nullptr, 0, first, count, src);
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

  void sweep(StageFn<T> fn,
             const BoundIo& io,
             const KernelArgs* args,
             std::size_t begin,
             std::size_t end) {
    if (fn == nullptr) {
      return;
    }
    // The bases of the matrices and blocks for this call: the bound ones, or the ones the call
    // moves. A handful of entries, rebuilt per call.
    const auto& matrices = binding_->matrices();
    const auto& blocks = binding_->blocks();
    matrixBases_.resize(matrices.size());
    for (std::size_t i = 0; i < matrices.size(); ++i) {
      const bool moved = args != nullptr && i < args->matrixCount && args->matrices != nullptr &&
                         args->matrices[i] != nullptr;
      matrixBases_[i] = moved ? args->matrices[i] : matrices[i].base;
    }
    blockDescriptors_.resize(blocks.size());
    for (std::size_t i = 0; i < blocks.size(); ++i) {
      const bool moved = args != nullptr && i < args->blockCount && args->blocks != nullptr &&
                         args->blocks[i] != nullptr;
      blockDescriptors_[i] = BlockDescriptor{moved ? args->blocks[i] : blocks[i].base,
                                             blocks[i].cellStride,
                                             blocks[i].modeStride,
                                             blocks[i].cellIndex};
    }
    const auto& permutation = binding_->permutation();
    // only the run stage reads and writes the states
    const bool stateful = fn == run_ && !binding_->states().empty();
    // the grids where their windows lie now
    for (std::size_t g = 0; g < grids_.size(); ++g) {
      if (grids_[g] != nullptr) {
        const reader::datafield::ArrayView* view = grids_[g]->arrayView();
        if (view == nullptr) {
          logError() << "expr: the grid" << gridPaths_[g] << "keeps no array to sample.";
        }
        descriptors_[g] = codegen::gridDescriptor(*view, view->data);
      }
    }

    for (std::size_t first = begin; first < end; first += tileSize_) {
      const std::size_t count = std::min(tileSize_, end - first);
      io.gather(first, count, inputTile_.data());
      if (stateful) {
        io.gatherState(first, count, stateTile_.data());
      }
      fn(inputTile_.data(),
         outputTile_.data(),
         stateTile_.data(),
         persistent(),
         binding_->numPoints(),
         first,
         count,
         matrixBases_.data(),
         blockDescriptors_.data(),
         permutation.empty() ? nullptr : permutation.data() + first,
         descriptors_.data());
      io.scatter(first, count, outputTile_.data());
      if (stateful) {
        io.scatterState(first, count, stateTile_.data());
      }
    }
    for (std::size_t g = 0; g < grids_.size(); ++g) {
      if (descriptors_[g].stale != 0) {
        reportCheck(CheckKind::Window, gridPaths_[g]);
      }
      if (descriptors_[g].outside != 0) {
        reportCheck(CheckKind::Outside, gridPaths_[g]);
      }
    }
  }

  T* persistent();

  Binding* binding_;
  LoweredProgram lowered_;
  /// the grids it samples (by GridId; null for one it does not), their paths for the messages, and
  /// their descriptors, filled before every call
  std::vector<const reader::datafield::Grid*> grids_;
  std::vector<std::string> gridPaths_;
  std::vector<codegen::GridDescriptor> descriptors_;
  StageFn<T> precompute_{nullptr};
  StageFn<T> run_{nullptr};
  std::size_t tileSize_{0};
  std::vector<T> inputTile_;
  std::vector<T> outputTile_;
  std::vector<T> stateTile_;
  std::vector<const void*> matrixBases_;
  std::vector<BlockDescriptor> blockDescriptors_;
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
  // -fopenmp-simd for the loops of a kernel in the shape of a tile; it needs no OpenMP runtime
  std::string flags = "-O3 -ffp-contract=off -fno-math-errno -fopenmp-simd";
  flags += options.arch.empty() || options.arch == "native" ? " -march=native"
                                                            : " -march=" + options.arch;
  return flags;
}

} // namespace

namespace {

/// A file name for what is generated for `program`.
std::string fingerprintName(const Program& program, const char* suffix) {
  std::ostringstream name;
  name << "expr-" << std::hex << program.fingerprint() << std::dec << "-"
       << (program.computeType() == ComputeType::F32 ? "f32" : "f64") << "." << suffix;
  return name.str();
}

} // namespace

namespace {

std::string emitCpuSourceFor(const Program& program,
                             const LoweredProgram& lowered,
                             const ContractLayout& layout,
                             const std::vector<reader::datafield::ElementType>& gridElements,
                             KernelShape shape) {
  codegen::ContractAddressing contract;
  for (const auto& matrix : program.matrices()) {
    contract.shapes.push_back(matrix.shape);
  }
  for (const auto type : layout.matrices) {
    contract.matrixTypes.emplace_back(codegen::elementTypeName(type));
  }
  for (const auto type : layout.blocks) {
    contract.blockTypes.emplace_back(codegen::elementTypeName(type));
  }
  contract.point = "pointIndex != 0 ? pointIndex[l] : first + l";
  contract.indexType = "unsigned long";
  contract.matrix = cpuMatrix;
  contract.blockBase = cpuBlockBase;
  contract.cellStride = cpuCellStride;
  contract.modeStride = cpuModeStride;
  contract.cellIndex = cpuCellIndex;

  std::ostringstream out;
  out << "// Generated by seissol::expr for program fingerprint 0x" << std::hex
      << program.fingerprint() << std::dec << ".\n"
      << "// The arithmetic below is stringified from SEISSOL_EXPR_PW_LIST, the same\n"
      << "// table the interpreter evaluates, so the two cannot disagree.\n"
      << "#include <cmath>\n\n"
      << "struct SeissolExprBlock {\n"
      << "  const char* base;\n"
      << "  unsigned long cellStride;\n"
      << "  unsigned long modeStride;\n"
      << "  const unsigned int* cellIndex;\n"
      << "};\n\n";
  // the grids, read through their descriptors; the arithmetic of their samples is compiled without
  // contraction of floating-point expressions anyway (-ffp-contract=off)
  std::vector<reader::datafield::ElementType> elements = gridElements;
  elements.resize(program.grids().size(), reader::datafield::ElementType::Double);
  codegen::emitIntegerType(out, false);
  codegen::emitGridDescriptor(
      out, elements, lowered.gridUses(), "unsigned long long", "static inline", "");
  if (lowered.hasStrict()) {
    codegen::emitStrictMacros(out, "((x) + (y))", "((x) - (y))", "((x) * (y))");
  }
  emitStage(out,
            "seissol_expr_precompute",
            lowered.precompute(),
            lowered.operands(),
            lowered.stateSlotCount(),
            lowered.dataType(),
            contract,
            shape);
  emitStage(out,
            "seissol_expr_run",
            lowered.run(),
            lowered.operands(),
            lowered.stateSlotCount(),
            lowered.dataType(),
            contract,
            shape);
  return out.str();
}

} // namespace

bool tiledStage(const StageCode& stage,
                const std::vector<std::int32_t>& operands,
                KernelShape shape) {
  if (shape == KernelShape::Point) {
    return false;
  }
  const bool contracts = std::any_of(stage.code.begin(), stage.code.end(), [](const auto& inst) {
    return inst.op == Opcode::Contract;
  });
  if (contracts || stage.code.empty()) {
    return false;
  }
  if (shape == KernelShape::Tile) {
    return true;
  }
  const bool calls = std::any_of(stage.code.begin(), stage.code.end(), [](const auto& inst) {
    return inst.op == Opcode::Pw && codegen::callsLibrary(inst.fn);
  });
  return calls && codegen::criticalPath(stage, operands) >= TileDepth;
}

std::string
    emitCpuSource(const Program& program, const LoweredProgram& lowered, KernelShape shape) {
  // Without a binding, the matrices and blocks are taken to be stored in the compute type.
  ContractLayout layout;
  const auto type = program.computeType() == ComputeType::F32 ? reader::scripting::DataType::F32
                                                              : reader::scripting::DataType::F64;
  layout.matrices.assign(program.matrices().size(), type);
  layout.blocks.assign(program.blocks().size(), type);
  return emitCpuSourceFor(program, lowered, layout, {}, shape);
}

std::string emitCpuSource(const Program& program,
                          const LoweredProgram& lowered,
                          const Binding& binding,
                          KernelShape shape,
                          const std::vector<reader::datafield::ElementType>& gridElements) {
  return emitCpuSourceFor(program, lowered, contractLayoutOf(binding), gridElements, shape);
}

std::size_t rtcCpuCacheSize() {
  const std::lock_guard<std::mutex> lock(cacheMutex());
  return cache().size();
}

std::unique_ptr<Kernel> makeRtcCpuKernel(const Program& program,
                                         Binding& binding,
                                         LoweredProgram lowered,
                                         const BackendOptions& options,
                                         reader::datafield::GridStore* grids,
                                         const std::vector<std::size_t>& gridIds) {
  // the grids it samples, and the types of their elements, which the source bakes in
  std::vector<const reader::datafield::Grid*> sampled(program.grids().size(), nullptr);
  std::vector<std::string> paths(program.grids().size());
  std::vector<reader::datafield::ElementType> elements(program.grids().size(),
                                                       reader::datafield::ElementType::Double);
  for (std::size_t g = 0; g < lowered.gridUses().size(); ++g) {
    if (lowered.gridUses()[g].dimensions == 0) {
      continue;
    }
    paths[g] = program.grids()[g].path;
    if (grids == nullptr || g >= gridIds.size()) {
      logWarning() << "expr: the compiled CPU kernel of a program that samples the grid" << paths[g]
                   << "has no store to read it from; using the interpreter.";
      return nullptr;
    }
    sampled[g] = &grids->get(gridIds[g]);
    const reader::datafield::ArrayView* view = sampled[g]->arrayView();
    if (view == nullptr) {
      logWarning() << "expr: the grid" << paths[g]
                   << "keeps no array for a compiled kernel to read; using the interpreter.";
      return nullptr;
    }
    elements[g] = view->type;
  }

  const std::string flags = defaultFlags(options);
  const ContractLayout layout = contractLayoutOf(binding);
  const CacheKey key{program.fingerprint(),
                     options.lowering.fingerprint(),
                     layout.fingerprint(),
                     gridFingerprint(elements),
                     program.computeType(),
                     flags,
                     options.shape};

  Artifact artifact;
  // how the kernel came about: compiled, from the cache across runs, or made before in this run
  const char* origin = "made before";
  {
    // Held across the compile on purpose. This runs from prepare(), i.e. at
    // init, where serialising a handful of 40 ms compiles costs nothing and
    // two threads compiling the same program would be pure waste. If this path
    // is ever reached during the timestep loop, that is the bug to fix rather
    // than the lock.
    const std::lock_guard<std::mutex> lock(cacheMutex());
    auto found = cache().find(key);
    if (found == cache().end()) {
      const std::string source =
          emitCpuSourceFor(program, lowered, layout, elements, options.shape);
      codegen::dump(fingerprintName(program, "cpu.cpp"), source);
      bool cached = false;
      found = cache().emplace(key, compileAndLoad(source, flags, cached)).first;
      origin = cached ? "from the cache" : "compiled";
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

  if (!options.quiet) {
    logInfo()
        << "expr: compiled CPU kernel --" << lowered.summary().c_str() << "--"
        << cost(program, lowered, program.computeType()).summary(program.computeType()).c_str()
        << "-- with" << flags << "(" << origin << ")";
  }

  // in the type the data moves in
  if (lowered.dataType() == ComputeType::F32) {
    return std::make_unique<RtcCpuKernel<float>>(
        binding, std::move(lowered), artifact, tileSize, std::move(sampled), std::move(paths));
  }
  return std::make_unique<RtcCpuKernel<double>>(
      binding, std::move(lowered), artifact, tileSize, std::move(sampled), std::move(paths));
}

} // namespace seissol::expr
