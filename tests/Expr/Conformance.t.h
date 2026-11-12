// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The conformance tests of the core format (the specification, "Conformance"): golden files, round
// trips, the values of the exact and the IEEE class at fixed points on every backend, the bounded
// class against a reference of high precision (MPFR), the check of a GPU backend against a stand-in
// and mutants, and the limits.

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Core.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Setup.h"

#ifdef ACL_DEVICE
#include "Memory/MemoryAllocator.h"

#include <Device/device.h>
#endif

#ifdef SEISSOL_TEST_MPFR
#include <mpfr.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace seissol::expr::test::conformance {

namespace df = reader::datafield;
namespace sc = reader::scripting;

/// The conformance programs and their golden files.
inline std::filesystem::path directory() {
  return std::filesystem::path(unit_test::TestSetup::SourcePath) / "Expr" / "Conformance";
}

inline std::string textOf(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::stringstream text;
  text << file.rdbuf();
  return text.str();
}

inline std::string hexOf(std::uint64_t bits) {
  std::array<char, 17> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%016llx", static_cast<unsigned long long>(bits));
  return buffer.data();
}

inline std::string hexOf(double value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return hexOf(bits);
}

inline double valueOf(const std::string& hex) {
  const std::uint64_t bits = std::stoull(hex, nullptr, 16);
  double value = 0.0;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

/// The same bits, or both NaN: the sign and the payload of a NaN are not kept (the core format).
inline bool same(double a, double b) {
  return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof(double)) == 0;
}

inline std::uint64_t splitMix(std::uint64_t& state) {
  std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31U);
}

/// Uniform in [0, 1).
inline double uniform(std::uint64_t& state) {
  return static_cast<double>(splitMix(state) >> 11U) * 0x1.0p-53;
}

/// `count` points of `inputs` inputs, [point][input]. Every input runs through special values
/// (signed zeros, infinities, NaN, the extremes of binary64, halves where rounding ties) at an
/// offset of its own, so that they meet each other; at every seventh point all inputs are equal; in
/// between, values from -4 to 4, every fourth a multiple of one half, so that comparisons and
/// selects take both branches.
inline std::vector<std::vector<double>> pointsFor(std::size_t inputs, std::size_t count) {
  constexpr double Inf = std::numeric_limits<double>::infinity();
  constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> specials = {0.0,
                                        -0.0,
                                        1.0,
                                        -1.0,
                                        0.5,
                                        -0.5,
                                        1.5,
                                        -2.5,
                                        2.0,
                                        3.0,
                                        0.1,
                                        -1e-3,
                                        Inf,
                                        -Inf,
                                        Nan,
                                        std::numeric_limits<double>::min(),
                                        std::numeric_limits<double>::denorm_min(),
                                        std::numeric_limits<double>::max(),
                                        -std::numeric_limits<double>::max(),
                                        1e-300,
                                        1e300,
                                        1.0 / 3.0,
                                        -2.0 / 3.0,
                                        7.0,
                                        -7.0};
  std::uint64_t state = 0x5eedULL;
  std::vector<std::vector<double>> points(count, std::vector<double>(inputs, 0.0));
  for (std::size_t k = 0; k < count; ++k) {
    for (std::size_t i = 0; i < inputs; ++i) {
      double value = 0.0;
      if ((k + i) % 3 == 0) {
        value = specials[(k * (i + 1) + 5 * i) % specials.size()];
      } else {
        value = 8.0 * uniform(state) - 4.0;
        if (k % 4 == 1) {
          value = std::round(2.0 * value) / 2.0;
        }
      }
      points[k][i] = k % 7 == 6 && i > 0 ? points[k][0] : value;
    }
  }
  return points;
}

/// The values of a program at fixed points, in one compute type.
struct Evaluation {
  std::string type;
  std::vector<std::string> inputs;
  std::vector<std::string> outputs;
  /// [point][input] and [point][output]
  std::vector<std::vector<double>> in;
  std::vector<std::vector<double>> out;
};

/// What a golden file says of a conformance program: the program as its fingerprint, its canonical
/// form and its normal form, and its values with the interpreter.
struct Golden {
  std::string fingerprint;
  std::string canonical;
  std::string normal;
  std::vector<Evaluation> evaluations;
};

/// The lines of `text`, each after `| ` (an empty one as `|`, which no editor strips further), with
/// a control character and the backslash as `\xhh` (the canonical form separates with 0x1f).
inline std::string blockOf(const std::string& text) {
  std::string out;
  std::istringstream lines(text);
  std::string line;
  while (std::getline(lines, line)) {
    std::string escaped;
    for (const char c : line) {
      if (static_cast<unsigned char>(c) < 0x20 || c == '\\' || c == 0x7f) {
        std::array<char, 5> code{};
        std::snprintf(code.data(), code.size(), "\\x%02x", static_cast<unsigned char>(c));
        escaped += code.data();
      } else {
        escaped += c;
      }
    }
    out += escaped.empty() ? "|\n" : "| " + escaped + "\n";
  }
  return out;
}

/// A line of a block, as blockOf() escaped it.
inline std::string unescaped(const std::string& line) {
  std::string out;
  for (std::size_t i = 0; i < line.size(); ++i) {
    if (line[i] == '\\' && i + 3 < line.size() && line[i + 1] == 'x') {
      out += static_cast<char>(std::stoi(line.substr(i + 2, 2), nullptr, 16));
      i += 3;
    } else {
      out += line[i];
    }
  }
  return out;
}

inline std::string goldenText(const std::string& name, const Golden& golden) {
  // the header of the file a golden file is (not that of this one)
  // REUSE-IgnoreStart
  std::string out = "# SPDX-FileCopyrightText: 2026 SeisSol Group\n"
                    "#\n"
                    "# SPDX-License-Identifier: BSD-3-Clause\n"
                    "# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/\n"
                    "#\n"
                    "# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff\n"
                    "#\n"
                    "# What the core format makes of " +
                    name +
                    ".sderiv (the specification, \"Conformance\"): the program as its\n"
                    "# fingerprint, its canonical form and its normal form, and its values at "
                    "fixed points with\n"
                    "# the interpreter, in hexadecimal, in either compute type. "
                    "SEISSOL_EXPR_GOLDEN=write writes it\n"
                    "# anew from what the interpreter of the build gives.\n";
  // REUSE-IgnoreEnd
  out += "fingerprint = " + golden.fingerprint + "\n";
  out += "canonical form\n" + blockOf(golden.canonical);
  out += "normal form\n" + blockOf(golden.normal);
  for (const auto& evaluation : golden.evaluations) {
    out += "evaluation " + evaluation.type + ":";
    for (const auto& input : evaluation.inputs) {
      out += " " + input;
    }
    out += " ->";
    for (const auto& output : evaluation.outputs) {
      out += " " + output;
    }
    out += "\n";
    for (std::size_t k = 0; k < evaluation.in.size(); ++k) {
      std::string line = "|";
      for (const double value : evaluation.in[k]) {
        line += " " + hexOf(value);
      }
      line += " ->";
      for (const double value : evaluation.out[k]) {
        line += " " + hexOf(value);
      }
      out += line + "\n";
    }
  }
  return out;
}

inline Golden goldenOf(const std::string& text) {
  Golden golden;
  std::istringstream lines(text);
  std::string line;
  std::string* block = nullptr;
  Evaluation* evaluation = nullptr;
  while (std::getline(lines, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line[0] == '|') {
      const std::string content = line.size() > 2 ? line.substr(2) : std::string();
      if (block != nullptr) {
        *block += unescaped(content) + "\n";
      } else if (evaluation != nullptr) {
        std::istringstream words(content);
        std::string word;
        bool outputs = false;
        std::vector<double> in;
        std::vector<double> out;
        while (words >> word) {
          if (word == "->") {
            outputs = true;
          } else {
            (outputs ? out : in).push_back(valueOf(word));
          }
        }
        evaluation->in.push_back(in);
        evaluation->out.push_back(out);
      }
      continue;
    }
    block = nullptr;
    evaluation = nullptr;
    if (line.rfind("fingerprint = ", 0) == 0) {
      golden.fingerprint = line.substr(14);
    } else if (line == "canonical form") {
      block = &golden.canonical;
    } else if (line == "normal form") {
      block = &golden.normal;
    } else if (line.rfind("evaluation ", 0) == 0) {
      const auto colon = line.find(':');
      golden.evaluations.emplace_back();
      evaluation = &golden.evaluations.back();
      evaluation->type = line.substr(11, colon - 11);
      std::istringstream words(line.substr(colon + 1));
      std::string word;
      bool outputs = false;
      while (words >> word) {
        if (word == "->") {
          outputs = true;
        } else {
          (outputs ? evaluation->outputs : evaluation->inputs).push_back(word);
        }
      }
    }
  }
  return golden;
}

#ifdef ACL_DEVICE
inline void setUpDevice() {
  static const bool done = [] {
    const auto& device = ::device::DeviceInstance::instance();
    device.api().setDevice(0);
    device.api().initialize();
    return true;
  }();
  static_cast<void>(done);
}
#endif

/// The outputs of `program` at `points` ([point][input]), [point][output], on `backend`: nothing
/// where this build has no such backend at hand (no compiler at run time, no device). The columns
/// lie where the backend reads them: on the device for a device backend.
inline std::optional<std::vector<std::vector<double>>> outputsAt(
    const Program& program, const std::vector<std::vector<double>>& points, BackendKind backend) {
  const std::size_t count = points.size();
  const std::size_t inputs = program.inputs().size();
  const std::size_t outputs = program.outputs().size();
  // the columns, of doubles, and the states in their types, from their initial values: kept in
  // the table, since a device backend fills none of its own
  std::vector<std::vector<double>> columns(inputs + outputs, std::vector<double>(count, -1.0));
  for (std::size_t k = 0; k < count; ++k) {
    for (std::size_t i = 0; i < inputs; ++i) {
      columns[i][k] = points[k][i];
    }
  }
  std::vector<std::vector<double>> wideStates;
  std::vector<std::vector<float>> narrowStates;
  std::vector<bool> narrow;
  for (const auto& state : program.state()) {
    narrow.push_back(state.type == ValueType::F32 ||
                     (state.type == ValueType::Real && program.computeType() == ComputeType::F32));
    wideStates.emplace_back(narrow.back() ? 0 : count, state.initial);
    narrowStates.emplace_back(narrow.back() ? count : 0, static_cast<float>(state.initial));
  }
  // where the backend reads and writes them: on the device for a device backend
  const bool device = backend != BackendKind::Interpreter && backend != BackendKind::RtcCpu;
  std::vector<void*> bases;
  for (auto& column : columns) {
    bases.push_back(column.data());
  }
  for (std::size_t s = 0; s < narrow.size(); ++s) {
    bases.push_back(narrow[s] ? static_cast<void*>(narrowStates[s].data())
                              : static_cast<void*>(wideStates[s].data()));
  }
  const auto bytesOf = [&](std::size_t base) {
    const bool state = base >= columns.size();
    return count * (state && narrow[base - columns.size()] ? sizeof(float) : sizeof(double));
  };
#ifdef ACL_DEVICE
  if (device) {
    setUpDevice();
    for (std::size_t b = 0; b < bases.size(); ++b) {
      void* onDevice = memory::allocate(
          std::max<std::size_t>(bytesOf(b), 8), 64, memory::Memkind::DeviceGlobalMemory);
      memory::memcopy(onDevice,
                      bases[b],
                      bytesOf(b),
                      memory::Memkind::DeviceGlobalMemory,
                      memory::Memkind::Standard);
      bases[b] = onDevice;
    }
  }
#else
  if (device) {
    return std::nullopt;
  }
#endif
  sc::DataTable table(count);
  for (std::size_t i = 0; i < inputs; ++i) {
    table.bindViewConst<double>(
        program.inputs()[i].name, sc::Direction::In, static_cast<const double*>(bases[i]));
  }
  for (std::size_t o = 0; o < outputs; ++o) {
    table.bindView<double>(
        program.outputs()[o].name, sc::Direction::Out, static_cast<double*>(bases[inputs + o]));
  }
  for (std::size_t s = 0; s < narrow.size(); ++s) {
    void* base = bases[columns.size() + s];
    if (narrow[s]) {
      table.bindState<float>(program.state()[s].name, static_cast<float*>(base), 1, 1, 1);
    } else {
      table.bindState<double>(program.state()[s].name, static_cast<double*>(base), 1, 1, 1);
    }
  }
  Binding binding = Binding::bind(program, table);
  df::GridStore store;
  BackendOptions options;
  options.quiet = true;
  options.preferred = backend;
  const auto kernel = makeKernel(program, binding, store, options);
  const bool ran = kernel != nullptr && kernel->kind() == backend;
  if (ran) {
    kernel->precompute(table);
    kernel->run(table);
  }
#ifdef ACL_DEVICE
  if (device) {
    for (std::size_t b = 0; b < bases.size(); ++b) {
      if (ran && b >= inputs && b < columns.size()) {
        memory::memcopy(columns[b].data(),
                        bases[b],
                        bytesOf(b),
                        memory::Memkind::Standard,
                        memory::Memkind::DeviceGlobalMemory);
      }
      memory::free(bases[b], memory::Memkind::DeviceGlobalMemory);
    }
  }
#endif
  if (!ran) {
    return std::nullopt;
  }
  std::vector<std::vector<double>> values(count, std::vector<double>(outputs, 0.0));
  for (std::size_t k = 0; k < count; ++k) {
    for (std::size_t o = 0; o < outputs; ++o) {
      values[k][o] = columns[inputs + o][k];
    }
  }
  return values;
}

/// The backends of this build: the interpreter, the compiled CPU kernels, and that of the device.
inline std::vector<BackendKind> backends() {
  std::vector<BackendKind> kinds = {BackendKind::Interpreter, BackendKind::RtcCpu};
  if (deviceBackendOfBuild() != BackendKind::RtcCpu) {
    kinds.push_back(deviceBackendOfBuild());
  }
  return kinds;
}

/// Whether the value of `root` holds a product into a sum, which a GPU may contract into a fused
/// multiply-add (the specification, "Accuracy"): its bits may differ there. Through a conversion
/// as well, which may be none in the compute type (`real(x * y) + z` in binary64).
inline bool contractible(const Program& program, NodeId root) {
  const Arena& arena = program.arena();
  std::vector<NodeId> stack = {root};
  std::set<NodeId> seen;
  const auto isProduct = [&](NodeId id) {
    while (id != NoNode && arena[id].kind == Kind::Cast) {
      id = arena[id].a;
    }
    return id != NoNode && arena[id].kind == Kind::PW && arena[id].fn == Fn::Mul;
  };
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (id == NoNode || !seen.insert(id).second) {
      continue;
    }
    const Node& node = arena[id];
    if (node.kind == Kind::PW && (node.fn == Fn::Add || node.fn == Fn::Sub) &&
        (isProduct(node.a) || isProduct(node.b))) {
      return true;
    }
    for (const NodeId child : {node.a, node.b, node.c}) {
      stack.push_back(child);
    }
  }
  return false;
}

inline ComputeType computeTypeOf(const std::string& type) {
  return type == "f32" ? ComputeType::F32 : ComputeType::F64;
}

TEST_SUITE("expr conformance") {
  TEST_CASE("golden files: the program of a core file, and its values on every backend") {
    const char* mode = std::getenv("SEISSOL_EXPR_GOLDEN");
    const bool write = mode != nullptr && std::string(mode) == "write";
    std::vector<std::filesystem::path> sources;
    for (const auto& entry : std::filesystem::directory_iterator(directory())) {
      if (entry.path().extension() == ".sderiv") {
        sources.push_back(entry.path());
      }
    }
    std::sort(sources.begin(), sources.end());
    REQUIRE(!sources.empty());
    for (const auto& source : sources) {
      const std::string stem = source.stem().string();
      CAPTURE(stem);
      const Program program = loadCore(textOf(source), source.string());
      Golden now;
      now.fingerprint = hexOf(program.fingerprint());
      now.canonical = program.canonicalForm();
      now.normal = printCore(program);

      // round trips: the normal form gives the same program, and prints as itself
      const Program again = loadCore(now.normal, source.string());
      CHECK(hexOf(again.fingerprint()) == now.fingerprint);
      CHECK(again.canonicalForm() == now.canonical);
      CHECK(printCore(again) == now.normal);

      const auto goldenPath = std::filesystem::path(source).replace_extension(".golden");
      if (write) {
        for (const char* type : {"f64", "f32"}) {
          Program typed = program;
          typed.setComputeType(computeTypeOf(type));
          Evaluation evaluation;
          evaluation.type = type;
          for (const auto& input : program.inputs()) {
            evaluation.inputs.push_back(input.name);
          }
          for (const auto& output : program.outputs()) {
            evaluation.outputs.push_back(output.name);
          }
          evaluation.in = pointsFor(program.inputs().size(), 64);
          const auto values = outputsAt(typed, evaluation.in, BackendKind::Interpreter);
          REQUIRE(values.has_value());
          evaluation.out = *values;
          now.evaluations.push_back(evaluation);
        }
        std::ofstream(goldenPath) << goldenText(stem, now);
        MESSAGE("wrote " << goldenPath.string());
        continue;
      }

      INFO("SEISSOL_EXPR_GOLDEN=write writes the golden file anew, where the change is intended");
      REQUIRE(std::filesystem::exists(goldenPath));
      const Golden golden = goldenOf(textOf(goldenPath));
      CHECK(now.fingerprint == golden.fingerprint);
      CHECK(now.canonical == golden.canonical);
      CHECK(now.normal == golden.normal);
      REQUIRE(golden.evaluations.size() == 2);
      for (const auto& evaluation : golden.evaluations) {
        CAPTURE(evaluation.type);
        Program typed = program;
        typed.setComputeType(computeTypeOf(evaluation.type));
        REQUIRE(evaluation.outputs.size() == program.outputs().size());
        for (const auto backend : backends()) {
          CAPTURE(name(backend));
          const auto values = outputsAt(typed, evaluation.in, backend);
          if (!values.has_value()) {
            MESSAGE("this build has no backend " << name(backend) << " at hand");
            continue;
          }
          const bool onDevice =
              backend != BackendKind::Interpreter && backend != BackendKind::RtcCpu;
          std::size_t differing = 0;
          std::string first;
          for (std::size_t o = 0; o < evaluation.outputs.size(); ++o) {
            // a GPU contracts a product into a sum; its first call is checked against the
            // perturbed reference instead (makeKernel())
            if (onDevice && contractible(program, program.roots()[o])) {
              continue;
            }
            for (std::size_t k = 0; k < evaluation.in.size(); ++k) {
              if (!same((*values)[k][o], evaluation.out[k][o])) {
                if (differing++ == 0) {
                  first = "`" + evaluation.outputs[o] + "` at point " + std::to_string(k) + ": " +
                          hexOf((*values)[k][o]) + " for " + hexOf(evaluation.out[k][o]);
                }
              }
            }
          }
          INFO(first);
          CHECK(differing == 0);
        }
      }
    }
  }

  TEST_CASE("the check of a GPU backend accepts a stand-in, and rejects mutants") {
    // u = x * y + z, which a GPU contracts; v = sin(x) * y, whose sine another math library moves;
    // w = select(lt(x * y, 1), v, u), which decides otherwise where x * y is 1 by a rounding
    const Program program = loadCore("#!sderiv 1\nin x, y, z\ndef p = x * y\ndef u = p + z\n"
                                     "def s = sin(x)\ndef v = s * y\ndef c = lt(p, 1)\n"
                                     "def w = select(c, v, u)\nout u = u\nout v = v\nout w = w\n",
                                     "");
    // the points: one where the product rounds to 1, the others of every magnitude
    constexpr std::size_t Points = 64;
    std::uint64_t state = 7;
    std::vector<double> x(Points);
    std::vector<double> y(Points);
    std::vector<double> z(Points);
    for (std::size_t k = 0; k < Points; ++k) {
      x[k] = (uniform(state) - 0.5) * 20.0;
      y[k] = std::ldexp(uniform(state) + 0.5, static_cast<int>(k % 40) - 20);
      z[k] = (uniform(state) - 0.5) * 4.0;
    }
    constexpr std::size_t Critical = 5;
    x[Critical] = 1.0 + std::ldexp(1.0, -27);
    y[Critical] = 1.0 - std::ldexp(1.0, -27);
    z[Critical] = 0.25;
    std::vector<double> inputs;
    for (const auto* column : {&x, &y, &z}) {
      inputs.insert(inputs.end(), column->begin(), column->end());
    }
    CallPool pool;
    for (std::size_t k = 0; k < Points; ++k) {
      pool.points.push_back(k);
      pool.identities.push_back(k);
    }
    // moved by `ulps` units in the last place
    const auto moved = [](double value, int ulps) {
      for (int i = 0; i < std::abs(ulps); ++i) {
        value = std::nextafter(value,
                               ulps > 0 ? std::numeric_limits<double>::infinity()
                                        : -std::numeric_limits<double>::infinity());
      }
      return value;
    };
    // what a backend writes, with the sine `sine` and the select decided by `decide`
    const auto written = [&](const std::function<double(std::size_t)>& sine,
                             const std::function<bool(std::size_t, double)>& decide) {
      std::vector<double> values(3 * Points);
      for (std::size_t k = 0; k < Points; ++k) {
        const double u = std::fma(x[k], y[k], z[k]);
        const double v = sine(k) * y[k];
        values[k] = u;
        values[Points + k] = v;
        values[2 * Points + k] = decide(k, x[k] * y[k]) ? v : u;
      }
      return values;
    };
    const auto asWritten = [](std::size_t /*k*/, double product) { return product < 1.0; };
    df::GridStore store;
    const auto judge = [&](const std::vector<double>& values) {
      const auto found = compareValues(program, store, pool, inputs, values, true);
      REQUIRE(found.has_value());
      return *found;
    };

    // the stand-in: contracted, a sine up to 3 ulp off (the bound k_fn = 4), and the select
    // decided otherwise where a rounding decides it -- accepted, the point unstable
    const auto standIn = judge(written(
        [&](std::size_t k) { return moved(std::sin(x[k]), static_cast<int>(k % 7) - 3); },
        [](std::size_t k, double product) { return k == Critical ? true : product < 1.0; }));
    CHECK(standIn.deviating == 0);
    CHECK(std::find(standIn.unstable.begin(), standIn.unstable.end(), Critical) !=
          standIn.unstable.end());
    CHECK(standIn.worstRatio <= 1.0);

    // mutants: a sine in binary32
    const auto narrow = judge(written(
        [&](std::size_t k) { return static_cast<double>(std::sin(static_cast<float>(x[k]))); },
        asWritten));
    CHECK(narrow.deviating > 0);
    // a sine moved by three times the bound (twice the bound lies on the tolerance, c = 2)
    const auto far =
        judge(written([&](std::size_t k) { return moved(std::sin(x[k]), 12); }, asWritten));
    CHECK(far.deviating > 0);
    // a deviation confined to the values far below the largest: |v| spans 2^-20 to 2^20, and only
    // those below 2^-10 are off, by a relative 1e-9, far below 1e-9 of the largest
    const auto small = judge(written(
        [&](std::size_t k) {
          const double sine = std::sin(x[k]);
          return std::fabs(sine * y[k]) < std::ldexp(1.0, -10) ? sine * (1.0 + 1e-9) : sine;
        },
        asWritten));
    CHECK(small.deviating > 0);
  }

  TEST_CASE("a file with a million definitions loads, prints and runs") {
    // in the normal form: one operation per definition, halving and adding x by turns
    constexpr int Definitions = 1000000;
    std::string text = "#!sderiv 1\nin x\ndef t0 = x + 1\n";
    for (int i = 1; i < Definitions; ++i) {
      text += "def t" + std::to_string(i) + " = t" + std::to_string(i - 1) +
              (i % 2 == 1 ? " * 0.5\n" : " + x\n");
    }
    text += "out u = t" + std::to_string(Definitions - 1) + "\n";
    const Program program = loadCore(text, "");
    CHECK(printCore(program) == text);
    const auto values = outputsAt(program, {{2.0}, {-1.0}}, BackendKind::Interpreter);
    REQUIRE(values.has_value());
    // by turns half of the last and that plus x: x, after a halving, where they settle
    CHECK((*values)[0][0] == 2.0);
    CHECK((*values)[1][0] == -1.0);
  }

#ifdef SEISSOL_TEST_MPFR
  /// A bounded operation and where it is measured.
  struct Bounded {
    const char* name;
    int (*one)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
    int (*two)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
    std::function<std::pair<double, double>(std::uint64_t&)> draw;
  };

  /// The error of `computed` in units in the last place of binary64 at `exact`.
  inline double ulpsOff(double computed, mpfr_srcptr exact) {
    if (mpfr_nan_p(exact) != 0) {
      return std::isnan(computed) ? 0.0 : std::numeric_limits<double>::infinity();
    }
    if (mpfr_inf_p(exact) != 0) {
      return std::isinf(computed) && (computed > 0) == (mpfr_sgn(exact) > 0)
                 ? 0.0
                 : std::numeric_limits<double>::infinity();
    }
    if (!std::isfinite(computed)) {
      return std::numeric_limits<double>::infinity();
    }
    // exact = m 2^e with 1/2 <= |m| < 1: the unit in the last place is 2^(e - 53), and at least
    // that of the subnormal numbers
    const long exponent =
        mpfr_zero_p(exact) != 0 ? -1074 : std::max(mpfr_get_exp(exact) - 53, -1074L);
    mpfr_t difference;
    mpfr_init2(difference, 2200);
    mpfr_set_d(difference, computed, MPFR_RNDN);
    mpfr_sub(difference, difference, exact, MPFR_RNDN);
    mpfr_mul_2si(difference, difference, -exponent, MPFR_RNDN);
    const double ulps = std::fabs(mpfr_get_d(difference, MPFR_RNDN));
    mpfr_clear(difference);
    return ulps;
  }

  TEST_CASE("the bounded operations lie within 4 ulp of the exact value, on every backend") {
    // the arguments: of every magnitude where the function is defined and finite
    const auto between = [](double low, double high) {
      return [low, high](std::uint64_t& state) {
        return std::pair<double, double>{low + (high - low) * uniform(state), 0.0};
      };
    };
    const auto magnitudes = [](int low, int high, bool negative) {
      return [low, high, negative](std::uint64_t& state) {
        const double value = std::exp2(low + (high - low) * uniform(state));
        return std::pair<double, double>{negative && uniform(state) < 0.5 ? -value : value, 0.0};
      };
    };
    const std::vector<Bounded> functions = {
        {"exp", mpfr_exp, nullptr, between(-700.0, 700.0)},
        {"log", mpfr_log, nullptr, magnitudes(-1000, 1000, false)},
        {"log2", mpfr_log2, nullptr, magnitudes(-1000, 1000, false)},
        {"log10", mpfr_log10, nullptr, magnitudes(-1000, 1000, false)},
        {"sin", mpfr_sin, nullptr, between(-100.0, 100.0)},
        {"cos", mpfr_cos, nullptr, between(-100.0, 100.0)},
        {"tan", mpfr_tan, nullptr, between(-100.0, 100.0)},
        {"asin", mpfr_asin, nullptr, between(-1.0, 1.0)},
        {"acos", mpfr_acos, nullptr, between(-1.0, 1.0)},
        {"atan", mpfr_atan, nullptr, magnitudes(-40, 40, true)},
        {"sinh", mpfr_sinh, nullptr, between(-700.0, 700.0)},
        {"cosh", mpfr_cosh, nullptr, between(-700.0, 700.0)},
        {"tanh", mpfr_tanh, nullptr, between(-20.0, 20.0)},
        {"asinh", mpfr_asinh, nullptr, magnitudes(-500, 500, true)},
        {"acosh",
         mpfr_acosh,
         nullptr,
         [](std::uint64_t& state) {
           return std::pair<double, double>{1.0 + std::exp2(-50.0 + 550.0 * uniform(state)), 0.0};
         }},
        {"atanh", mpfr_atanh, nullptr, between(-0.999, 0.999)},
        {"erf", mpfr_erf, nullptr, between(-6.0, 6.0)},
        {"pow",
         nullptr,
         mpfr_pow,
         [](std::uint64_t& state) {
           return std::pair<double, double>{std::exp2(-20.0 + 40.0 * uniform(state)),
                                            -30.0 + 60.0 * uniform(state)};
         }},
        {"atan2",
         nullptr,
         mpfr_atan2,
         [](std::uint64_t& state) {
           const auto value = [&state] {
             const double magnitude = std::exp2(-50.0 + 100.0 * uniform(state));
             return uniform(state) < 0.5 ? -magnitude : magnitude;
           };
           const double a = value();
           return std::pair<double, double>{a, value()};
         }},
    };
    // one program of all of them, each with inputs of its own
    constexpr std::size_t Points = 4096;
    std::string text = "#!sderiv 1\nin ";
    for (std::size_t f = 0; f < functions.size(); ++f) {
      text += std::string(f == 0 ? "" : ", ") + "a" + std::to_string(f) + ", b" + std::to_string(f);
    }
    text += "\n";
    for (std::size_t f = 0; f < functions.size(); ++f) {
      const std::string a = "a" + std::to_string(f);
      const std::string b = "b" + std::to_string(f);
      text += "def t" + std::to_string(f) + " = " + functions[f].name + "(" + a +
              (functions[f].two != nullptr ? ", " + b : std::string()) + ")\n";
    }
    for (std::size_t f = 0; f < functions.size(); ++f) {
      text += "out u_" + std::string(functions[f].name) + " = t" + std::to_string(f) + "\n";
    }
    const Program program = loadCore(text, "");
    std::uint64_t state = 0xacc0ULL;
    std::vector<std::vector<double>> points(Points, std::vector<double>(2 * functions.size()));
    for (std::size_t k = 0; k < Points; ++k) {
      for (std::size_t f = 0; f < functions.size(); ++f) {
        const auto [a, b] = functions[f].draw(state);
        points[k][2 * f] = a;
        points[k][2 * f + 1] = b;
      }
    }
    // the exact values, of 256 bits
    // (mpfr_t is an array; the vectors hold what it is an array of)
    std::vector<std::vector<__mpfr_struct>> exact(functions.size(),
                                                  std::vector<__mpfr_struct>(Points));
    mpfr_t a;
    mpfr_t b;
    mpfr_init2(a, 53);
    mpfr_init2(b, 53);
    for (std::size_t f = 0; f < functions.size(); ++f) {
      for (std::size_t k = 0; k < Points; ++k) {
        mpfr_ptr value = &exact[f][k];
        mpfr_init2(value, 256);
        mpfr_set_d(a, points[k][2 * f], MPFR_RNDN);
        mpfr_set_d(b, points[k][2 * f + 1], MPFR_RNDN);
        if (functions[f].one != nullptr) {
          functions[f].one(value, a, MPFR_RNDN);
        } else {
          functions[f].two(value, a, b, MPFR_RNDN);
        }
      }
    }
    mpfr_clear(a);
    mpfr_clear(b);

    for (const auto backend : backends()) {
      CAPTURE(name(backend));
      const auto values = outputsAt(program, points, backend);
      if (!values.has_value()) {
        MESSAGE("this build has no backend " << name(backend) << " at hand");
        continue;
      }
      std::string report;
      for (std::size_t f = 0; f < functions.size(); ++f) {
        CAPTURE(functions[f].name);
        double worst = 0.0;
        std::size_t at = 0;
        for (std::size_t k = 0; k < Points; ++k) {
          const double ulps = ulpsOff((*values)[k][f], &exact[f][k]);
          if (!(ulps <= worst)) {
            worst = ulps;
            at = k;
          }
        }
        std::array<char, 64> line{};
        std::snprintf(line.data(), line.size(), " %s %.3f", functions[f].name, worst);
        report += line.data();
        INFO("at " << hexOf(points[at][2 * f]) << ", " << hexOf(points[at][2 * f + 1]) << ": "
                   << hexOf((*values)[at][f]));
        CHECK(worst <= 4.0);
      }
      MESSAGE(name(backend) << ", the largest error in ulp:" << report);
    }
    for (auto& values : exact) {
      for (auto& value : values) {
        mpfr_clear(&value);
      }
    }
  }

  TEST_CASE("the IEEE operations are correctly rounded, on every backend and in both types") {
    // in binary64 against MPFR with the exponent range and the subnormals of binary64; in binary32
    // against the correctly rounded binary64 result rounded once more, which for these operations
    // is the correctly rounded binary32 one (53 >= 2 * 24 + 2)
    const std::vector<std::string> operations = {
        "a + b", "a - b", "a * b", "a / b", "sqrt(a)", "rcp(b)", "mod(a, b)"};
    std::string text = "#!sderiv 1\nin a, b\n";
    for (std::size_t o = 0; o < operations.size(); ++o) {
      text += "def t" + std::to_string(o) + " = " + operations[o] + "\n";
    }
    for (std::size_t o = 0; o < operations.size(); ++o) {
      text += "out u" + std::to_string(o) + " = t" + std::to_string(o) + "\n";
    }
    const Program program = loadCore(text, "");
    constexpr std::size_t Points = 8192;
    std::uint64_t state = 0x1eeeULL;
    std::vector<std::vector<double>> wide(Points, std::vector<double>(2));
    std::vector<std::vector<double>> narrow(Points, std::vector<double>(2));
    for (std::size_t k = 0; k < Points; ++k) {
      for (std::size_t i = 0; i < 2; ++i) {
        const double sign = uniform(state) < 0.5 ? -1.0 : 1.0;
        // every magnitude, the subnormals too, and pairs of close magnitudes for cancellation
        const double mantissa = 1.0 + uniform(state);
        wide[k][i] = sign * std::ldexp(mantissa, static_cast<int>(uniform(state) * 2098.0) - 1074);
        narrow[k][i] = sign * std::ldexp(static_cast<double>(static_cast<float>(mantissa)),
                                         static_cast<int>(uniform(state) * 300.0) - 150);
        if (k % 3 == 0 && i == 1) {
          wide[k][1] = -wide[k][0] * (1.0 + std::ldexp(uniform(state), -30));
          narrow[k][1] =
              -narrow[k][0] *
              static_cast<double>(1.0F + std::ldexp(static_cast<float>(uniform(state)), -12));
        }
        narrow[k][i] = static_cast<double>(static_cast<float>(narrow[k][i]));
      }
    }
    // the correctly rounded value of operation `o` in binary64
    const auto rounded = [](std::size_t o, double x, double y) {
      const auto oldMin = mpfr_get_emin();
      const auto oldMax = mpfr_get_emax();
      mpfr_set_emin(-1073);
      mpfr_set_emax(1024);
      mpfr_t a;
      mpfr_t b;
      mpfr_t r;
      mpfr_inits2(53, a, b, r, static_cast<mpfr_ptr>(nullptr));
      mpfr_set_d(a, x, MPFR_RNDN);
      mpfr_set_d(b, y, MPFR_RNDN);
      int inexact = 0;
      switch (o) {
      case 0:
        inexact = mpfr_add(r, a, b, MPFR_RNDN);
        break;
      case 1:
        inexact = mpfr_sub(r, a, b, MPFR_RNDN);
        break;
      case 2:
        inexact = mpfr_mul(r, a, b, MPFR_RNDN);
        break;
      case 3:
        inexact = mpfr_div(r, a, b, MPFR_RNDN);
        break;
      case 4:
        inexact = mpfr_sqrt(r, a, MPFR_RNDN);
        break;
      case 5:
        inexact = mpfr_ui_div(r, 1, b, MPFR_RNDN);
        break;
      default: {
        // the sign of b, as in Lua: the remainder, exact, and b added to it where the signs differ
        inexact = mpfr_fmod(r, a, b, MPFR_RNDN);
        if (mpfr_zero_p(r) == 0 && mpfr_nan_p(r) == 0 && (mpfr_sgn(r) < 0) != (y < 0)) {
          inexact = mpfr_add(r, r, b, MPFR_RNDN);
        }
      }
      }
      inexact = mpfr_subnormalize(r, inexact, MPFR_RNDN);
      static_cast<void>(inexact);
      const double value = mpfr_get_d(r, MPFR_RNDN);
      mpfr_clears(a, b, r, static_cast<mpfr_ptr>(nullptr));
      mpfr_set_emin(oldMin);
      mpfr_set_emax(oldMax);
      return value;
    };
    for (const char* type : {"f64", "f32"}) {
      CAPTURE(type);
      Program typed = program;
      typed.setComputeType(computeTypeOf(type));
      const bool f32 = std::string(type) == "f32";
      const auto& points = f32 ? narrow : wide;
      for (const auto backend : backends()) {
        CAPTURE(name(backend));
        const auto values = outputsAt(typed, points, backend);
        if (!values.has_value()) {
          MESSAGE("this build has no backend " << name(backend) << " at hand");
          continue;
        }
        for (std::size_t o = 0; o < operations.size(); ++o) {
          CAPTURE(operations[o]);
          // OpenCL C allows division and square roots of binary32 that are not correctly rounded
          if (f32 && backend == BackendKind::RtcOpenCl && (o == 3 || o == 4 || o == 5)) {
            continue;
          }
          std::size_t differing = 0;
          std::string first;
          for (std::size_t k = 0; k < Points; ++k) {
            const double x = points[k][0];
            const double y = points[k][1];
            const double expected =
                f32 ? static_cast<double>(static_cast<float>(rounded(o, x, y))) : rounded(o, x, y);
            if (!same((*values)[k][o], expected) && differing++ == 0) {
              first = hexOf(x) + ", " + hexOf(y) + ": " + hexOf((*values)[k][o]) + " for " +
                      hexOf(expected);
            }
          }
          INFO(first);
          CHECK(differing == 0);
        }
      }
    }
  }
#endif
}

} // namespace seissol::expr::test::conformance
