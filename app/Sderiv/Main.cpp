// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// sderiv, the command line tool of the programs SeisSol evaluates (docs/sderiv-spec.rst): checks a
// program, prints it in the normal form of the core format, evaluates it at points, compares two
// programs, and reads the evaluation records of a run: evaluates them again with the interpreter
// of this build, merges those of the ranks into the sample of all points, and compares two runs.

#include "Expr/Backend.h"
#include "Expr/Core.h"
#include "Expr/Program.h"
#include "Expr/Record.h"
#include "Expr/SderivFrontend.h"
#include "Parallel/MPI.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using seissol::expr::Program;
using seissol::expr::Record;

constexpr const char* Usage = R"(usage: sderiv <command> [arguments]

Programs (sderiv in the full language or the core format, or a Lua model that traces):
  check FILE                  load and validate a program, and describe it
  print FILE                  the program in the normal form of the core format
  eval FILE [NAME=VALUE...] [--points POINTS] [--hex]
                              evaluate the program with the interpreter: at one point, its inputs
                              given by name, or at the points of a file (a first line with the names
                              of the columns, then a point per line, separated by spaces or commas)
  compare FILE FILE [--points POINTS | --record RECORD]
                              whether two programs are the same, and where they are not, their
                              values at the points of a file or of the sample of a record

Evaluation records (written into <prefix>-records by every rank of a run):
  replay RECORD...            evaluate the sample of each record again with the interpreter of this
                              build, bit for bit against the values the record holds
  merge [-o DIRECTORY] RECORD...
                              the record of the sample of all points, from the records of the ranks,
                              per program and backend
  diff RECORD RECORD          two records of one program, of two runs: their grids (the files, by
                              their hashes, and the geometry), and the points they share, layer by
                              layer (inputs, the interpreter's values, the backend's)

The exit status is 0 where everything agrees, 1 where something differs or fails.
)";

/// MPI, which a grid needs to load (AsagiLite reads collectively), started on first use; and the
/// grids of the programs, closed before MPI.
bool mpiStarted = false;
std::unique_ptr<seissol::reader::datafield::GridStore> grids;

int finish(int status) {
  grids.reset();
  if (mpiStarted) {
    seissol::Mpi::finalize();
    mpiStarted = false;
  }
  return status;
}

[[noreturn]] void fail(const std::string& what) {
  std::cerr << "sderiv: " << what << "\n";
  std::exit(finish(1));
}

bool endsWith(const std::string& text, const std::string& suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string readFile(const std::string& path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    fail("cannot read " + path);
  }
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

/// The grids of programs, loaded where a program samples one.
seissol::reader::datafield::GridStore& gridsFor(const Program& program) {
  if (!program.grids().empty() && !mpiStarted) {
    int argc = 0;
    char** argv = nullptr;
    seissol::Mpi::mpi.init(argc, argv);
    mpiStarted = true;
  }
  if (grids == nullptr) {
    grids = std::make_unique<seissol::reader::datafield::GridStore>();
  }
  return *grids;
}

/// The program of the file `path`.
Program loadProgram(const std::string& path) {
  if (endsWith(path, ".lua") || path.rfind("lua:", 0) == 0) {
    std::string reason;
    auto program = seissol::reader::scripting::buildProgram(path, &reason);
    if (!program.has_value()) {
      fail(path + ": the Lua model does not trace to a program: " + reason);
    }
    return *program;
  }
  if (endsWith(path, ".yaml") || endsWith(path, ".yml") || path.rfind("easi:", 0) == 0) {
    fail(path + ": an easi model does not compile to a program");
  }
  const std::string file = path.rfind("sderiv:", 0) == 0 ? path.substr(7) : path;
  const std::string text = readFile(file);
  seissol::expr::SderivOptions options;
  options.file = file;
  try {
    return seissol::expr::compileSderivModule(text, options);
  } catch (const seissol::expr::CoreError& error) {
    fail(error.what());
  } catch (const seissol::expr::SderivError& error) {
    if (error.position() >= 0) {
      const auto where =
          seissol::expr::positionOf(text, static_cast<std::size_t>(error.position()), file);
      fail(file + ":" + std::to_string(where.line) + ":" + std::to_string(where.column) + ": " +
           error.what());
    }
    fail(file + ": " + error.what());
  } catch (const std::exception& error) {
    fail(file + ": " + error.what());
  }
}

std::string fingerprintOf(const Program& program) {
  std::array<char, 32> buffer{};
  std::snprintf(buffer.data(),
                buffer.size(),
                "%016llx",
                static_cast<unsigned long long>(program.fingerprint()));
  return buffer.data();
}

std::string valueText(double value, bool hex) {
  std::array<char, 64> buffer{};
  std::snprintf(buffer.data(), buffer.size(), hex ? "%a" : "%.17g", value);
  return buffer.data();
}

/// A number as sderiv spells it, signed: decimal or hexadecimal, `inf` or `nan`.
double numberOf(const std::string& spelling) {
  const bool negative = !spelling.empty() && spelling[0] == '-';
  const std::string body =
      negative || (!spelling.empty() && spelling[0] == '+') ? spelling.substr(1) : spelling;
  double value = 0.0;
  if (body == "nan") {
    value = std::numeric_limits<double>::quiet_NaN();
  } else if (body == "inf") {
    value = std::numeric_limits<double>::infinity();
  } else {
    value = seissol::expr::numberOf(body);
  }
  return negative ? -value : value;
}

bool sameBits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0 || (std::isnan(a) && std::isnan(b));
}

/// Points by the names of their columns: columns[name][point].
struct Points {
  std::vector<std::string> names;
  std::vector<std::vector<double>> columns;
  std::size_t count{0};

  [[nodiscard]] const std::vector<double>* column(const std::string& name) const {
    const auto found = std::find(names.begin(), names.end(), name);
    return found == names.end() ? nullptr : &columns[found - names.begin()];
  }
};

std::vector<std::string> fieldsOf(const std::string& line) {
  std::string spaced = line;
  std::replace(spaced.begin(), spaced.end(), ',', ' ');
  std::istringstream stream(spaced);
  std::vector<std::string> fields;
  std::string field;
  while (stream >> field) {
    fields.push_back(field);
  }
  return fields;
}

Points readPoints(const std::string& path) {
  std::istringstream lines(readFile(path));
  Points points;
  std::string line;
  std::size_t number = 0;
  while (std::getline(lines, line)) {
    ++number;
    const auto hash = line.find('#');
    const auto fields = fieldsOf(hash == std::string::npos ? line : line.substr(0, hash));
    if (fields.empty()) {
      continue;
    }
    if (points.names.empty()) {
      points.names = fields;
      points.columns.resize(fields.size());
      continue;
    }
    if (fields.size() != points.names.size()) {
      fail(path + ":" + std::to_string(number) + ": " + std::to_string(fields.size()) +
           " values for " + std::to_string(points.names.size()) + " columns");
    }
    for (std::size_t c = 0; c < fields.size(); ++c) {
      try {
        points.columns[c].push_back(numberOf(fields[c]));
      } catch (const std::exception&) {
        fail(path + ":" + std::to_string(number) + ": `" + fields[c] + "` is no number");
      }
    }
    ++points.count;
  }
  return points;
}

/// The points of the sample of a record.
Points pointsOf(const Record& record) {
  Points points;
  points.names = record.inputs;
  points.columns.resize(record.inputs.size());
  for (const auto& point : record.sample) {
    for (std::size_t i = 0; i < record.inputs.size(); ++i) {
      points.columns[i].push_back(point.inputs[i]);
    }
  }
  points.count = record.sample.size();
  return points;
}

/// What the interpreter gives for `program` at `points`, [output][point].
std::vector<double>
    evaluate(const Program& program, const Points& points, const std::string& what) {
  if (!program.state().empty()) {
    fail(what + " has state, which its consumer advances from call to call; it is evaluated in a "
                "run");
  }
  std::vector<double> inputs;
  inputs.reserve(program.inputs().size() * points.count);
  for (const auto& input : program.inputs()) {
    const auto* column = points.column(input.name);
    if (column == nullptr) {
      fail(what + " reads the input `" + input.name + "`, which the points do not give");
    }
    inputs.insert(inputs.end(), column->begin(), column->end());
  }
  return seissol::expr::interpretValues(program, gridsFor(program), inputs, points.count);
}

Record readRecord(const std::string& path) {
  try {
    return seissol::expr::parseRecord(readFile(path), path);
  } catch (const std::exception& error) {
    fail(error.what());
  }
}

std::string headOf(const Record& record, const std::string& key) {
  for (const auto& [name, value] : record.head) {
    if (name == key) {
      return value;
    }
  }
  return {};
}

// ================================================================================== commands =====

int check(const std::vector<std::string>& args) {
  if (args.size() != 1) {
    fail("check takes one file");
  }
  const Program program = loadProgram(args[0]);
  std::cout << "program " << fingerprintOf(program) << ", "
            << (program.computeType() == seissol::expr::ComputeType::F32 ? "f32" : "f64") << "\n";
  const auto list = [](const auto& items, const auto& nameOf) {
    std::string text;
    for (const auto& item : items) {
      text += (text.empty() ? "" : ", ") + nameOf(item);
    }
    return text.empty() ? std::string("none") : text;
  };
  std::cout << "  inputs:  " << list(program.inputs(), [](const auto& input) { return input.name; })
            << "\n";
  std::cout << "  outputs: "
            << list(program.outputs(), [](const auto& output) { return output.name; }) << "\n";
  if (!program.state().empty()) {
    std::cout << "  states:  " << list(program.state(), [](const auto& state) {
      return state.name + " = " + valueText(state.initial, false);
    }) << "\n";
  }
  for (const auto& grid : program.grids()) {
    std::cout << "  grid " << (grid.name.empty() ? std::string("(unnamed)") : grid.name) << ": "
              << grid.path << (grid.variable.empty() ? "" : ", variable " + grid.variable)
              << (grid.timeDependent() ? ", time-dependent" : "") << "\n";
  }
  if (!program.checks().empty()) {
    std::cout << "  checks:  " << program.checks().size() << "\n";
  }
  try {
    static_cast<void>(seissol::expr::printCore(program));
  } catch (const std::exception& error) {
    std::cout << "  not in the core format: " << error.what() << "\n";
  }
  return 0;
}

int print(const std::vector<std::string>& args) {
  if (args.size() != 1) {
    fail("print takes one file");
  }
  const Program program = loadProgram(args[0]);
  try {
    std::cout << seissol::expr::printCore(program);
  } catch (const std::exception& error) {
    fail(args[0] + ": " + error.what());
  }
  return 0;
}

int eval(const std::vector<std::string>& args) {
  if (args.empty()) {
    fail("eval takes a file, and the inputs of a point or --points");
  }
  const Program program = loadProgram(args[0]);
  bool hex = false;
  std::optional<std::string> pointsFile;
  Points points;
  for (std::size_t a = 1; a < args.size(); ++a) {
    if (args[a] == "--hex") {
      hex = true;
    } else if (args[a] == "--points") {
      if (a + 1 == args.size()) {
        fail("--points takes a file");
      }
      pointsFile = args[++a];
    } else {
      const auto equals = args[a].find('=');
      if (equals == std::string::npos) {
        fail("expected NAME=VALUE, not `" + args[a] + "`");
      }
      points.names.push_back(args[a].substr(0, equals));
      try {
        points.columns.push_back({numberOf(args[a].substr(equals + 1))});
      } catch (const std::exception&) {
        fail("`" + args[a].substr(equals + 1) + "` is no number");
      }
      points.count = 1;
    }
  }
  if (pointsFile.has_value()) {
    if (!points.names.empty()) {
      fail("the inputs of a point and --points exclude each other");
    }
    points = readPoints(*pointsFile);
  }
  if (points.count == 0 && !program.inputs().empty()) {
    fail("no point: give the inputs of one as NAME=VALUE, or a file of them with --points");
  }
  if (program.inputs().empty()) {
    points.count = std::max<std::size_t>(points.count, 1);
  }
  const auto values = evaluate(program, points, args[0]);
  const std::size_t n = points.count;
  if (!pointsFile.has_value()) {
    for (std::size_t o = 0; o < program.outputs().size(); ++o) {
      std::cout << program.outputs()[o].name << " = " << valueText(values[o * n], hex) << "\n";
    }
    return 0;
  }
  // a table: the inputs the program reads, then its outputs
  for (const auto& input : program.inputs()) {
    std::cout << input.name << " ";
  }
  for (std::size_t o = 0; o < program.outputs().size(); ++o) {
    std::cout << program.outputs()[o].name << (o + 1 < program.outputs().size() ? " " : "");
  }
  std::cout << "\n";
  for (std::size_t k = 0; k < n; ++k) {
    for (const auto& input : program.inputs()) {
      std::cout << valueText((*points.column(input.name))[k], hex) << " ";
    }
    for (std::size_t o = 0; o < program.outputs().size(); ++o) {
      std::cout << valueText(values[o * n + k], hex)
                << (o + 1 < program.outputs().size() ? " " : "");
    }
    std::cout << "\n";
  }
  return 0;
}

int compare(const std::vector<std::string>& args) {
  if (args.size() != 2 && args.size() != 4) {
    fail("compare takes two files, and --points or --record");
  }
  const Program first = loadProgram(args[0]);
  const Program second = loadProgram(args[1]);
  if (first.fingerprint() == second.fingerprint()) {
    std::cout << "the same program, " << fingerprintOf(first) << "\n";
    return 0;
  }
  std::cout << "two programs: " << fingerprintOf(first) << " and " << fingerprintOf(second) << "\n";
  std::set<std::string> firstOutputs;
  for (const auto& output : first.outputs()) {
    firstOutputs.insert(output.name);
  }
  std::set<std::string> secondOutputs;
  for (const auto& output : second.outputs()) {
    secondOutputs.insert(output.name);
  }
  for (const auto& name : firstOutputs) {
    if (secondOutputs.count(name) == 0) {
      std::cout << "  `" << name << "` only in the first\n";
    }
  }
  for (const auto& name : secondOutputs) {
    if (firstOutputs.count(name) == 0) {
      std::cout << "  `" << name << "` only in the second\n";
    }
  }
  if (args.size() == 2) {
    return 1;
  }
  Points points;
  if (args[2] == "--points") {
    points = readPoints(args[3]);
  } else if (args[2] == "--record") {
    points = pointsOf(readRecord(args[3]));
  } else {
    fail("expected --points or --record, not `" + args[2] + "`");
  }
  const auto a = evaluate(first, points, args[0]);
  const auto b = evaluate(second, points, args[1]);
  const std::size_t n = points.count;
  std::size_t differing = 0;
  for (std::size_t o = 0; o < first.outputs().size(); ++o) {
    const auto& name = first.outputs()[o].name;
    const auto other = std::find_if(second.outputs().begin(),
                                    second.outputs().end(),
                                    [&](const auto& output) { return output.name == name; });
    if (other == second.outputs().end()) {
      continue;
    }
    const auto p = static_cast<std::size_t>(other - second.outputs().begin());
    std::size_t here = 0;
    for (std::size_t k = 0; k < n; ++k) {
      if (!sameBits(a[o * n + k], b[p * n + k])) {
        if (here++ == 0) {
          std::cout << "  `" << name << "` differs, first at point " << k << ": "
                    << valueText(a[o * n + k], true) << " and " << valueText(b[p * n + k], true)
                    << "\n";
        }
      }
    }
    if (here > 0) {
      std::cout << "  `" << name << "` differs at " << here << " of " << n << " points\n";
    }
    differing += here;
  }
  if (differing == 0) {
    std::cout << "the same bits at all " << n << " points\n";
  }
  return 1;
}

int replay(const std::vector<std::string>& args) {
  if (args.empty()) {
    fail("replay takes records");
  }
  int status = 0;
  for (const auto& path : args) {
    const Record record = readRecord(path);
    if (record.program.empty()) {
      std::cout << path << ": no program in the core format (" << record.notCore << ")\n";
      continue;
    }
    seissol::expr::Replay found;
    try {
      const Program program = seissol::expr::programOf(record);
      found = seissol::expr::replayRecord(record, gridsFor(program));
    } catch (const std::exception& error) {
      std::cout << path << ": " << error.what() << "\n";
      status = 1;
      continue;
    }
    if (found.differing == 0) {
      std::cout << path << ": the same bits at " << found.points << " points\n";
      continue;
    }
    status = 1;
    std::cout << path << ": " << found.differing << " values differ at " << found.points
              << " points\n";
    for (const auto& difference : found.differences) {
      std::cout << "  " << difference << "\n";
    }
  }
  return status;
}

int merge(const std::vector<std::string>& args) {
  std::optional<std::string> directory;
  std::vector<std::string> paths;
  for (std::size_t a = 0; a < args.size(); ++a) {
    if (args[a] == "-o") {
      if (a + 1 == args.size()) {
        fail("-o takes a directory");
      }
      directory = args[++a];
    } else {
      paths.push_back(args[a]);
    }
  }
  if (paths.empty()) {
    fail("merge takes records");
  }
  // per program and backend
  std::map<std::pair<std::string, std::string>, std::vector<Record>> groups;
  for (const auto& path : paths) {
    Record record = readRecord(path);
    groups[{record.fingerprint, record.backend}].push_back(std::move(record));
  }
  if (!directory.has_value() && groups.size() > 1) {
    fail("the records are of " + std::to_string(groups.size()) +
         " programs and backends; -o names the directory of their merged records");
  }
  if (directory.has_value()) {
    std::filesystem::create_directories(*directory);
  }
  int status = 0;
  for (const auto& [key, records] : groups) {
    const std::string name = key.first + " on " + key.second;
    try {
      const Program program = seissol::expr::programOf(records.front());
      const Record merged = seissol::expr::mergeRecords(records, gridsFor(program));
      const std::string text = seissol::expr::printRecord(merged);
      std::cerr << name << ": " << records.size() << " records, a sample of "
                << merged.sample.size() << " points; selects " << headOf(merged, "selects");
      for (const auto* note : {"duplicates", "inconsistent", "groups and extremes"}) {
        const auto value = headOf(merged, note);
        if (!value.empty()) {
          std::cerr << "; " << note << ": " << value;
        }
      }
      std::cerr << "\n";
      if (!headOf(merged, "inconsistent").empty()) {
        status = 1;
      }
      if (directory.has_value()) {
        const auto path =
            std::filesystem::path(*directory) / (key.first + "-" + key.second + "-merged.txt");
        std::ofstream file(path);
        file << text;
        if (!file) {
          fail("cannot write " + path.string());
        }
      } else {
        std::cout << text;
      }
    } catch (const std::exception& error) {
      std::cerr << name << ": " << error.what() << "\n";
      status = 1;
    }
  }
  return status;
}

int diff(const std::vector<std::string>& args) {
  if (args.size() != 2) {
    fail("diff takes two records");
  }
  const Record first = readRecord(args[0]);
  const Record second = readRecord(args[1]);
  int status = 0;
  if (first.fingerprint != second.fingerprint) {
    std::cout << "two programs: " << first.fingerprint << " and " << second.fingerprint << "\n";
    status = 1;
  }
  if (!first.scope.globalIdentities || !second.scope.globalIdentities) {
    std::cout << "the points of a record are identified by their index in the call, which another "
                 "run need not keep\n";
  }
  // the grids first: what the program samples
  for (const auto& grid : first.grids) {
    const auto other =
        std::find_if(second.grids.begin(), second.grids.end(), [&](const auto& candidate) {
          return candidate.grid == grid.grid;
        });
    if (other == second.grids.end()) {
      std::cout << "the grid " << grid.grid << " only in the first\n";
      status = 1;
      continue;
    }
    if (grid.sha256 != other->sha256) {
      std::cout << "the file of the grid " << grid.grid
                << " differs: " << (grid.sha256.empty() ? std::string("no hash") : grid.sha256)
                << " (" << grid.path << ") and "
                << (other->sha256.empty() ? std::string("no hash") : other->sha256) << " ("
                << other->path << ")\n";
      status = 1;
    }
    if (grid.variable != other->variable || grid.components != other->components) {
      std::cout << "the grid " << grid.grid << " reads other variables or components\n";
      status = 1;
    }
    for (std::size_t d = 0; d < std::max(grid.axes.size(), other->axes.size()); ++d) {
      if (d >= grid.axes.size() || d >= other->axes.size() || !(grid.axes[d] == other->axes[d])) {
        const auto describe = [](const std::vector<seissol::expr::GridRecord::Axis>& axes,
                                 std::size_t d) -> std::string {
          if (d >= axes.size()) {
            return "none";
          }
          const auto& axis = axes[d];
          return "`" + axis.file.dimension + "` from " + valueText(axis.min, false) + " by " +
                 valueText(axis.spacing, false) + " over " + std::to_string(axis.nodes) + " nodes" +
                 (axis.file.units.empty() ? "" : ", in " + axis.file.units) +
                 (axis.file.positive.empty() ? "" : ", positive " + axis.file.positive);
        };
        std::cout << "axis " << d << " of the grid " << grid.grid
                  << " differs: " << describe(grid.axes, d) << ", and " << describe(other->axes, d)
                  << "\n";
        status = 1;
      }
    }
  }
  for (const auto& grid : second.grids) {
    if (std::none_of(first.grids.begin(), first.grids.end(), [&](const auto& candidate) {
          return candidate.grid == grid.grid;
        })) {
      std::cout << "the grid " << grid.grid << " only in the second\n";
      status = 1;
    }
  }

  std::map<std::uint64_t, std::size_t> others;
  for (std::size_t k = 0; k < second.sample.size(); ++k) {
    others.emplace(second.sample[k].identity, k);
  }
  // the columns of the second by those of the first, by name
  const auto mapOf = [](const std::vector<std::string>& from, const std::vector<std::string>& to) {
    std::vector<std::optional<std::size_t>> map;
    for (const auto& name : from) {
      const auto found = std::find(to.begin(), to.end(), name);
      map.push_back(found == to.end()
                        ? std::nullopt
                        : std::optional<std::size_t>(static_cast<std::size_t>(found - to.begin())));
    }
    return map;
  };
  const auto inputMap = mapOf(first.inputs, second.inputs);
  const auto outputMap = mapOf(first.outputs, second.outputs);

  std::size_t shared = 0;
  std::array<std::size_t, 3> differing{};
  std::array<std::string, 3> firstDifference;
  const std::array<const char*, 3> layers = {
      "the inputs", "the interpreter's values", "the backend's values"};
  const auto note = [&](std::size_t layer, const std::string& what) {
    if (differing[layer]++ == 0) {
      firstDifference[layer] = what;
    }
  };
  for (const auto& point : first.sample) {
    const auto found = others.find(point.identity);
    if (found == others.end()) {
      continue;
    }
    ++shared;
    const auto& other = second.sample[found->second];
    const std::string at = " at the point of identity " + std::to_string(point.identity) + ": ";
    for (std::size_t i = 0; i < inputMap.size(); ++i) {
      if (inputMap[i].has_value() && !sameBits(point.inputs[i], other.inputs[*inputMap[i]])) {
        note(0,
             "`" + first.inputs[i] + "`" + at + valueText(point.inputs[i], true) + " and " +
                 valueText(other.inputs[*inputMap[i]], true));
      }
    }
    for (std::size_t o = 0; o < outputMap.size(); ++o) {
      if (!outputMap[o].has_value()) {
        continue;
      }
      const std::size_t p = *outputMap[o];
      if (!sameBits(point.interpreter[o], other.interpreter[p])) {
        note(1,
             "`" + first.outputs[o] + "`" + at + valueText(point.interpreter[o], true) + " and " +
                 valueText(other.interpreter[p], true));
      }
      if (first.backendValues && second.backendValues &&
          !sameBits(point.backend[o], other.backend[p])) {
        note(2,
             "`" + first.outputs[o] + "`" + at + valueText(point.backend[o], true) + " and " +
                 valueText(other.backend[p], true));
      }
    }
  }
  std::cout << shared << " points in both, " << first.sample.size() - shared
            << " only in the first, " << second.sample.size() - shared << " only in the second\n";
  for (std::size_t layer = 0; layer < layers.size(); ++layer) {
    if (differing[layer] > 0) {
      std::cout << layers[layer] << " differ at " << differing[layer] << " values; the first, "
                << firstDifference[layer] << "\n";
      status = 1;
    }
  }
  if (status == 0) {
    std::cout << "the same bits at the points in both\n";
  }
  return status;
}

} // namespace

int main(int argc, char* argv[]) {
  const std::vector<std::string> all(argv + 1, argv + argc);
  if (all.empty() || all[0] == "-h" || all[0] == "--help" || all[0] == "help") {
    std::cout << Usage;
    return all.empty() ? 1 : 0;
  }
  const std::string command = all[0];
  const std::vector<std::string> args(all.begin() + 1, all.end());
  int status = 1;
  if (command == "check") {
    status = check(args);
  } else if (command == "print") {
    status = print(args);
  } else if (command == "eval") {
    status = eval(args);
  } else if (command == "compare") {
    status = compare(args);
  } else if (command == "replay") {
    status = replay(args);
  } else if (command == "merge") {
    status = merge(args);
  } else if (command == "diff") {
    status = diff(args);
  } else {
    std::cerr << "sderiv: unknown command `" << command << "`\n\n" << Usage;
    return 1;
  }
  return finish(status);
}
