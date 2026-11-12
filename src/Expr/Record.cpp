// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Record.h"

#include "Expr/Backend.h"
#include "Expr/Core.h"
#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

constexpr const char* Title = "# The evaluation record of a program (docs/sderiv-spec.rst, "
                              "\"Evaluation record\"), written by SeisSol at the first call of its "
                              "kernel.";
constexpr const char* ProgramMarker = "# the program, in the normal form of the core format";
constexpr const char* NotCorePrefix = "# not in the core format: ";
constexpr const char* SamplePrefix = "# the sample";
constexpr const char* UnstablePrefix = "# unstable points";
constexpr const char* ClampComment =
    "# per grid and axis, the points of the closed pool beyond its volume: below, above, NaN";
/// The points of the smallest hash in a sample (as Backend.cpp chooses them).
constexpr std::size_t SmallestHashes = 64;

std::string hexOf(double x) {
  std::array<char, 64> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%a", x);
  return buffer.data();
}

std::string hashText(std::uint64_t hash) {
  std::array<char, 32> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "0x%016llx", static_cast<unsigned long long>(hash));
  return buffer.data();
}

std::string fingerprintText(const Program& program) {
  std::array<char, 32> buffer{};
  std::snprintf(buffer.data(),
                buffer.size(),
                "%016llx",
                static_cast<unsigned long long>(program.fingerprint()));
  return buffer.data();
}

/// The value of the spelling `%a` gives: hexadecimal floating point, `inf` and `nan`, signed.
double valueOf(const std::string& spelling) {
  const bool negative = !spelling.empty() && spelling[0] == '-';
  const std::string body =
      negative || (!spelling.empty() && spelling[0] == '+') ? spelling.substr(1) : spelling;
  double value = 0.0;
  if (body == "nan") {
    value = std::numeric_limits<double>::quiet_NaN();
  } else if (body == "inf") {
    value = std::numeric_limits<double>::infinity();
  } else {
    value = numberOf(body);
  }
  return negative ? -value : value;
}

/// Whether two values have the same bits, any NaN equal to any other.
bool sameBits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0 || (std::isnan(a) && std::isnan(b));
}

std::vector<std::string> wordsOf(const std::string& line) {
  std::istringstream stream(line);
  std::vector<std::string> words;
  std::string word;
  while (stream >> word) {
    words.push_back(word);
  }
  return words;
}

bool startsWith(const std::string& text, const std::string& prefix) {
  return text.compare(0, prefix.size(), prefix) == 0;
}

class Reader {
  public:
  Reader(const std::string& text, std::string file) : file_(std::move(file)) {
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
      lines_.push_back(line);
    }
  }

  [[nodiscard]] bool done() const { return next_ >= lines_.size(); }
  [[nodiscard]] const std::string& peek() const { return lines_[next_]; }
  const std::string& take() { return lines_[next_++]; }
  [[nodiscard]] std::size_t lineNumber() const { return next_; }

  [[noreturn]] void fail(const std::string& what) const {
    throw std::runtime_error((file_.empty() ? std::string("record") : file_) + ":" +
                             std::to_string(next_) + ": " + what);
  }

  private:
  std::string file_;
  std::vector<std::string> lines_;
  std::size_t next_{0};
};

std::size_t countOf(Reader& reader, const std::string& value, const std::string& unit) {
  const auto words = wordsOf(value);
  if (words.size() != 2 || words[1] != unit) {
    reader.fail("expected `<number> " + unit + "`, not `" + value + "`");
  }
  try {
    return static_cast<std::size_t>(std::stoull(words[0]));
  } catch (const std::exception&) {
    reader.fail("expected a number, not `" + words[0] + "`");
  }
}

void readPool(Reader& reader, const std::string& value, PoolScope& scope) {
  if (value == "all points") {
    scope.closedAll = true;
    return;
  }
  // N points of the smallest hashes, up to 0x...
  const auto words = wordsOf(value);
  if (words.size() != 9 || words[1] != "points" || words[7] != "to") {
    reader.fail("expected `all points` or `<number> points of the smallest hashes, up to <hash>`");
  }
  scope.closedAll = false;
  try {
    scope.closed = static_cast<std::size_t>(std::stoull(words[0]));
    scope.bound = std::stoull(words[8], nullptr, 16);
  } catch (const std::exception&) {
    reader.fail("expected a number and a hash in `" + value + "`");
  }
}

/// `0 12 0, 0 0 0, 31 0 0 of 4096 points`, the clamps of the grid `grid`.
ClampCount readClamps(Reader& reader, const std::string& grid, const std::string& value) {
  ClampCount count;
  count.grid = grid;
  const auto of = value.rfind(" of ");
  if (of == std::string::npos) {
    reader.fail("expected `<below> <above> <nan>, ... of <number> points`, not `" + value + "`");
  }
  count.of = countOf(reader, value.substr(of + 4), "points");
  std::istringstream axes(value.substr(0, of));
  std::string axis;
  while (std::getline(axes, axis, ',')) {
    const auto words = wordsOf(axis);
    if (words.size() != 3) {
      reader.fail("expected three counts per axis, not `" + axis + "`");
    }
    std::array<std::size_t, 3> sides{};
    try {
      for (std::size_t i = 0; i < 3; ++i) {
        sides[i] = static_cast<std::size_t>(std::stoull(words[i]));
      }
    } catch (const std::exception&) {
      reader.fail("expected three counts per axis, not `" + axis + "`");
    }
    count.axes.push_back(sides);
  }
  return count;
}

std::string poolText(const PoolScope& scope) {
  if (scope.closedAll) {
    return "all points";
  }
  return std::to_string(scope.closed) + " points of the smallest hashes, up to " +
         hashText(scope.bound);
}

/// The index of the input `group`, if the record has one.
std::optional<std::size_t> groupColumn(const Record& record) {
  const auto found = std::find(record.inputs.begin(), record.inputs.end(), "group");
  if (found == record.inputs.end()) {
    return std::nullopt;
  }
  return static_cast<std::size_t>(found - record.inputs.begin());
}

} // namespace

Record recordOf(const Program& program,
                const std::string& backend,
                const BackendComparison& found,
                const std::vector<std::string>& build,
                int rank) {
  Record record;
  record.fingerprint = fingerprintText(program);
  record.backend = backend;
  record.rank = rank;
  record.scope = found.scope;
  record.compared = found.compared;
  record.compute = program.computeType() == ComputeType::F32 ? "f32" : "f64";
  for (const auto& line : build) {
    const auto equals = line.find(" = ");
    if (equals == std::string::npos) {
      record.head.emplace_back(line, "");
    } else {
      record.head.emplace_back(line.substr(0, equals), line.substr(equals + 3));
    }
  }
  const bool compiled = backend != name(BackendKind::Interpreter);
  if (compiled) {
    std::string check =
        found.deviating == 0
            ? (found.perturbed && found.worstRatio > 0.0 ? "within the tolerance" : "the same bits")
            : "deviates at " + std::to_string(found.deviating) + " values";
    if (found.perturbed && (found.worstRatio > 0.0 || !found.unstable.empty())) {
      check += ", at most " + std::to_string(found.worstRatio) + " of the tolerance, " +
               std::to_string(found.unstable.size()) + " unstable points";
    }
    if (!found.worst.empty()) {
      check += "; the largest: " + found.worst;
    }
    record.head.emplace_back("check", check);
  }
  try {
    record.program = printCore(program);
    // what the tool checks by the fingerprint: the normal form reads back into the program
    Program again = loadCore(record.program, "");
    again.setComputeType(program.computeType());
    if (again.fingerprint() != program.fingerprint()) {
      record.head.emplace_back("normal form",
                               "reads back as program " + fingerprintText(again) +
                                   ", which the fingerprint of the record does not check");
    }
  } catch (const std::invalid_argument& error) {
    record.notCore = error.what();
  }
  for (const auto& input : program.inputs()) {
    record.inputs.push_back(input.name);
  }
  for (const auto& output : program.outputs()) {
    record.outputs.push_back(output.name);
  }
  record.backendValues = compiled;
  record.tolerances = compiled && found.perturbed;
  record.clamps = found.clamps;

  const std::size_t n = found.points.size();
  const std::size_t inputs = record.inputs.size();
  const std::size_t outputs = record.outputs.size();
  for (std::size_t k = 0; k < n; ++k) {
    Record::Point point;
    point.rank = rank;
    point.point = found.points[k];
    point.identity = k < found.identities.size() ? found.identities[k] : found.points[k];
    for (std::size_t i = 0; i < inputs && i * n + k < found.inputs.size(); ++i) {
      point.inputs.push_back(found.inputs[i * n + k]);
    }
    for (std::size_t o = 0; o < outputs; ++o) {
      point.interpreter.push_back(found.interpreter[o * n + k]);
      if (record.backendValues) {
        point.backend.push_back(found.backend[o * n + k]);
      }
      if (record.tolerances) {
        point.tolerance.push_back(found.tolerance[o * n + k]);
      }
    }
    record.sample.push_back(std::move(point));
  }
  for (std::size_t u = 0; u < found.unstable.size(); ++u) {
    record.unstable.emplace_back(found.unstable[u],
                                 u < found.unstableIdentities.size() ? found.unstableIdentities[u]
                                                                     : found.unstable[u]);
  }
  return record;
}

std::string printRecord(const Record& record) {
  std::ostringstream out;
  out << Title << "\n";
  out << "record = " << RecordVersion << "\n";
  out << "fingerprint = " << record.fingerprint << "\n";
  out << "backend = " << record.backend << "\n";
  if (record.rank.has_value()) {
    out << "rank = " << *record.rank << "\n";
  } else {
    out << "merged =";
    for (const int rank : record.ranks) {
      out << " " << rank;
    }
    out << "\n";
  }
  out << "call = " << record.scope.callPoints << " points\n";
  out << "identities = " << (record.scope.globalIdentities ? "global" : "local") << "\n";
  out << "pool = " << poolText(record.scope) << "\n";
  out << "scan = " << (record.scope.scanned ? "all points" : "the pool") << "\n";
  out << "compared = " << record.compared << " points\n";
  out << "compute = " << record.compute << "\n";
  for (const auto& [key, value] : record.head) {
    out << key << (value.empty() ? "" : " = " + value) << "\n";
  }
  if (!record.clamps.empty()) {
    out << ClampComment << "\n";
  }
  for (const auto& count : record.clamps) {
    out << "clamped " << count.grid << " =";
    for (std::size_t d = 0; d < count.axes.size(); ++d) {
      out << (d == 0 ? " " : ", ") << count.axes[d][0] << " " << count.axes[d][1] << " "
          << count.axes[d][2];
    }
    out << " of " << count.of << " points\n";
  }

  out << "\n" << ProgramMarker << "\n";
  if (record.program.empty()) {
    out << NotCorePrefix << record.notCore << "\n";
  } else {
    out << record.program;
  }

  // the sample: per point its inputs and, per output, the value of the interpreter and of the
  // backend (and its tolerance, against the perturbed reference)
  out << "\n"
      << SamplePrefix
      << ": per point its inputs, and per output the interpreter's value, the backend's"
      << (record.tolerances ? " and its tolerance" : "") << ", in hexadecimal\n";
  out << (record.rank.has_value() ? "" : "rank ") << "point identity";
  for (const auto& input : record.inputs) {
    out << " " << input;
  }
  for (const auto& output : record.outputs) {
    out << " " << output << ":interpreter";
    if (record.backendValues) {
      out << " " << output << ":backend";
      if (record.tolerances) {
        out << " " << output << ":tolerance";
      }
    }
  }
  out << "\n";
  for (const auto& point : record.sample) {
    if (!record.rank.has_value()) {
      out << point.rank << " ";
    }
    out << point.point << " " << point.identity;
    for (const double value : point.inputs) {
      out << " " << hexOf(value);
    }
    for (std::size_t o = 0; o < record.outputs.size(); ++o) {
      out << " " << hexOf(point.interpreter[o]);
      if (record.backendValues) {
        out << " " << hexOf(point.backend[o]);
        if (record.tolerances) {
          out << " " << hexOf(point.tolerance[o]);
        }
      }
    }
    out << "\n";
  }
  if (!record.unstable.empty()) {
    out << "\n"
        << UnstablePrefix
        << ", where a perturbed evaluation decides a select otherwise: the point and its "
           "identity\n";
    for (const auto& [point, identity] : record.unstable) {
      out << point << " " << identity << "\n";
    }
  }
  return out.str();
}

Record parseRecord(const std::string& text, const std::string& file) {
  Reader reader(text, file);
  Record record;

  // the head: `key = value` up to the first empty line
  bool seenFingerprint = false;
  while (!reader.done() && !reader.peek().empty()) {
    const std::string& line = reader.take();
    if (startsWith(line, "#")) {
      continue;
    }
    const auto equals = line.find(" = ");
    const std::string key = equals == std::string::npos ? line : line.substr(0, equals);
    const std::string value = equals == std::string::npos ? "" : line.substr(equals + 3);
    if (key == "record") {
      int version = 0;
      try {
        version = std::stoi(value);
      } catch (const std::exception&) {
        reader.fail("expected the version of the record, not `" + value + "`");
      }
      if (version > RecordVersion) {
        reader.fail("a record of version " + value + ", which this build reads only up to " +
                    std::to_string(RecordVersion));
      }
    } else if (key == "fingerprint") {
      record.fingerprint = value;
      seenFingerprint = true;
    } else if (startsWith(key, "clamped ")) {
      record.clamps.push_back(readClamps(reader, key.substr(8), value));
    } else if (key == "backend") {
      record.backend = value;
    } else if (key == "rank") {
      try {
        record.rank = std::stoi(value);
      } catch (const std::exception&) {
        reader.fail("expected the number of a rank, not `" + value + "`");
      }
    } else if (key == "merged" || line == "merged =") {
      record.rank.reset();
      for (const auto& word : wordsOf(value)) {
        record.ranks.push_back(std::stoi(word));
      }
    } else if (key == "call") {
      record.scope.callPoints = countOf(reader, value, "points");
    } else if (key == "identities") {
      if (value != "global" && value != "local") {
        reader.fail("expected `global` or `local`, not `" + value + "`");
      }
      record.scope.globalIdentities = value == "global";
    } else if (key == "pool") {
      readPool(reader, value, record.scope);
    } else if (key == "scan") {
      record.scope.scanned = value == "all points";
    } else if (key == "compared") {
      record.compared = countOf(reader, value, "points");
    } else if (key == "compute") {
      if (value != "f64" && value != "f32") {
        reader.fail("expected the compute type `f64` or `f32`, not `" + value + "`");
      }
      record.compute = value;
    } else {
      record.head.emplace_back(key, value);
    }
  }
  if (!seenFingerprint) {
    reader.fail("a record starts with its head, which names the fingerprint of its program");
  }

  // the program, up to the sample
  while (!reader.done() && reader.peek() != ProgramMarker) {
    reader.take();
  }
  if (reader.done()) {
    reader.fail("the record has no program");
  }
  reader.take();
  std::vector<std::string> program;
  while (!reader.done() && !startsWith(reader.peek(), SamplePrefix)) {
    program.push_back(reader.take());
  }
  while (!program.empty() && program.back().empty()) {
    program.pop_back();
  }
  if (program.size() == 1 && startsWith(program[0], NotCorePrefix)) {
    record.notCore = program[0].substr(std::string(NotCorePrefix).size());
  } else {
    for (const auto& line : program) {
      record.program += line + "\n";
    }
  }
  if (reader.done()) {
    reader.fail("the record has no sample");
  }
  reader.take();

  // the columns of the sample
  if (reader.done()) {
    reader.fail("the sample has no columns");
  }
  const auto columns = wordsOf(reader.take());
  const bool merged = !columns.empty() && columns[0] == "rank";
  const std::size_t lead = merged ? 3 : 2;
  if (columns.size() < lead || columns[lead - 2] != "point" || columns[lead - 1] != "identity") {
    reader.fail("the columns of a sample start with `point identity`");
  }
  if (merged == record.rank.has_value()) {
    reader.fail(merged ? "a record of one rank has no column `rank`"
                       : "a merged record has the column `rank`");
  }
  std::size_t c = lead;
  for (; c < columns.size() && columns[c].find(':') == std::string::npos; ++c) {
    record.inputs.push_back(columns[c]);
  }
  for (; c < columns.size(); ++c) {
    const auto colon = columns[c].rfind(':');
    const std::string name = columns[c].substr(0, colon);
    const std::string kind = columns[c].substr(colon + 1);
    if (kind == "interpreter") {
      record.outputs.push_back(name);
    } else if (kind == "backend") {
      record.backendValues = true;
    } else if (kind == "tolerance") {
      record.tolerances = true;
    } else {
      reader.fail("unknown column `" + columns[c] + "`");
    }
  }
  const std::size_t perOutput = 1 + (record.backendValues ? 1 : 0) + (record.tolerances ? 1 : 0);
  if (columns.size() != lead + record.inputs.size() + perOutput * record.outputs.size()) {
    reader.fail("the columns of the outputs are not complete");
  }

  // the points
  while (!reader.done() && !reader.peek().empty() && !startsWith(reader.peek(), "#")) {
    const auto words = wordsOf(reader.take());
    if (words.size() != columns.size()) {
      reader.fail("expected " + std::to_string(columns.size()) + " values, not " +
                  std::to_string(words.size()));
    }
    Record::Point point;
    try {
      std::size_t w = 0;
      if (merged) {
        point.rank = std::stoi(words[w++]);
      } else {
        point.rank = *record.rank;
      }
      point.point = static_cast<std::size_t>(std::stoull(words[w++]));
      point.identity = std::stoull(words[w++]);
      for (std::size_t i = 0; i < record.inputs.size(); ++i) {
        point.inputs.push_back(valueOf(words[w++]));
      }
      for (std::size_t o = 0; o < record.outputs.size(); ++o) {
        point.interpreter.push_back(valueOf(words[w++]));
        if (record.backendValues) {
          point.backend.push_back(valueOf(words[w++]));
        }
        if (record.tolerances) {
          point.tolerance.push_back(valueOf(words[w++]));
        }
      }
    } catch (const std::exception& error) {
      reader.fail(std::string("a value of the sample does not read: ") + error.what());
    }
    record.sample.push_back(std::move(point));
  }

  // the unstable points
  while (!reader.done() && !startsWith(reader.peek(), UnstablePrefix)) {
    reader.take();
  }
  if (!reader.done()) {
    reader.take();
    while (!reader.done() && !reader.peek().empty()) {
      const auto words = wordsOf(reader.take());
      if (words.empty() || words.size() > 2) {
        reader.fail("expected a point and its identity");
      }
      const auto point = static_cast<std::size_t>(std::stoull(words[0]));
      record.unstable.emplace_back(point, words.size() == 2 ? std::stoull(words[1]) : point);
    }
  }
  return record;
}

Program programOf(const Record& record) {
  if (record.program.empty()) {
    throw std::runtime_error("the record of program " + record.fingerprint +
                             " has no program in the core format" +
                             (record.notCore.empty() ? "" : ": " + record.notCore));
  }
  Program program = loadCore(record.program, "");
  program.setComputeType(record.compute == "f32" ? ComputeType::F32 : ComputeType::F64);
  const std::string fingerprint = fingerprintText(program);
  if (fingerprint != record.fingerprint) {
    throw std::runtime_error("the program of the record has the fingerprint " + fingerprint +
                             ", not " + record.fingerprint +
                             " as recorded: the core format does not read it back as it was");
  }
  return program;
}

Record mergeRecords(const std::vector<Record>& records, reader::datafield::GridStore& grids) {
  if (records.empty()) {
    throw std::runtime_error("no records to merge");
  }
  const Record& first = records.front();
  const auto rankName = [](const Record& record) {
    return record.rank.has_value() ? "rank " + std::to_string(*record.rank)
                                   : std::string("a merged record");
  };
  for (const auto& record : records) {
    if (record.fingerprint != first.fingerprint || record.backend != first.backend) {
      throw std::runtime_error(
          "the records are of different programs or backends: " + first.fingerprint + " on " +
          first.backend + ", " + record.fingerprint + " on " + record.backend);
    }
    if (record.inputs != first.inputs || record.outputs != first.outputs ||
        record.backendValues != first.backendValues || record.compute != first.compute) {
      throw std::runtime_error("the records of program " + first.fingerprint +
                               " do not have the same columns or compute type");
    }
    if (!record.scope.globalIdentities) {
      throw std::runtime_error(
          "the record of " + rankName(record) +
          " identifies its points by their index in the call, which no partition keeps: its "
          "consumer gives no identities, and its sample cannot be merged with those of other "
          "ranks");
    }
  }
  const Program program = programOf(first);

  Record merged;
  merged.fingerprint = first.fingerprint;
  merged.backend = first.backend;
  merged.compute = first.compute;
  merged.program = first.program;
  merged.inputs = first.inputs;
  merged.outputs = first.outputs;
  merged.backendValues = first.backendValues;
  merged.tolerances = std::all_of(
      records.begin(), records.end(), [](const Record& record) { return record.tolerances; });
  merged.scope.globalIdentities = true;
  merged.scope.closedAll = true;
  merged.scope.bound = std::numeric_limits<std::uint64_t>::max();
  merged.scope.scanned = true;
  std::set<int> ranks;
  for (const auto& record : records) {
    if (record.rank.has_value()) {
      ranks.insert(*record.rank);
    }
    ranks.insert(record.ranks.begin(), record.ranks.end());
    merged.scope.callPoints += record.scope.callPoints;
    merged.scope.closed += record.scope.closedAll ? record.scope.callPoints : record.scope.closed;
    merged.scope.scanned = merged.scope.scanned && record.scope.scanned;
    if (!record.scope.closedAll) {
      merged.scope.closedAll = false;
      merged.scope.bound = std::min(merged.scope.bound, record.scope.bound);
    }
    merged.compared += record.compared;
  }
  merged.ranks.assign(ranks.begin(), ranks.end());
  if (merged.scope.closedAll) {
    merged.scope.bound = 0;
  }

  // the clamps of every grid, summed
  for (const auto& record : records) {
    for (const auto& count : record.clamps) {
      auto found = std::find_if(merged.clamps.begin(),
                                merged.clamps.end(),
                                [&](const ClampCount& kept) { return kept.grid == count.grid; });
      if (found == merged.clamps.end()) {
        merged.clamps.push_back(ClampCount{count.grid, {}, 0});
        found = merged.clamps.end() - 1;
      }
      found->axes.resize(std::max(found->axes.size(), count.axes.size()), {0, 0, 0});
      for (std::size_t d = 0; d < count.axes.size(); ++d) {
        for (std::size_t i = 0; i < 3; ++i) {
          found->axes[d][i] += count.axes[d][i];
        }
      }
      found->of += count.of;
    }
  }

  // the head: the build as the first record has it, the checks of all
  std::vector<std::string> checks;
  for (const auto& [key, value] : first.head) {
    if (key != "check" && key != "selects" && key != "duplicates" && key != "inconsistent" &&
        key != "clamped points" && key != "groups and extremes") {
      merged.head.emplace_back(key, value);
    }
  }
  for (const auto& record : records) {
    for (const auto& [key, value] : record.head) {
      if (key == "check") {
        checks.push_back(value);
      }
    }
  }

  // the union of the samples, a point evaluated on two ranks once
  std::vector<Record::Point> points;
  std::vector<bool> closed;
  std::map<std::uint64_t, std::size_t> byIdentity;
  std::size_t duplicates = 0;
  std::size_t inconsistent = 0;
  std::string firstInconsistency;
  for (const auto& record : records) {
    for (const auto& point : record.sample) {
      const bool inClosed =
          record.scope.closedAll || splitMix64(point.identity) <= record.scope.bound;
      const auto [found, inserted] = byIdentity.emplace(point.identity, points.size());
      if (inserted) {
        points.push_back(point);
        closed.push_back(inClosed);
        continue;
      }
      ++duplicates;
      const auto& other = points[found->second];
      closed[found->second] = closed[found->second] || inClosed;
      bool same = true;
      for (std::size_t i = 0; i < point.inputs.size(); ++i) {
        same = same && sameBits(point.inputs[i], other.inputs[i]);
      }
      for (std::size_t o = 0; o < point.interpreter.size(); ++o) {
        same = same && sameBits(point.interpreter[o], other.interpreter[o]);
      }
      if (!same) {
        if (inconsistent++ == 0) {
          firstInconsistency = "the point of identity " + std::to_string(point.identity) +
                               " on ranks " + std::to_string(other.rank) + " and " +
                               std::to_string(point.rank);
        }
      }
    }
  }
  std::set<std::uint64_t> unstable;
  for (const auto& record : records) {
    for (const auto& entry : record.unstable) {
      unstable.insert(entry.second);
    }
  }

  // the rules of the sample over the union, in the order of the hashes
  const std::size_t n = points.size();
  std::vector<std::uint64_t> hashes(n);
  std::vector<std::size_t> byHash(n);
  for (std::size_t k = 0; k < n; ++k) {
    hashes[k] = splitMix64(points[k].identity);
    byHash[k] = k;
  }
  std::sort(byHash.begin(), byHash.end(), [&](std::size_t a, std::size_t b) {
    return hashes[a] < hashes[b];
  });
  Record view;
  view.backendValues = merged.backendValues;
  view.sample = points;
  std::vector<bool> chosen(n, false);
  for (std::size_t i = 0; i < std::min(SmallestHashes, n); ++i) {
    chosen[byHash[i]] = true;
  }
  if (const auto group = groupColumn(merged)) {
    std::map<double, int> taken;
    for (const std::size_t k : byHash) {
      const double key = points[k].inputs[*group];
      if (!std::isnan(key) && taken[key]++ < 2) {
        chosen[k] = true;
      }
    }
  }
  const std::size_t outputs = merged.outputs.size();
  for (std::size_t o = 0; o < outputs; ++o) {
    std::size_t smallest = n;
    std::size_t largest = n;
    std::size_t worst = n;
    double worstDeviation = 0.0;
    for (const std::size_t k : byHash) {
      const double value = view.written(k, o);
      if (!std::isnan(value)) {
        if (smallest == n || value < view.written(smallest, o)) {
          smallest = k;
        }
        if (largest == n || value > view.written(largest, o)) {
          largest = k;
        }
      }
      if (merged.backendValues && !sameBits(points[k].backend[o], points[k].interpreter[o])) {
        const double deviation = std::fabs(points[k].backend[o] - points[k].interpreter[o]);
        const double measure =
            std::isnan(deviation) ? std::numeric_limits<double>::infinity() : deviation;
        if (worst == n || measure > worstDeviation) {
          worst = k;
          worstDeviation = measure;
        }
      }
    }
    for (const std::size_t k : {smallest, largest, worst}) {
      if (k < n) {
        chosen[k] = true;
      }
    }
  }

  // the selects, and the grids beyond their volume, from the points of the closed pools, as the
  // interpreter evaluates them: exact where every point of a smaller hash lies in the closed pool
  // of its rank; a branch no closed pool takes may be taken nowhere, or beyond them
  std::size_t uncertain = 0;
  std::size_t branchCount = 0;
  std::size_t untaken = 0;
  std::size_t uncertainClamps = 0;
  const auto exact = [&](std::size_t found) {
    return merged.scope.closedAll || hashes[found] <= merged.scope.bound;
  };
  if (n > 0) {
    std::vector<double> inputs(merged.inputs.size() * n);
    for (std::size_t k = 0; k < n; ++k) {
      for (std::size_t i = 0; i < merged.inputs.size(); ++i) {
        inputs[i * n + k] = points[k].inputs[i];
      }
    }
    std::map<std::size_t, std::vector<std::uint8_t>> branches;
    static_cast<void>(interpretValues(program, grids, inputs, n, &branches));
    for (const auto& [select, taken] : branches) {
      static_cast<void>(select);
      for (const std::uint8_t branch : {std::uint8_t{1}, std::uint8_t{0}}) {
        ++branchCount;
        std::size_t found = n;
        for (const std::size_t k : byHash) {
          if (closed[k] && k < taken.size() && taken[k] == branch) {
            found = k;
            break;
          }
        }
        if (found == n) {
          ++untaken;
        } else {
          chosen[found] = true;
          uncertain += exact(found) ? 0 : 1;
        }
      }
    }
    const auto clamps = clampsAt(program, grids, inputs, n);
    for (const auto& grid : clamps) {
      int taken = 0;
      for (const std::size_t k : byHash) {
        bool beyond = false;
        for (const auto& axis : grid.axes) {
          beyond = beyond || (k < axis.size() && axis[k] != 0);
        }
        if (closed[k] && beyond && taken < 2) {
          chosen[k] = true;
          ++taken;
          uncertainClamps += exact(k) ? 0 : 1;
        }
      }
    }
  }
  std::size_t unstableChosen = 0;
  for (const std::size_t k : byHash) {
    if (unstable.count(points[k].identity) != 0 && unstableChosen++ < SmallestHashes) {
      chosen[k] = true;
    }
  }

  // the sample, by rank and point
  for (std::size_t k = 0; k < n; ++k) {
    if (chosen[k]) {
      merged.sample.push_back(points[k]);
    }
  }
  std::sort(merged.sample.begin(), merged.sample.end(), [](const auto& a, const auto& b) {
    return std::make_pair(a.rank, a.point) < std::make_pair(b.rank, b.point);
  });
  for (const auto& record : records) {
    for (const auto& entry : record.unstable) {
      if (std::none_of(merged.unstable.begin(), merged.unstable.end(), [&](const auto& kept) {
            return kept.second == entry.second;
          })) {
        merged.unstable.push_back(entry);
      }
    }
  }

  if (!checks.empty()) {
    const bool same = std::all_of(
        checks.begin(), checks.end(), [&](const std::string& check) { return check == checks[0]; });
    merged.head.emplace_back("check", same ? checks[0] : "differs between the ranks");
    if (!same) {
      for (std::size_t r = 0; r < records.size(); ++r) {
        for (const auto& [key, value] : records[r].head) {
          if (key == "check") {
            merged.head.emplace_back("check of " + rankName(records[r]), value);
          }
        }
      }
    }
  }
  std::string selects =
      uncertain == 0
          ? std::string("exact")
          : std::to_string(uncertain) + " branches uncertain, beyond the closed pool of a rank";
  if (untaken > 0) {
    selects += "; " + std::to_string(untaken) + " of " + std::to_string(branchCount) +
               " branches taken in no closed pool";
  }
  merged.head.emplace_back("selects", selects);
  if (uncertainClamps > 0) {
    merged.head.emplace_back("clamped points",
                             std::to_string(uncertainClamps) +
                                 " uncertain, beyond the closed pool of a rank");
  }
  if (!merged.scope.scanned) {
    merged.head.emplace_back("groups and extremes", "uncertain, a rank chose them from its pool");
  }
  if (duplicates > 0) {
    merged.head.emplace_back("duplicates",
                             std::to_string(duplicates) + " points sampled on two ranks");
  }
  if (inconsistent > 0) {
    merged.head.emplace_back("inconsistent",
                             std::to_string(inconsistent) +
                                 " points with other inputs or values on another rank; the "
                                 "first, " +
                                 firstInconsistency);
  }
  return merged;
}

Replay replayRecord(const Record& record, reader::datafield::GridStore& grids) {
  const Program program = programOf(record);
  const std::size_t n = record.sample.size();
  const std::size_t inputCount = record.inputs.size();
  std::vector<double> inputs(inputCount * n);
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < inputCount; ++i) {
      inputs[i * n + k] = record.sample[k].inputs[i];
    }
  }
  Replay replay;
  replay.points = n;
  if (n == 0) {
    return replay;
  }
  const auto values = interpretValues(program, grids, inputs, n);
  constexpr std::size_t Described = 8;
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t o = 0; o < record.outputs.size(); ++o) {
      const double now = values[o * n + k];
      const double then = record.sample[k].interpreter[o];
      if (sameBits(now, then)) {
        continue;
      }
      if (replay.differing++ < Described) {
        replay.differences.push_back("`" + record.outputs[o] + "` at the point of identity " +
                                     std::to_string(record.sample[k].identity) + ": " + hexOf(now) +
                                     ", recorded " + hexOf(then));
      }
    }
  }
  return replay;
}

} // namespace seissol::expr
