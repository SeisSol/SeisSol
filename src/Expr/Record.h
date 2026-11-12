// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_RECORD_H_
#define SEISSOL_SRC_EXPR_RECORD_H_

// The evaluation record of a program (docs/sderiv-spec.rst, "Evaluation record"): what SeisSol
// writes at the first call of a kernel, one file per program, backend and rank
// (setRecordDirectory()), and what the command line tool reads, merges over the ranks and
// evaluates again.
//
//   # The evaluation record of a program ...
//   record = 1
//   fingerprint = 3c351edd263b9464
//   backend = rtc-cpu
//   rank = 0
//   call = 10000 points
//   identities = global
//   pool = 4096 points of the smallest hashes, up to 0x68f2c3a7d2b1e0f4
//   scan = all points
//   compared = 4103 points
//   compute = f64
//   ...                                  the build, as SeisSol gives it
//   check = the same bits
//   clamped m = 0 12 0, 0 0 0, 31 0 0 of 4096 points
//   grid m = path "vmodel.nc", sha256 "5a1c...", variable "data", components "rho" "mu"
//   grid m axis 0 = dimension "x", units "m", min -1000, spacing 100, nodes 21
//
//   # the program, in the normal form of the core format
//   #!sderiv 1
//   ...
//
//   # the sample: ...
//   point identity x group u:interpreter u:backend
//   12 9876543 0x1.8p-2 0x1p+0 0x1.8p-2 0x1.8p-2
//
//   # unstable points, ...: the point and its identity
//   12 9876543
//
// A merged record has `merged = 0 1 2`, the ranks it merges, in place of `rank`, and the rank of
// each point of its sample in a first column.

#include "Expr/Backend.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace seissol::expr {

/// The newest version of the record that this build writes and reads; it reads every older one.
inline constexpr int RecordVersion = 1;

struct Record {
  std::string fingerprint;
  std::string backend;
  /// the rank that wrote the record; none for a merged one, which names its ranks
  std::optional<int> rank;
  std::vector<int> ranks;
  PoolScope scope;
  std::size_t compared{0};
  /// the compute type: "f64" or "f32"
  std::string compute{"f64"};
  /// the other lines of the head, in their order: the build, the check, what a merge found
  std::vector<std::pair<std::string, std::string>> head;
  /// per grid, the points of the closed pool it samples beyond its volume
  std::vector<ClampCount> clamps;
  /// the grids of the program: files, hashes, geometry
  std::vector<GridRecord> grids;
  /// the program in the normal form of the core format; empty for one the format cannot express,
  /// with the reason in `notCore`
  std::string program;
  std::string notCore;

  std::vector<std::string> inputs;
  std::vector<std::string> outputs;
  /// whether the sample has the backend's values (for a compiled backend), and their tolerances
  /// (against the perturbed reference)
  bool backendValues{false};
  bool tolerances{false};

  struct Point {
    int rank{0};
    std::size_t point{0};
    std::uint64_t identity{0};
    std::vector<double> inputs;
    std::vector<double> interpreter;
    std::vector<double> backend;
    std::vector<double> tolerance;
  };
  std::vector<Point> sample;
  /// where a perturbed evaluation decides a select otherwise: the point and its identity
  std::vector<std::pair<std::size_t, std::uint64_t>> unstable;

  /// What the backend wrote at the point `k` of the sample for output `o`: the interpreter's value
  /// for the interpreter.
  [[nodiscard]] double written(std::size_t k, std::size_t o) const {
    return backendValues ? sample[k].backend[o] : sample[k].interpreter[o];
  }
};

/// The record of the check of the first call of a kernel of `program` on `backend`, written by
/// `rank`, with the lines `build` of the build ("key = value").
[[nodiscard]] Record recordOf(const Program& program,
                              const std::string& backend,
                              const BackendComparison& found,
                              const std::vector<std::string>& build,
                              int rank);

/// The text of `record`, which parseRecord() reads back.
[[nodiscard]] std::string printRecord(const Record& record);

/// The record of the text `text`; `file` names it in errors. Throws std::runtime_error.
[[nodiscard]] Record parseRecord(const std::string& text, const std::string& file = "");

/// The program of `record`, in its compute type. Throws std::runtime_error where it has none, or
/// where the program does not have the fingerprint of the record.
[[nodiscard]] Program programOf(const Record& record);

/// The record of the sample of all points, from the records of the ranks of one program on one
/// backend (the specification's "Evaluation record"): the rules of the sample over the union of
/// their samples, a point evaluated on two ranks (a ghost cell) once. Exact for the rules of the
/// smallest hash, of the groups and of the extreme values, where every rank scanned all points of
/// its call; for the selects, where the point of the smallest hash that takes a branch lies in the
/// closed pool of every rank (else the head says `selects = ...`). The interpreter evaluates the
/// program at the points of the union for the branches of its selects, with the grids of `grids`.
/// Throws std::runtime_error for records of different programs or backends, or of points without
/// identities that hold across the partition.
[[nodiscard]] Record mergeRecords(const std::vector<Record>& records,
                                  reader::datafield::GridStore& grids);

/// What the interpreter of this build gives at the points of a sample, against what the record
/// says the interpreter gave there.
struct Replay {
  std::size_t points{0};
  /// the values that differ in their bits (NaN equals NaN)
  std::size_t differing{0};
  /// the first of them, described
  std::vector<std::string> differences;
};

/// Evaluates the program of `record` with the interpreter at the inputs of its sample, with the
/// grids of `grids`, and compares bit for bit (the specification's "Evaluation record": what a
/// change of the interpreter or of the maths library moves).
[[nodiscard]] Replay replayRecord(const Record& record, reader::datafield::GridStore& grids);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_RECORD_H_
