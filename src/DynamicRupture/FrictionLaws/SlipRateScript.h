// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SLIPRATESCRIPT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SLIPRATESCRIPT_H_

#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace seissol::dr::friction_law {

/// The slip of the imposed slip rates of a script (FL 36).
///
/// The script reads the time `t`, the point `x`, `y`, `z`, the simulation `sim`, and any other
/// name as a parameter of the fault, read from the fault parameter file (as FL 33 reads
/// `strike_slip`). It gives either the slip at `t`, `slip_strike` and `slip_dip`, or the slip rate
/// at `t`, `slip_rate_strike` and `slip_rate_dip`: along strike and dip, with the signs of FL 33.
/// A script of the rate may also read `dt`, the length of the sub-step that ends at `t`.
///
/// The rate of a sub-step is taken from the slip as its increment over the sub-step divided by
/// the length of the sub-step, so the slip SeisSol accumulates is that of the script at every end
/// of a sub-step, whatever the shape of the function. A rate is taken as given at the end of the
/// sub-step, as FL 33 does.
///
/// Either is turned into one program here, whose inputs are rows of values per point (`rows()`,
/// filled by the initializer) and the times of a sub-step, the same at all points (passed per
/// call), and whose outputs are the slip rates along the two directions of the fault coordinate
/// system of the face (`Output1`, `Output2`).
class SlipRateScript {
  public:
  /// What fills a row.
  enum class Row : std::uint8_t {
    /// a fault parameter of that name
    Parameter,
    X,
    Y,
    Z,
    /// the simulation of the point
    Simulation,
    /// the rotation from strike and dip into the fault coordinate system of the face
    Cos,
    Sin
  };

  /// The inputs given per call: the end, the length and the start of a sub-step. Names a script
  /// cannot use are those the rewrite introduces.
  static constexpr const char* Time = "t";
  static constexpr const char* TimeStep = "dt";
  static constexpr const char* StartTime = "$t_start";
  static constexpr const char* RotationCos = "$cos";
  static constexpr const char* RotationSin = "$sin";
  static constexpr const char* Output1 = "$slip_rate_1";
  static constexpr const char* Output2 = "$slip_rate_2";

  /// Loads the script at `path` (as `wavefieldscript`: an sderiv module or a Lua model that
  /// traces) and checks it; one that cannot be used is a configuration error.
  explicit SlipRateScript(const std::string& path);

  /// From the program of a script, which `path` only names in messages. Throws
  /// std::invalid_argument for a program that cannot be used.
  SlipRateScript(const expr::Program& script, const std::string& path);

  /// The program of the slip rates.
  [[nodiscard]] const expr::Program& program() const { return program_; }

  /// The rows, in order: the names of their inputs of program(), and what fills them.
  [[nodiscard]] const std::vector<std::string>& rowNames() const { return rowNames_; }
  [[nodiscard]] const std::vector<Row>& rowKinds() const { return rowKinds_; }
  [[nodiscard]] std::size_t rows() const { return rowNames_.size(); }

  /// Whether the script gives the slip (and not the slip rate).
  [[nodiscard]] bool cumulative() const { return cumulative_; }

  [[nodiscard]] const std::string& path() const { return path_; }

  private:
  expr::Program program_;
  std::vector<std::string> rowNames_;
  std::vector<Row> rowKinds_;
  bool cumulative_{true};
  std::string path_;
};

/// Evaluates a slip rate script on the points of one layer of faces, for the sub-steps of a time
/// step, into the slip rates the friction law imposes.
///
/// Both arrays of a layer are laid out row-major over all of its points (face by face, the
/// points of a face in order): row r of the parameters at `parameters + r * points`, and the slip
/// rate along direction d (0, 1) of sub-step i at `slipRates + (2 * i + d) * points`, in the type
/// given to the constructor. A row is thus one column of the program, which a device kernel reads
/// like any other.
///
/// On the host, one kernel per thread evaluates a range of points per call; on a device, one
/// kernel evaluates the whole layer per sub-step, on the stream of the friction law. A build
/// without a runtime compiler for its device evaluates on the host and says so.
class SlipRateEvaluator {
  public:
  struct LayerData {
    const double* parameters{nullptr};
    void* slipRates{nullptr};
  };

  SlipRateEvaluator(std::shared_ptr<const SlipRateScript> script,
                    reader::scripting::DataType slipRateType,
                    std::size_t timeSteps);
  ~SlipRateEvaluator();

  SlipRateEvaluator(const SlipRateEvaluator&) = delete;
  SlipRateEvaluator& operator=(const SlipRateEvaluator&) = delete;
  SlipRateEvaluator(SlipRateEvaluator&&) = delete;
  SlipRateEvaluator& operator=(SlipRateEvaluator&&) = delete;

  /// The layer, with `points` points. `device` holds the addresses on the device if the friction
  /// law runs there; it is not used otherwise.
  void setLayer(std::size_t points, const LayerData& host, const LayerData& device);

  /// The slip rates of the sub-steps of the step that starts at `startTime`, with the sub-steps
  /// `deltaT` long. With `onDevice`, they are written on the device on `stream` if it can, which
  /// it returns; false means they were written on the host, and the caller copies them.
  bool evaluate(double startTime, const std::vector<double>& deltaT, void* stream, bool onDevice);

  private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SLIPRATESCRIPT_H_
