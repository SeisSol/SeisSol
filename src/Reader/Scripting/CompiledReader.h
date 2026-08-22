// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_COMPILEDREADER_H_
#define SEISSOL_SRC_READER_SCRIPTING_COMPILEDREADER_H_

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataReader.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

struct CompiledReaderOptions {
  expr::BackendOptions backend;

  /// Points COMPARED against the reference reader in prepare(), not points
  /// evaluated: the interpreted reader has no subset entry point, so the check
  /// costs one full reference pass whatever this is set to. 0 disables it.
  std::size_t differentialSamples{64};

  /// Relative deviation above which the compiled kernel is distrusted.
  ///
  /// Not bit equality, even though LuaTracer.t.h uses that bar for the same
  /// comparison. That test pins one libm on one machine; this is a RUNTIME gate,
  /// and making a production run fall back to the slow path because a vendor
  /// libm rounds pow() one ulp differently would be a performance cliff nobody
  /// could diagnose from a log. The worst observed deviation is logged either
  /// way, so a real lowering bug is visible long before it reaches this bound.
  double differentialTolerance{1e-10};
};

/// A DataReader that evaluates a compiled expr::Program.
///
/// Frontend-agnostic on purpose: it takes a Program, not a script. Tracing or
/// compiling happens in buildReader, because the decision of what to do when a
/// frontend fails -- fall back, or refuse -- belongs to whoever knows which
/// frontend it was.
class CompiledReader : public DataReader {
  public:
  /// `reference` is both the differential oracle and the fallback, and it may
  /// be null. One object rather than two because they must agree: falling back
  /// to a reader other than the one that was compared is a leap of faith.
  ///
  /// Null is the sderiv case -- there is no interpreter for that language, so
  /// there is nothing to compare against and nothing to fall back TO. A failure
  /// there is a configuration error, which is what makes the frontend's parse
  /// diagnostics load-bearing rather than a convenience.
  CompiledReader(expr::Program program,
                 std::unique_ptr<DataReader> reference,
                 CompiledReaderOptions options = {},
                 datafield::GridStore* grids = nullptr);

  ~CompiledReader() override;

  CompiledReader(const CompiledReader&) = delete;
  CompiledReader& operator=(const CompiledReader&) = delete;
  CompiledReader(CompiledReader&&) = delete;
  CompiledReader& operator=(CompiledReader&&) = delete;

  const std::vector<std::string>& inputVars() override { return inputs_; }
  const std::vector<std::string>& outputVars() override { return outputs_; }

  void prepare(const DataTable& table) override;
  void call(const DataTable& table) override;

  /// Whether the differential check rejected the kernel and calls are being
  /// served by the reference reader. For tests and for a status line.
  [[nodiscard]] bool usingFallback() const { return fallback_; }

  [[nodiscard]] const expr::Program& program() const { return program_; }

  private:
  /// Runs the reference reader, then the kernel, and compares a sample.
  /// Returns false when the kernel is not to be trusted.
  bool differentialCheck(const DataTable& table);

  expr::Program program_;
  std::unique_ptr<DataReader> reference_;
  CompiledReaderOptions options_;
  std::vector<std::string> inputs_;
  std::vector<std::string> outputs_;

  datafield::GridStore* grids_;
  expr::Binding binding_;
  std::unique_ptr<expr::Kernel> kernel_;

  /// Identity of the table the binding resolves against. A binding is a set of
  /// column indices into ONE table; using it with another is how a program
  /// silently reads the wrong column, which is the failure this whole layer is
  /// built to prevent.
  const DataTable* preparedFor_{nullptr};
  bool fallback_{false};
};

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_COMPILEDREADER_H_
