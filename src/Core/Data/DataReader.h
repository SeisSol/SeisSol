// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_CORE_DATA_DATAREADER_H_
#define SEISSOL_SRC_CORE_DATA_DATAREADER_H_

#include "Core/Data/DataTable.h"

#include <exception>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::data {

/**
    A data source.
    Anything data that's piped in externally.
    */
class DataReader {
  public:
  virtual ~DataReader() = default;

  virtual const std::vector<std::string>& inputVars() = 0;
  virtual const std::vector<std::string>& outputVars() = 0;

  /// Resolve everything that depends on `table` but not on its values: column
  /// lookup, direction and type checks, and for a compiled reader the binding,
  /// the kernel and its precompute stage.
  ///
  /// Calling it is optional: call() prepares itself for a table it has not seen.
  /// It exists so a consumer that evaluates the same table repeatedly can move
  /// the cost out of the first timestep, which is what the analytic boundary
  /// condition will want.
  virtual void prepare(const DataTable& table) = 0;

  /**
      The actual call. Should be optimized for dispatches in SIMD environments.
      */
  virtual void call(const DataTable& table) = 0;
};

/// Opens the model in the file `path`, which reads the columns `inputs` (e.g. x, y and z). The
/// composition decides which readers there are; a part that reads a model is handed one of these,
/// and does not know the formats.
using ModelLoader = std::function<std::unique_ptr<DataReader>(
    const std::string& path, const std::vector<std::string>& inputs)>;

/// The columns that `model` supplies.
inline std::set<std::string> suppliedParameters(DataReader& model) {
  const auto& outputs = model.outputVars();
  return {outputs.begin(), outputs.end()};
}

/// Evaluates `model` for the points of `table`; an error stops the run, naming `hint`.
inline void evaluate(DataReader& model, const DataTable& table, const std::string& hint) {
  try {
    model.call(table);
  } catch (const std::exception& error) {
    logError() << "Error while evaluating a model for" << hint << ":" << std::string(error.what());
  }
}

} // namespace seissol::data

#endif // SEISSOL_SRC_CORE_DATA_DATAREADER_H_
