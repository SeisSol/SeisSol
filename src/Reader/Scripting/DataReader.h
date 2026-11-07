// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_DATAREADER_H_
#define SEISSOL_SRC_READER_SCRIPTING_DATAREADER_H_

#include "Reader/Scripting/DataTable.h"

#include <string>
#include <vector>
namespace seissol::reader::scripting {

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
  virtual void call(const scripting::DataTable& table) = 0;
};

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_DATAREADER_H_
