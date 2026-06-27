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

  virtual void prepare() = 0;

  /**
      The actual call. Should be optimized for dispatches in SIMD environments.
      */
  virtual void call(const scripting::DataTable& table) = 0;
};

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_DATAREADER_H_
