// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_READER_SCRIPTING_EASIREADER_H_
#define SEISSOL_SRC_READER_SCRIPTING_EASIREADER_H_

#include "Reader/Scripting/DataReader.h"

#include <memory>

namespace easi {
class Component;
class AsagiReader;
class YAMLParser;
} // namespace easi

namespace seissol::reader::scripting {

/**
    Reads data through an easi file.
 */
class EasiReader : public DataReader {
  public:
  ~EasiReader() override;

  EasiReader(const std::string& script, const std::vector<std::string>& inVars);

  const std::vector<std::string>& inputVars() override { return inVars_; }

  const std::vector<std::string>& outputVars() override { return outVars_; }

  void prepare() override {}

  /**
      The actual call. Should be optimized for dispatches in SIMD environments.
      */
  void call(const scripting::DataTable& table) override;

  private:
  std::unique_ptr<easi::Component> components_;
  std::unique_ptr<easi::AsagiReader> asagiReader_;
  std::unique_ptr<easi::YAMLParser> parser_;

  std::vector<std::string> inVars_;
  std::vector<std::string> outVars_;

  std::string script_;
};

} // namespace seissol::reader::scripting
#endif // SEISSOL_SRC_READER_SCRIPTING_EASIREADER_H_
