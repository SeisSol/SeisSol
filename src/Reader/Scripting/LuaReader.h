// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_LUAREADER_H_
#define SEISSOL_SRC_READER_SCRIPTING_LUAREADER_H_

#include "Reader/Scripting/DataReader.h"

extern "C" {
struct lua_State;
}

namespace seissol::reader::scripting {

/**
    A Lua script reader.
 */
class LuaReader : public DataReader {
  public:
  ~LuaReader() override;

  explicit LuaReader(const std::string& code);

  const std::vector<std::string>& inputVars() override { return inVars_; }

  const std::vector<std::string>& outputVars() override { return outVars_; }

  void prepare() override {}

  /**
    The actual call. Should be optimized for dispatches in SIMD environments.
    */
  void call(const scripting::DataTable& table) override;

  private:
  std::vector<std::string> inVars_;
  std::vector<std::string> outVars_;

  std::string code_;

  std::vector<lua_State*> luaStates_;
};

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_LUAREADER_H_
