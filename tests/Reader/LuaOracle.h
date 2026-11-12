// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_TESTS_READER_LUAORACLE_H_
#define SEISSOL_TESTS_READER_LUAORACLE_H_

#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

extern "C" {
struct lua_State;
}

namespace seissol::unit_test {

struct LuaStateState {
  lua_State* luaState{nullptr};
  int32_t refModule{-1};
  int32_t refEvaluate{-1};
  int32_t refFields{-1};
};

/// The interpreted Lua of the tests: a Lua model (a `lua:` module) run by Lua 5.4 with lmathx, as
/// easi and SeisSol load it, point by point -- the interpreted reader SeisSol had, kept as the
/// oracle the Lua frontend (LuaFrontend.h) is compared with, bit for bit. It has no grids.
class LuaOracle : public reader::scripting::DataReader {
  public:
  ~LuaOracle() override;

  LuaOracle(const LuaOracle&) = delete;
  LuaOracle& operator=(const LuaOracle&) = delete;
  LuaOracle(LuaOracle&& other) noexcept = delete;
  LuaOracle& operator=(LuaOracle&& other) noexcept = delete;

  explicit LuaOracle(const std::string& code);

  const std::vector<std::string>& inputVars() override { return inputs_; }

  const std::vector<std::string>& outputVars() override { return outputs_; }

  void prepare(const reader::scripting::DataTable& /*table*/) override {}

  void call(const reader::scripting::DataTable& table) override;

  private:
  /// An argument of evaluate after `fields`: an input column, or the value of a state.
  struct Argument {
    bool state{false};
    std::size_t index{0};
  };

  void readMetadata(const LuaStateState& luaState);

  std::vector<std::string> inputs_;
  std::vector<std::string> outputs_;

  std::vector<Argument> arguments_;

  // M.state, sorted by name, and the value of every state at every point of the last call's
  // point set; a point set of another size starts over from the initial values.
  std::vector<std::string> stateNames_;
  std::vector<double> stateInitial_;
  std::vector<std::vector<double>> stateValues_;
  std::size_t statePoints_{0};

  std::string code_;

  std::string version_;
  std::string sourceFile_;

  // state per thread
  std::vector<LuaStateState> luaStates_;
};

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_READER_LUAORACLE_H_
