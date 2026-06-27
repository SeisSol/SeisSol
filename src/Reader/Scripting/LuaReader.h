// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_LUAREADER_H_
#define SEISSOL_SRC_READER_SCRIPTING_LUAREADER_H_

#include "Reader/Scripting/DataReader.h"

#include <memory>

extern "C" {
struct lua_State;
}

namespace seissol::reader::scripting {

/// One data field (e.g. ASAGI) declared by the model. Populated from M.field_specs.
struct FieldSpec {
  ///< Internal field handle, e.g. "field_0001"
  std::string name;

  ///< currently "asagi" or "scec". For new types, add here.
  std::string kind;

  ///< Path to the file. Relative to the script.
  std::string file;

  ///< interpolation mode; either "linear" or "nearest"
  std::string interpolation;

  ///< Output names in order (ignored for SCEC)
  std::vector<std::string> parameters;
};

/// Signature of a user-supplied field sampler.
///   out:    pointer to the output buffer (write `parameters.size()` values)
///   coords: pointer to the coordinate vector at which to sample
using FieldSampler = std::function<void(double* out, const double* coords)>;

struct LuaStateState {
  lua_State* luaState{nullptr};
  int32_t refModule{-1};
  int32_t refEvaluate{-1};
  int32_t refFields{-1};
};

/**
    A Lua script reader.
 */
class LuaReader : public DataReader {
  public:
  ~LuaReader() override;

  // Non-copyable, movable
  LuaReader(const LuaReader&) = delete;
  LuaReader& operator=(const LuaReader&) = delete;
  LuaReader(LuaReader&& other) noexcept;
  LuaReader& operator=(LuaReader&& other) noexcept;

  explicit LuaReader(const std::string& code);

  const std::vector<std::string>& inputVars() override { return inputs_; }

  const std::vector<std::string>& outputVars() override { return outputs_; }

  void prepare() override {}

  /**
    The actual call. Should be optimized for dispatches in SIMD environments.
    */
  void call(const scripting::DataTable& table) override;

  /// Internal field closure, public so the extern "C" Lua trampoline can
  /// see the type. Not part of the user-facing API.
  struct FieldClosure {
    LuaReader* owner{}; // back-pointer for context (currently unused)
    FieldSampler sampler;
    std::size_t nOutputs{};
  };

  private:
  void readMetadata(const LuaStateState& luaState);
  void bindField(const std::string& name, FieldSampler sampler);
  void registerFields(const LuaStateState& luaState);

  std::vector<std::string> inputs_;
  std::vector<std::string> outputs_;

  std::string code_;

  std::string version_;
  std::string sourceFile_;

  // state per thread
  std::vector<LuaStateState> luaStates_;

  // constant module data
  std::vector<FieldSpec> fieldSpecs_;
  std::unordered_map<std::string, std::unique_ptr<FieldClosure>> fieldClosures_;
};

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_LUAREADER_H_
