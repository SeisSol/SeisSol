// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCRIPTSCOPE_H_
#define SEISSOL_SRC_PHYSICS_SCRIPTSCOPE_H_

#include "Common/ConfigRegistry.h"
#include "Expr/Program.h"
#include "Reader/Scripting/DataReader.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace seissol::physics {

/// The names of the components of the state of the configuration `config`: its quantities, then
/// the memory variables of its mechanisms (model::stateNames).
std::vector<std::string> stateNamesOf(ConfigId config);

/// The names of the quantities of the material of the configuration `config`.
std::vector<std::string> quantityNamesOf(ConfigId config);

/// The names of the parameters of the material of the configuration `config`, as scripts read
/// them (model::scriptParameterNames).
std::vector<std::string> parameterNamesOf(ConfigId config);

/// Whether the configuration `config` knows the name `name`, as a consumer of a script defines it.
using KnowsName = std::function<bool(ConfigId config, const std::string& name)>;

/// Whether a configuration of `configs` has `name` among the names `namesOf` lists for it.
KnowsName knowsAmong(const std::vector<ConfigId>& configs,
                     const std::function<std::vector<std::string>(ConfigId)>& namesOf);

/// A transformation of the program of a script, which its consumer applies after compiling it.
using ProgramFilter = std::function<expr::Program(const expr::Program&)>;

/// A script that serves all configurations of a run -- an initial condition, a boundary condition,
/// an output -- restricted to each of them.
///
/// The configurations of a run may compute different materials, which have names of their own: an
/// elastic material has the stresses s_xx, ..., where an acoustic one has the pressure pprime, and
/// a viscoelastic one has the parameters Qp and Qs, which an elastic one has not. A script that
/// serves them all gives and reads the names of each. In a configuration, an output is left out
/// that is a name only other configurations of the run take, or that reads a name only other
/// configurations provide, directly or through a state (expr::restrictProgram); the consumer takes
/// what is left out as not given. Any other name stays, also one that no configuration knows, for
/// the consumer to report as it does in a run of one configuration -- which therefore runs the
/// script as it is.
class ScriptScope {
  public:
  /// The configurations `configs` of a run, which provide the inputs `provides` says, and take the
  /// outputs `takes` says.
  ScriptScope(std::vector<ConfigId> configs, KnowsName provides, KnowsName takes);

  /// Whether the script keeps what reads the input `name` in the configuration `config`.
  [[nodiscard]] bool keepsInput(ConfigId config, const std::string& name) const;

  /// Whether the script keeps the output `name` in the configuration `config`.
  [[nodiscard]] bool keepsOutput(ConfigId config, const std::string& name) const;

  /// `program`, the script `name`, restricted to the configuration `config`; the log says what is
  /// left out, if anything and if `log`.
  [[nodiscard]] expr::Program restrict(const expr::Program& program,
                                       ConfigId config,
                                       const std::string& name,
                                       bool log = true) const;

  /// restrict() for the configuration `config`, as the filter of the program of the script `name`.
  [[nodiscard]] ProgramFilter filter(ConfigId config, std::string name, bool log = true) const;

  /// The reader of the script at `path` (reader::scripting::buildReader with `defaultInArgs`) for
  /// the configuration `config`: of the restricted program if the script compiles and something is
  /// left out, else the reader as it is built.
  [[nodiscard]] std::unique_ptr<reader::scripting::DataReader>
      reader(const std::string& path,
             const std::vector<std::string>& defaultInArgs,
             ConfigId config) const;

  /// Stops with an error for an output of `program`, the script `name`, that is left out in every
  /// configuration of the run: one that reads the names of several materials, say.
  void checkKept(const expr::Program& program, const std::string& name) const;

  private:
  std::vector<ConfigId> configs_;
  KnowsName provides_;
  KnowsName takes_;
};

/// The scope of a field of the state, as the initial condition and the analytical boundary are: it
/// reads the parameters of the material of a configuration and gives the components of its state.
ScriptScope fieldScope(const std::vector<ConfigId>& configs);

/// The scope of a ghost state given by the state inside, as the nonlinear Dirichlet boundary is: it
/// reads the quantities and the parameters of the material of a configuration and gives its
/// quantities.
ScriptScope ghostStateScope(const std::vector<ConfigId>& configs);

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_SCRIPTSCOPE_H_
