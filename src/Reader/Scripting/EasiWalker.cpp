// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "EasiWalker.h"

#include "Expr/Core.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/AsagiLite.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/ScecGrid.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/EasiSpecial.h"
#include "Reader/Scripting/ImpalaToLua.h"
#include "Reader/Scripting/LuaFrontend.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::reader::scripting {

namespace {

namespace fs = std::filesystem;
using expr::Fn;
using expr::NodeId;
using expr::NoNode;

/// A value a component gives: its node, and where in the files the components are that compute it.
struct Value {
  NodeId node{NoNode};
  std::vector<expr::SourcePosition> from;
};

/// The values of the dimensions or parameters of a component, by name: in the order of easi's sets.
using Values = std::map<std::string, Value>;

std::string listOf(const std::set<std::string>& names) {
  std::string text;
  for (const auto& name : names) {
    text += (text.empty() ? "" : ", ") + name;
  }
  return text;
}

std::set<std::string> keysOf(const Values& values) {
  std::set<std::string> keys;
  for (const auto& [name, value] : values) {
    keys.insert(name);
  }
  return keys;
}

class Walker {
  public:
  explicit Walker(EasiQuery query) : query_(std::move(query)) {
    for (const auto& dimension : query_.dimensions) {
      rootDimensions_.insert(dimension);
    }
  }

  expr::Program translate(const std::string& path) {
    const std::string file = loadPath(path, "");
    const YAML::Node root = load(file);
    program_.addMeta(expr::MetaEntry{"generator", {expr::MetaValue{"seissol easi walker", true}}});
    // the dimensions, in the order of easi's set, then the group where a filter reads it
    Values in;
    for (const auto& dimension : rootDimensions_) {
      program_.addInput(dimension, DataType::F64);
      in[dimension] = Value{arena().field(dimension), {}};
    }
    const std::set<std::string> wanted(query_.parameters.begin(), query_.parameters.end());
    files_.push_back(file);
    const Values values = evaluate(root, in, wanted, NoNode);
    files_.pop_back();
    for (const auto& parameter : query_.parameters) {
      const auto found = values.find(parameter);
      const NodeId node = found != values.end() && found->second.node != NoNode
                              ? found->second.node
                              : arena().konst(std::numeric_limits<double>::quiet_NaN());
      program_.addOutput(parameter, DataType::F64, node);
      if (found != values.end()) {
        for (const auto& position : found->second.from) {
          program_.addSource(
              expr::SourceOf{expr::SourceOf::What::Output, NoNode, parameter, position});
        }
      }
    }
    for (auto& [condition, message] : checks_) {
      program_.addCheck(condition, message);
    }
    return std::move(program_);
  }

  std::set<std::string> parameters(const std::string& path) {
    const std::string file = loadPath(path, "");
    files_.push_back(file);
    const auto supplied = suppliedBy(load(file), rootDimensions_);
    files_.pop_back();
    return supplied;
  }

  private:
  expr::Arena& arena() { return program_.arena(); }
  NodeId konst(double value) { return arena().konst(value); }
  NodeId pw(Fn fn, NodeId a, NodeId b) { return arena().pw(fn, a, b); }
  NodeId pw(Fn fn, NodeId a, NodeId b, NodeId c) { return arena().pw(fn, a, b, c); }

  /// `a and b`, where NoNode is true.
  NodeId both(NodeId a, NodeId b) {
    if (a == NoNode) {
      return b;
    }
    if (b == NoNode) {
      return a;
    }
    return pw(Fn::And, a, b);
  }
  /// `a or b`, where NoNode is true.
  NodeId either(NodeId a, NodeId b) {
    if (a == NoNode || b == NoNode) {
      return NoNode;
    }
    return pw(Fn::Or, a, b);
  }
  NodeId negation(NodeId a) { return pw(Fn::Eq, a, konst(0.0)); }

  [[nodiscard]] const std::string& file() const { return files_.back(); }

  /// easi's reference to a node: the file and the line.
  [[nodiscard]] std::string where(const YAML::Node& node) const {
    return file() + "@" + std::to_string(node.Mark().line + 1);
  }

  [[nodiscard]] expr::SourcePosition positionOf(const YAML::Node& node) const {
    return expr::SourcePosition{file(), node.Mark().line + 1, node.Mark().column + 1};
  }

  [[noreturn]] void fail(const YAML::Node& node, const std::string& what) const {
    throw EasiModelError(where(node) + ": " + what);
  }

  /// That every point which reaches here (`reach`, NoNode for every point) meets `condition`:
  /// a check of the program, with easi's message where it does not.
  void require(NodeId reach, NodeId condition, const std::string& message) {
    const NodeId holds = reach == NoNode ? condition : pw(Fn::Or, negation(reach), condition);
    checks_.emplace_back(holds, message);
  }

  /// The file `name` an !Include names, as easi finds it: next to the file that includes it, else
  /// as it is given.
  [[nodiscard]] std::string loadPath(const std::string& name, const std::string& current) const {
    const fs::path next(name);
    const fs::path next0 =
        next.is_relative() && !current.empty() ? fs::path(current).parent_path() / next : next;
    std::error_code error;
    if (fs::exists(next0, error)) {
      return fs::canonical(next0, error).string();
    }
    if (fs::exists(next, error)) {
      return fs::canonical(next, error).string();
    }
    return next.string();
  }

  /// The data file `name` of a component (ASAGI), as easi finds it: next to the current file,
  /// where it is there, else as it is given.
  [[nodiscard]] std::string dataPath(const std::string& name) const {
    const fs::path entry(name);
    if (entry.is_relative()) {
      const fs::path next = fs::path(file()).parent_path() / entry;
      std::error_code error;
      if (fs::exists(next, error)) {
        return fs::canonical(next, error).string();
      }
    }
    return entry.string();
  }

  YAML::Node load(const std::string& path) const {
    try {
      return YAML::LoadFile(path);
    } catch (const YAML::Exception& error) {
      throw EasiModelError(path + ": " + error.what());
    }
  }

  template <typename T>
  T as(const YAML::Node& node, const char* key) const {
    const YAML::Node value = node[key];
    if (!value) {
      fail(node, std::string("the component has no `") + key + "`");
    }
    try {
      return value.as<T>();
    } catch (const YAML::Exception& error) {
      fail(value, std::string("`") + key + "` is not of the form easi reads: " + error.what());
    }
  }

  /// The components of a composite: its `components`, one or a list.
  static std::vector<YAML::Node> childrenOf(const YAML::Node& node) {
    std::vector<YAML::Node> children;
    if (node.IsMap() && node["components"]) {
      const YAML::Node components = node["components"];
      if (components.IsSequence()) {
        for (const auto& child : components) {
          children.push_back(child);
        }
      } else {
        children.push_back(components);
      }
    }
    return children;
  }

  /// The root of the file an !Include names, and the file.
  std::pair<YAML::Node, std::string> included(const YAML::Node& node) {
    const std::string path = loadPath(node.as<std::string>(), file());
    return {load(path), path};
  }

  /// The values a map gives at `in`, its own outputs, before its components: those of a
  /// ConstantMap, an AffineMap, a PolynomialMap, an EvalModel and an ASAGI grid.
  Values mapped(const YAML::Node& node, const Values& in, NodeId reach) {
    const std::string tag = node.Tag();
    Values out;
    if (tag == "!ConstantMap") {
      for (const auto& [name, value] : as<std::map<std::string, double>>(node, "map")) {
        out[name] = Value{konst(value), {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!AffineMap") {
      const auto matrix = as<std::map<std::string, std::vector<double>>>(node, "matrix");
      const auto translation = as<std::map<std::string, double>>(node, "translation");
      if (matrix.size() != translation.size()) {
        fail(node, "Matrix and translation must have the same size in an affine map.");
      }
      // easi's product x * M from zero, a dimension after the other, then the translation; a
      // translation is paired with the row of the same position
      auto shift = translation.begin();
      for (const auto& [name, row] : matrix) {
        if (row.size() != in.size()) {
          fail(node, "Number of matrix entries does not match number of input parameters.");
        }
        NodeId value = konst(0.0);
        std::size_t k = 0;
        for (const auto& [dimension, coordinate] : in) {
          value = pw(Fn::Add, value, pw(Fn::Mul, coordinate.node, konst(row[k++])));
        }
        out[name] = Value{pw(Fn::Add, value, konst((shift++)->second)), {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!PolynomialMap") {
      if (in.size() != 1) {
        fail(node, "Polynomial map requires 1D input.");
      }
      const NodeId x = in.begin()->second.node;
      const auto coefficients = as<std::map<std::string, std::vector<double>>>(node, "map");
      std::size_t count = std::numeric_limits<std::size_t>::max();
      for (const auto& [name, values] : coefficients) {
        if (count != std::numeric_limits<std::size_t>::max() && values.size() != count) {
          fail(node,
               "All parameters in a polynomial map must have the same number of coefficients.");
        }
        count = values.size();
      }
      for (const auto& [name, values] : coefficients) {
        if (values.empty()) {
          fail(node, "a polynomial map has coefficients");
        }
        // Horner's scheme as easi has it: c_k + y * x
        NodeId value = konst(values[0]);
        for (std::size_t k = 1; k < values.size(); ++k) {
          value = pw(Fn::Add, konst(values[k]), pw(Fn::Mul, value, x));
        }
        out[name] = Value{value, {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!EvalModel") {
      if (!node["model"]) {
        fail(node, "the component has no `model`");
      }
      const auto names = as<std::vector<std::string>>(node, "parameters");
      const std::set<std::string> wanted(names.begin(), names.end());
      Values out = evaluate(node["model"], in, wanted, reach);
      // a parameter the model does not give fails a check there; NaN here
      for (const auto& name : wanted) {
        if (out.count(name) == 0) {
          out[name] = Value{konst(std::numeric_limits<double>::quiet_NaN()), {}};
        }
      }
      return out;
    }
    if (tag == "!ASAGI") {
      return grid(node, in);
    }
    if (tag == "!SCECFile") {
      return scecGrid(node, in);
    }
    if (tag == "!FunctionMap") {
      return impalaMapped(node, in, reach);
    }
    if (tag == "!LuaMap") {
      return luaMapped(node, in, reach);
    }
    if (const auto* special = easiSpecialOf(tag); special != nullptr) {
      return specialMapped(node, *special, in);
    }
    if (tag == "!Include") {
      fail(node, "an !Include is no map");
    }
    fail(node, "Unknown tag " + tag);
  }

  /// A special component (EasiSpecial.h): its inputs from its constants, the others from `in`,
  /// which are to be exactly those, as easi has it.
  Values specialMapped(const YAML::Node& node, const EasiSpecial& special, const Values& in) {
    std::map<std::string, double> constants;
    if (node["constants"]) {
      if (!node["constants"].IsMap()) {
        fail(node["constants"], "`constants` is not of the form easi reads: a map");
      }
      constants = as<std::map<std::string, double>>(node, "constants");
    }
    std::set<std::string> rest(special.in.begin(), special.in.end());
    std::map<std::string, NodeId> inputs;
    for (const auto& [name, value] : constants) {
      if (rest.erase(name) == 0) {
        fail(node, "Unknown constant " + name + ".");
      }
      inputs[name] = konst(value);
    }
    if (rest != keysOf(in)) {
      fail(node,
           "The input parameters of a component ({" + listOf(rest) +
               "}) added to a composite are incompatible (should be {" + listOf(keysOf(in)) +
               "}).");
    }
    for (const auto& name : rest) {
      inputs[name] = in.at(name).node;
    }
    Values out;
    for (const auto& [name, value] : easiSpecialValues(arena(), special, inputs)) {
      out[name] = Value{value, {positionOf(node)}};
    }
    return out;
  }

  /// The indicator of the block scalar `scalar` (`|` or `>`), or 0 for another scalar.
  [[nodiscard]] char blockIndicatorOf(const YAML::Node& scalar) {
    const auto mark = scalar.Mark();
    auto text = texts_.find(file());
    if (text == texts_.end()) {
      std::ifstream stream(file());
      text = texts_
                 .emplace(file(),
                          std::string((std::istreambuf_iterator<char>(stream)),
                                      std::istreambuf_iterator<char>()))
                 .first;
    }
    const bool block = mark.pos >= 0 && static_cast<std::size_t>(mark.pos) < text->second.size() &&
                       (text->second[mark.pos] == '|' || text->second[mark.pos] == '>');
    return block ? text->second[mark.pos] : '\0';
  }

  /// The line of the file the first line of the text of `scalar` is on: for a block scalar (`|`,
  /// `>`), the line after its indicator.
  [[nodiscard]] int firstLineOf(const YAML::Node& scalar) {
    return scalar.Mark().line + (blockIndicatorOf(scalar) != '\0' ? 2 : 1);
  }

  /// The Lua `code` of a component, at `origin`, translated by the Lua frontend as easi runs the
  /// Lua of an !LuaMap (Lua 5.4 with lmathx, the chunk once and `f` at every point): the values
  /// `returns` of `f` at `in`. `what` names the Lua where the frontend refuses it.
  lua::LuaMapTranslation translatedLua(const std::string& code,
                                       const lua::LuaOrigin& origin,
                                       const Values& in,
                                       const std::set<std::string>& returns,
                                       NodeId reach,
                                       const std::string& what) {
    std::vector<std::pair<std::string, NodeId>> inputs;
    for (const auto& [name, value] : in) {
      inputs.emplace_back(name, value.node);
    }
    try {
      auto translation = lua::translateLuaMap(program_, code, inputs, returns, reach, origin);
      for (auto& check : translation.checks) {
        checks_.emplace_back(check.holds, std::move(check.message));
      }
      return translation;
    } catch (const lua::LuaUnsupported& error) {
      throw EasiUnsupported(origin.at(error.line()) + ": the walker does not translate " + what +
                            ": " + error.message() + "; easi evaluates the file");
    } catch (const lua::LuaError& error) {
      throw EasiModelError(origin.at(error.line()) + ": " + error.message());
    }
  }

  /// The values of an !LuaMap at `in`: its function `f`.
  Values luaMapped(const YAML::Node& node, const Values& in, NodeId reach) {
    std::set<std::string> returns;
    const YAML::Node given = node["returns"];
    if (!given) {
      fail(node, "the component has no `returns`");
    }
    try {
      if (given.IsScalar()) {
        returns.insert(given.as<std::string>());
      } else {
        const auto names = given.as<std::vector<std::string>>();
        returns.insert(names.begin(), names.end());
      }
    } catch (const YAML::Exception& error) {
      fail(given, std::string("`returns` is not of the form easi reads: ") + error.what());
    }
    const auto code = as<std::string>(node, "function");
    const lua::LuaOrigin origin{file(), firstLineOf(node["function"]), '@'};
    const auto translation =
        translatedLua(code, origin, in, returns, reach, "the Lua of this !LuaMap");
    std::vector<expr::SourcePosition> from{positionOf(node)};
    for (const int line : translation.lines) {
      from.push_back(expr::SourcePosition{file(), origin.firstLine + line - 1, 1});
    }
    Values out;
    for (const auto& [name, value] : translation.values) {
      out[name] = Value{value, from};
    }
    return out;
  }

  /// The values of a !FunctionMap at `in`: the Lua easi makes of its Impala (ImpalaToLua.h),
  /// translated as that of an !LuaMap, with a warning where easi's grouping changes what a chain
  /// means.
  Values impalaMapped(const YAML::Node& node, const Values& in, NodeId reach) {
    const YAML::Node map = node["map"];
    if (!map) {
      fail(node, "the component has no `map`");
    }
    if (!map.IsMap()) {
      fail(map, "`map` is not of the form easi reads: the code of every parameter");
    }
    std::vector<ImpalaEntry> entries;
    std::map<std::string, expr::SourcePosition> positions;
    for (const auto& entry : map) {
      ImpalaEntry impala;
      try {
        impala.name = entry.first.as<std::string>();
        impala.code = entry.second.as<std::string>();
      } catch (const YAML::Exception& error) {
        fail(entry.first, std::string("`map` is not of the form easi reads: ") + error.what());
      }
      impala.line = firstLineOf(entry.second);
      impala.literal = blockIndicatorOf(entry.second) == '|';
      positions.insert_or_assign(impala.name, positionOf(entry.second));
      entries.push_back(std::move(impala));
    }
    ImpalaLua lua;
    try {
      lua = impalaToLua(entries, keysOf(in));
    } catch (const ImpalaError& error) {
      throw EasiModelError(file() + "@" + std::to_string(error.line()) + ": " + error.what());
    }
    for (const auto& chain : lua.chains) {
      warnGrouping(chain);
    }
    std::set<std::string> returns;
    for (const auto& [name, position] : positions) {
      returns.insert(name);
    }
    const auto translation = translatedLua(lua.code,
                                           lua::LuaOrigin{file(), lua.firstLine, '@'},
                                           in,
                                           returns,
                                           reach,
                                           "the Lua easi makes of the Impala of this !FunctionMap");
    Values out;
    for (const auto& [name, value] : translation.values) {
      out[name] = Value{value, {positionOf(node), positions.at(name)}};
    }
    return out;
  }

  /// A warning where easi's grouping changes what a chain of Impala means, once per chain and run.
  void warnGrouping(const ImpalaChain& chain) const {
    static std::mutex mutex;
    static std::set<std::string> warned;
    const std::string at = file() + "@" + std::to_string(chain.line);
    const std::lock_guard<std::mutex> lock(mutex);
    if (!warned.insert(at + "/" + chain.name + "/" + chain.op).second) {
      return;
    }
    logWarning() << "easi:" << at << ": the Impala of" << chain.name << "has a chain where"
                 << chain.op << "is followed by" << chain.following
                 << "; easi's conversion into Lua groups it to the right (a - b + c as a - (b + "
                    "c)), and the walker computes what easi computes.";
  }

  /// The file and the variable of an !ASAGI, and the axes ASAGI reads there, as many as `in` has
  /// dimensions.
  std::vector<datafield::AsagiAxis> asagiAxes(const YAML::Node& node, const Values& in) {
    const std::string path = dataPath(as<std::string>(node, "file"));
    const std::string variable = node["var"] ? as<std::string>(node, "var") : std::string("data");
    const auto axes = datafield::asagiAxesOf(path, variable);
    if (!axes.has_value()) {
      fail(node, "Could not open " + path + " with ASAGI.");
    }
    if (axes->size() != in.size()) {
      fail(node,
           "ASAGI requires " + std::to_string(axes->size()) + "D input (got " + listOf(keysOf(in)) +
               ").");
    }
    for (const auto& axis : *axes) {
      if (axis.nodes > 1 && !(axis.first < axis.last)) {
        throw EasiUnsupported(where(node) + ": a grid of ASAGI whose coordinates decrease");
      }
    }
    return *axes;
  }

  /// A grid of ASAGI, sampled at the dimensions of `in` in their order, as easi samples it.
  Values grid(const YAML::Node& node, const Values& in) {
    static_cast<void>(asagiAxes(node, in));
    const auto parameters = as<std::vector<std::string>>(node, "parameters");
    const std::string interpolation = gridInterpolation(node);
    datafield::GridSpec spec;
    spec.kind = "asagi";
    spec.path = dataPath(as<std::string>(node, "file"));
    spec.variable = node["var"] ? as<std::string>(node, "var") : std::string("data");
    spec.components = parameters;
    spec.interpolation = {interpolation};
    spec.boundary = {"clamp"};
    datafield::GridDesc desc;
    try {
      desc = datafield::gridDescOf(spec);
    } catch (const std::invalid_argument& error) {
      fail(node, error.what());
    }
    const auto id = program_.internGrid(desc);
    program_.nameGrid(id, expr::GridNames{"asagi" + std::to_string(id), parameters});
    std::vector<NodeId> coordinates;
    for (const auto& [dimension, value] : in) {
      coordinates.push_back(value.node);
    }
    Values out;
    for (std::size_t p = 0; p < parameters.size(); ++p) {
      out[parameters[p]] =
          Value{arena().lookup(id, static_cast<std::int32_t>(p), coordinates), {positionOf(node)}};
    }
    return out;
  }

  /// The interpolation of a grid of easi (`interpolation`, by default linear): nearest, linear or
  /// cubic, as easi 1.7 reads it.
  std::string gridInterpolation(const YAML::Node& node) const {
    if (!node["interpolation"]) {
      return "linear";
    }
    const auto given = as<std::string>(node, "interpolation");
    std::string interpolation = given;
    std::transform(interpolation.begin(),
                   interpolation.end(),
                   interpolation.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (interpolation != "nearest" && interpolation != "linear" && interpolation != "cubic") {
      fail(node,
           "Invalid interpolation type " + given + " (expected one of nearest, linear, cubic).");
    }
    return interpolation;
  }

  /// The parameters of an !SCECFile, by easi's names, and the columns of the file they are.
  static const std::map<std::string, std::int32_t>& scecParameters() {
    static const std::map<std::string, std::int32_t> Parameters = {{"cohesion", 8},
                                                                   {"d_c", 7},
                                                                   {"forced_rupture_time", 9},
                                                                   {"mu_d", 6},
                                                                   {"mu_s", 5},
                                                                   {"s_dip", 2},
                                                                   {"s_normal", 0},
                                                                   {"s_strike", 1}};
    return Parameters;
  }

  /// An !SCECFile: the grid of its file (ScecGrid.h), sampled at the two dimensions of `in` in
  /// their order, as easi samples it -- the first along strike. Its parameters are eight of the
  /// ten values of a node, by easi's names; the file is read when the grid is loaded, its header
  /// here.
  Values scecGrid(const YAML::Node& node, const Values& in) {
    if (in.size() != 2) {
      fail(node, "SCECFile requires 2D input (got " + listOf(keysOf(in)) + ").");
    }
    const std::string path = dataPath(as<std::string>(node, "file"));
    try {
      static_cast<void>(datafield::readScecHeader(path));
    } catch (const std::runtime_error& error) {
      fail(node, error.what());
    }
    datafield::GridSpec spec;
    spec.kind = "scec";
    spec.path = path;
    spec.components = datafield::scecComponentNames();
    spec.interpolation = {gridInterpolation(node)};
    spec.boundary = {"clamp"};
    datafield::GridDesc desc;
    try {
      desc = datafield::gridDescOf(spec);
    } catch (const std::invalid_argument& error) {
      fail(node, error.what());
    }
    const auto id = program_.internGrid(desc);
    program_.nameGrid(id, expr::GridNames{"scec" + std::to_string(id), spec.components});
    std::vector<NodeId> coordinates;
    for (const auto& [dimension, value] : in) {
      coordinates.push_back(value.node);
    }
    Values out;
    for (const auto& [name, column] : scecParameters()) {
      out[name] = Value{arena().lookup(id, column, coordinates), {positionOf(node)}};
    }
    return out;
  }

  /// A layered model: the nodes of its parameters along the one dimension its map gives, a grid
  /// of them (Grid.h) -- lower and upper the steps previous and next, linear linear.
  Values layers(const YAML::Node& node,
                const Values& in,
                const std::set<std::string>& wanted,
                NodeId reach) {
    if (!node["map"]) {
      fail(node, "the component has no `map`");
    }
    const YAML::Node map = node["map"];
    std::string interpolation = as<std::string>(node, "interpolation");
    std::transform(interpolation.begin(),
                   interpolation.end(),
                   interpolation.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const auto nodes = as<std::map<double, std::vector<double>>>(node, "nodes");
    const auto parameters = as<std::vector<std::string>>(node, "parameters");
    if (interpolation != "lower" && interpolation != "upper" && interpolation != "linear") {
      fail(node, "Invalid interpolation type " + interpolation + ".");
    }
    if (nodes.empty()) {
      fail(node, "a layered model has nodes");
    }
    std::vector<double> knots;
    std::vector<double> data;
    bool infinite = false;
    for (const auto& [depth, values] : nodes) {
      if (values.size() != parameters.size()) {
        fail(node, "Number of parameters must match number of node entries.");
      }
      if (std::isnan(depth)) {
        fail(node, "Not a number nodes are forbidden.");
      }
      if (std::isinf(depth)) {
        if (interpolation == "linear") {
          fail(node, "Infinite nodes are forbidden for interpolation = linear.");
        }
        infinite = true;
      }
      knots.push_back(depth);
      data.insert(data.end(), values.begin(), values.end());
    }

    // the one dimension of its map: its own, or what its components supply
    Value along;
    const auto mapChildren = childrenOf(map);
    if (mapChildren.empty()) {
      const Values out = mapped(map, in, reach);
      if (out.size() != 1) {
        fail(node, "Layered model requires 1D input.");
      }
      along = out.begin()->second;
    } else {
      const auto supplied = suppliedBy(map, keysOf(in));
      if (supplied.size() != 1) {
        fail(node, "Layered model requires 1D input.");
      }
      along = evaluate(map, in, supplied, reach).begin()->second;
    }
    // easi's layers accept every depth but NaN
    const YAML::Node at = map;
    require(reach,
            pw(Fn::Eq, along.node, along.node),
            "Could not find model for point (" + where(at) + ").");

    if (infinite) {
      return layerSteps(node, knots, data, parameters, interpolation, along, wanted, reach);
    }
    if (interpolation == "linear") {
      return layerLines(node, knots, data, parameters, along, wanted, reach);
    }
    datafield::GridSpec spec;
    spec.kind = "inline";
    spec.components = parameters;
    spec.knots = {knots};
    spec.values = data;
    spec.interpolation = {interpolation == "lower" ? "previous" : "next"};
    spec.boundary = {"clamp"};
    datafield::GridDesc desc;
    try {
      desc = datafield::gridDescOf(spec);
    } catch (const std::invalid_argument& error) {
      fail(node, error.what());
    }
    const auto id = program_.internGrid(desc);
    program_.nameGrid(id, expr::GridNames{"layers" + std::to_string(id), parameters});
    Values out;
    std::set<std::string> missing;
    for (const auto& parameter : wanted) {
      const auto found = std::find(parameters.begin(), parameters.end(), parameter);
      if (found == parameters.end()) {
        missing.insert(parameter);
        continue;
      }
      out[parameter] = Value{
          arena().lookup(id, static_cast<std::int32_t>(found - parameters.begin()), {along.node}),
          {positionOf(node)}};
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Component does not supply all required parameters. Missing: {" +
                  listOf(missing) + "}.");
    }
    return out;
  }

  /// A layered model that interpolates linearly, as easi builds it: a layer between two nodes is a
  /// polynomial of the slope and the intercept easi computes from them, which it evaluates as
  /// easi evaluates a polynomial, the intercept plus the slope times the depth; up to the first
  /// node a depth has the values of the first, beyond the last those of the last, and on a node
  /// those of the layer below it. The slopes and the intercepts are a grid of the nodes, each with
  /// those of the layer up to it, whose step `next` finds the layer of a depth.
  Values layerLines(const YAML::Node& node,
                    const std::vector<double>& knots,
                    const std::vector<double>& data,
                    const std::vector<std::string>& parameters,
                    const Value& along,
                    const std::set<std::string>& wanted,
                    NodeId reach) {
    const std::size_t count = parameters.size();
    const std::size_t nodes = knots.size();
    // per node and parameter the slope and the intercept of the layer up to it, as easi's
    // LayeredModelBuilder computes them; the first node has those of the first layer, which no
    // depth reads there
    std::vector<double> coefficients(2 * count * nodes, 0.0);
    for (std::size_t k = 1; k < nodes; ++k) {
      const double lower = knots[k - 1];
      const double upper = knots[k];
      for (std::size_t p = 0; p < count; ++p) {
        const double below = data[(k - 1) * count + p];
        const double above = data[k * count + p];
        coefficients[2 * (k * count + p)] = (below - above) / (lower - upper);
        coefficients[2 * (k * count + p) + 1] = (-upper * below + lower * above) / (lower - upper);
      }
    }
    if (nodes > 1) {
      std::copy_n(coefficients.begin() + 2 * count, 2 * count, coefficients.begin());
    }
    std::vector<std::string> names;
    for (const auto& parameter : parameters) {
      names.push_back(parameter + "_slope");
      names.push_back(parameter + "_intercept");
    }
    expr::GridId id = expr::NoGrid;
    if (nodes > 1) {
      datafield::GridSpec spec;
      spec.kind = "inline";
      spec.components = names;
      spec.knots = {knots};
      spec.values = coefficients;
      spec.interpolation = {"next"};
      spec.boundary = {"clamp"};
      datafield::GridDesc desc;
      try {
        desc = datafield::gridDescOf(spec);
      } catch (const std::invalid_argument& error) {
        fail(node, error.what());
      }
      id = program_.internGrid(desc);
      program_.nameGrid(id, expr::GridNames{"layers" + std::to_string(id), names});
    }
    const NodeId d = along.node;
    Values out;
    std::set<std::string> missing;
    for (const auto& parameter : wanted) {
      const auto found = std::find(parameters.begin(), parameters.end(), parameter);
      if (found == parameters.end()) {
        missing.insert(parameter);
        continue;
      }
      const auto p = static_cast<std::size_t>(found - parameters.begin());
      const NodeId first = konst(data[p]);
      if (nodes == 1) {
        out[parameter] = Value{first, {positionOf(node)}};
        continue;
      }
      const NodeId slope = arena().lookup(id, static_cast<std::int32_t>(2 * p), {d});
      const NodeId intercept = arena().lookup(id, static_cast<std::int32_t>(2 * p + 1), {d});
      const NodeId line = pw(Fn::Add, intercept, pw(Fn::Mul, slope, d));
      const NodeId last = konst(data[(nodes - 1) * count + p]);
      out[parameter] = Value{pw(Fn::Select,
                                pw(Fn::Le, d, konst(knots.front())),
                                first,
                                pw(Fn::Select, pw(Fn::Lt, konst(knots.back()), d), last, line)),
                             {positionOf(node)}};
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Component does not supply all required parameters. Missing: {" +
                  listOf(missing) + "}.");
    }
    return out;
  }

  /// The steps of a layered model with a node at an infinite depth, which a grid has not: easi's
  /// layers one after another, each from a node to the next (the first from minus infinity, the
  /// last to infinity), a depth on a node in the first that contains it, as a select between them.
  Values layerSteps(const YAML::Node& node,
                    const std::vector<double>& knots,
                    const std::vector<double>& data,
                    const std::vector<std::string>& parameters,
                    const std::string& interpolation,
                    const Value& along,
                    const std::set<std::string>& wanted,
                    NodeId reach) {
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    const std::size_t count = parameters.size();
    // per layer its limits and the node whose values it takes
    struct Layer {
      double lower;
      double upper;
      std::size_t node;
    };
    std::vector<Layer> layers;
    for (std::size_t n = 0; n <= knots.size(); ++n) {
      const double lower = n == 0 ? -Infinity : knots[n - 1];
      const double upper = n == knots.size() ? Infinity : knots[n];
      const std::size_t taken = n == 0                     ? 0
                                : n == knots.size()        ? knots.size() - 1
                                : interpolation == "lower" ? n - 1
                                                           : n;
      layers.push_back(Layer{lower, upper, taken});
    }
    Values out;
    std::set<std::string> missing;
    const NodeId d = along.node;
    for (const auto& parameter : wanted) {
      const auto found = std::find(parameters.begin(), parameters.end(), parameter);
      if (found == parameters.end()) {
        missing.insert(parameter);
        continue;
      }
      const auto p = static_cast<std::size_t>(found - parameters.begin());
      NodeId value = konst(std::numeric_limits<double>::quiet_NaN());
      for (std::size_t l = layers.size(); l-- > 0;) {
        const NodeId within = pw(
            Fn::And, pw(Fn::Le, konst(layers[l].lower), d), pw(Fn::Le, d, konst(layers[l].upper)));
        value = pw(Fn::Select, within, konst(data[layers[l].node * count + p]), value);
      }
      out[parameter] = Value{value, {positionOf(node)}};
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Component does not supply all required parameters. Missing: {" +
                  listOf(missing) + "}.");
    }
    return out;
  }

  /// Whether a component, as a component of a composite, accepts a point at `in`: NoNode where it
  /// accepts every point.
  NodeId accepts(const YAML::Node& node, const Values& in) {
    const std::string tag = node.Tag();
    if (tag == "!Include") {
      const auto [root, path] = included(node);
      files_.push_back(path);
      const NodeId accepted = accepts(root, in);
      files_.pop_back();
      return accepted;
    }
    if (tag == "!GroupFilter") {
      const YAML::Node groups = node["groups"];
      if (!groups) {
        fail(node, "the component has no `groups`");
      }
      std::vector<int> listed;
      try {
        if (groups.IsScalar()) {
          listed.push_back(groups.as<int>());
        } else {
          listed = groups.as<std::vector<int>>();
        }
      } catch (const YAML::Exception& error) {
        fail(groups, std::string("`groups` is not of the form easi reads: ") + error.what());
      }
      NodeId accepted = konst(0.0);
      for (const int group : std::set<int>(listed.begin(), listed.end())) {
        accepted = pw(Fn::Or, accepted, pw(Fn::Eq, groupNode(), konst(group)));
      }
      return accepted;
    }
    if (tag == "!AxisAlignedCuboidalDomainFilter") {
      const auto limits = as<std::map<std::string, std::pair<double, double>>>(node, "limits");
      matchDimensions(node, limits, in);
      NodeId accepted = NoNode;
      for (const auto& [dimension, range] : limits) {
        const NodeId x = in.at(dimension).node;
        accepted = both(
            accepted,
            pw(Fn::And, pw(Fn::Le, konst(range.first), x), pw(Fn::Le, x, konst(range.second))));
      }
      return accepted;
    }
    if (tag == "!SphericalDomainFilter") {
      const auto centre = as<std::map<std::string, double>>(node, "centre");
      const auto radius = as<double>(node, "radius");
      matchDimensions(node, centre, in);
      NodeId distance = konst(0.0);
      for (const auto& [dimension, coordinate] : centre) {
        const NodeId difference = pw(Fn::Sub, in.at(dimension).node, konst(coordinate));
        distance = pw(Fn::Add, distance, pw(Fn::Mul, difference, difference));
      }
      return pw(Fn::Le, distance, konst(radius * radius));
    }
    if (tag == "!ASAGI") {
      // within the grid, as ASAGI bounds it: from the first coordinate by its scaling, a dimension
      // of one node unbounded (Level::open())
      const auto axes = asagiAxes(node, in);
      NodeId accepted = NoNode;
      std::size_t d = 0;
      for (const auto& [dimension, value] : in) {
        const auto& axis = axes[d++];
        double lo = -std::numeric_limits<double>::infinity();
        double hi = std::numeric_limits<double>::infinity();
        if (axis.nodes > 1) {
          const double scaling = (axis.last - axis.first) / static_cast<double>(axis.nodes - 1);
          const double far = axis.first + scaling * static_cast<double>(axis.nodes - 1);
          lo = std::min(axis.first, far);
          hi = std::max(axis.first, far);
        }
        accepted =
            both(accepted,
                 pw(Fn::And, pw(Fn::Le, konst(lo), value.node), pw(Fn::Le, value.node, konst(hi))));
      }
      return accepted;
    }
    // the other maps, !Switch, !Any and !IdentityMap
    return NoNode;
  }

  /// The group of the point: the input `group`, or 0 for points without one.
  NodeId groupNode() {
    if (!query_.group) {
      return konst(0.0);
    }
    if (!groupDeclared_) {
      program_.addInput("group", DataType::F64);
      groupDeclared_ = true;
    }
    return arena().field("group");
  }

  /// That the dimensions a filter names are those of its input.
  template <typename Map>
  void matchDimensions(const YAML::Node& node, const Map& named, const Values& in) const {
    std::set<std::string> dimensions;
    for (const auto& [dimension, value] : named) {
      dimensions.insert(dimension);
    }
    if (dimensions != keysOf(in)) {
      fail(node,
           "The input parameters of a component ({" + listOf(dimensions) +
               "}) added to a composite are incompatible (should be {" + listOf(keysOf(in)) +
               "}).");
    }
  }

  /// The parameters `wanted` that the component gives at `in`, for the points that reach it
  /// (`reach`, NoNode for every point).
  Values evaluate(const YAML::Node& node,
                  const Values& in,
                  const std::set<std::string>& wanted,
                  NodeId reach) {
    const std::string tag = node.Tag();
    if (tag == "!Include") {
      // easi reads an included file with the dimensions of the root
      const auto [root, path] = included(node);
      if (keysOf(in) != rootDimensions_) {
        fail(node,
             "The input parameters of a component ({" + listOf(rootDimensions_) +
                 "}) added to a composite are incompatible (should be {" + listOf(keysOf(in)) +
                 "}).");
      }
      files_.push_back(path);
      Values values = evaluate(root, in, wanted, reach);
      files_.pop_back();
      return values;
    }
    if (tag == "!Switch") {
      return switched(node, in, wanted, reach);
    }
    if (tag == "!LayeredModel") {
      return layers(node, in, wanted, reach);
    }
    if (tag == "!Any" || tag == "!IdentityMap" || tag == "!GroupFilter" ||
        tag == "!AxisAlignedCuboidalDomainFilter" || tag == "!SphericalDomainFilter") {
      return composite(node, in, wanted, reach);
    }
    if (tag == "!ConstantMap" || tag == "!AffineMap" || tag == "!PolynomialMap" ||
        tag == "!EvalModel" || tag == "!ASAGI" || tag == "!SCECFile" || tag == "!FunctionMap" ||
        tag == "!LuaMap" || easiSpecialOf(tag) != nullptr) {
      return composite(node, mapped(node, in, reach), wanted, reach);
    }
    fail(node, "Unknown tag " + tag);
  }

  /// A !Switch: each parameter from the component that names it.
  Values switched(const YAML::Node& node,
                  const Values& in,
                  const std::set<std::string>& wanted,
                  NodeId reach) {
    std::set<std::string> all;
    std::vector<std::pair<std::set<std::string>, YAML::Node>> branches;
    for (const auto& entry : node) {
      std::set<std::string> restrictions;
      try {
        if (entry.first.IsSequence()) {
          const auto names = entry.first.as<std::vector<std::string>>();
          restrictions.insert(names.begin(), names.end());
        } else {
          restrictions.insert(entry.first.as<std::string>());
        }
      } catch (const YAML::Exception& error) {
        fail(node, std::string("a branch of a switch names parameters: ") + error.what());
      }
      std::set<std::string> overlap;
      std::set_intersection(restrictions.begin(),
                            restrictions.end(),
                            all.begin(),
                            all.end(),
                            std::inserter(overlap, overlap.end()));
      if (!overlap.empty()) {
        fail(node,
             "Switch requires a partition of parameters. The parameters {" + listOf(overlap) +
                 "} were specified in multiple branches.");
      }
      all.insert(restrictions.begin(), restrictions.end());
      branches.emplace_back(std::move(restrictions), entry.second);
    }
    Values out;
    std::set<std::string> missing;
    for (const auto& parameter : wanted) {
      if (all.count(parameter) == 0) {
        missing.insert(parameter);
      }
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Switch is not complete with respect to request. Missing: {" +
                  listOf(missing) + "}.");
    }
    for (const auto& [restrictions, child] : branches) {
      std::set<std::string> sub;
      std::set_intersection(restrictions.begin(),
                            restrictions.end(),
                            wanted.begin(),
                            wanted.end(),
                            std::inserter(sub, sub.end()));
      if (sub.empty()) {
        continue;
      }
      for (auto& [name, value] : evaluate(child, in, sub, reach)) {
        out[name] = std::move(value);
      }
    }
    return out;
  }

  /// A composite of the values `values` (its input mapped): its own where it has no components,
  /// else those of the first of its components that accepts the point.
  Values composite(const YAML::Node& node,
                   const Values& values,
                   const std::set<std::string>& wanted,
                   NodeId reach) {
    const auto children = childrenOf(node);
    if (children.empty()) {
      Values out;
      std::set<std::string> missing;
      for (const auto& parameter : wanted) {
        const auto found = values.find(parameter);
        if (found == values.end()) {
          missing.insert(parameter);
        } else {
          out[parameter] = found->second;
        }
      }
      if (!missing.empty()) {
        require(reach,
                konst(0.0),
                where(node) + ": Component does not supply all required parameters. Missing: {" +
                    listOf(missing) + "}.");
      }
      return out;
    }
    const std::string message = "Could not find model for point (" + where(node) + ").";
    // the first component that accepts a point takes it; the later ones see the rest
    std::vector<NodeId> accepted;
    std::vector<Values> results;
    NodeId before = NoNode; // that no component before accepts
    NodeId any = konst(0.0);
    bool always = false;
    for (const auto& child : children) {
      const NodeId accepts = this->accepts(child, values);
      const NodeId reaches = both(reach, both(before, accepts));
      accepted.push_back(accepts);
      results.push_back(evaluate(child, values, wanted, reaches));
      if (accepts == NoNode) {
        always = true;
        break;
      }
      any = pw(Fn::Or, any, accepts);
      before = both(before, negation(accepts));
    }
    if (!always) {
      require(reach, any, message);
    }
    // from the last component back: a select of each that accepts, before the rest, a value it
    // does not give NaN (where a check fails)
    const NodeId nan = konst(std::numeric_limits<double>::quiet_NaN());
    Values out;
    for (const auto& parameter : wanted) {
      NodeId node = nan;
      std::vector<expr::SourcePosition> from;
      bool given = false;
      for (std::size_t c = results.size(); c-- > 0;) {
        const auto found = results[c].find(parameter);
        const bool gives = found != results[c].end() && found->second.node != NoNode;
        const NodeId here = gives ? found->second.node : nan;
        if (gives) {
          given = true;
          for (const auto& position : found->second.from) {
            if (std::find(from.begin(), from.end(), position) == from.end()) {
              from.push_back(position);
            }
          }
        }
        node = accepted[c] == NoNode ? here : pw(Fn::Select, accepted[c], here, node);
      }
      if (given) {
        out[parameter] = Value{node, std::move(from)};
      }
    }
    return out;
  }

  /// The parameters a component supplies at `in` (easi's suppliedParameters()): those of its
  /// components, or its own.
  std::set<std::string> suppliedBy(const YAML::Node& node, const std::set<std::string>& in) {
    const std::string tag = node.Tag();
    if (tag == "!Include") {
      const auto [root, path] = included(node);
      files_.push_back(path);
      auto supplied = suppliedBy(root, rootDimensions_);
      files_.pop_back();
      return supplied;
    }
    if (tag == "!Switch") {
      std::set<std::string> supplied;
      for (const auto& entry : node) {
        const auto sub = suppliedBy(entry.second, in);
        supplied.insert(sub.begin(), sub.end());
      }
      return supplied;
    }
    if (tag == "!LayeredModel") {
      const auto parameters = as<std::vector<std::string>>(node, "parameters");
      return {parameters.begin(), parameters.end()};
    }
    // the outputs of the component itself
    std::set<std::string> own;
    if (tag == "!Any" || tag == "!IdentityMap" || tag == "!GroupFilter") {
      own = in;
    } else if (tag == "!AxisAlignedCuboidalDomainFilter" || tag == "!SphericalDomainFilter") {
      own = in;
    } else if (tag == "!ConstantMap" || tag == "!PolynomialMap" || tag == "!FunctionMap") {
      for (const auto& [name, value] : as<std::map<std::string, YAML::Node>>(node, "map")) {
        own.insert(name);
      }
    } else if (tag == "!AffineMap") {
      for (const auto& [name, value] : as<std::map<std::string, YAML::Node>>(node, "matrix")) {
        own.insert(name);
      }
    } else if (tag == "!EvalModel" || tag == "!ASAGI") {
      const auto names = as<std::vector<std::string>>(node, "parameters");
      own.insert(names.begin(), names.end());
    } else if (tag == "!LuaMap") {
      const YAML::Node returns = node["returns"];
      if (returns && returns.IsScalar()) {
        own.insert(returns.as<std::string>());
      } else if (returns) {
        const auto names = returns.as<std::vector<std::string>>();
        own.insert(names.begin(), names.end());
      }
    } else if (const auto* special = easiSpecialOf(tag); special != nullptr) {
      own.insert(special->out.begin(), special->out.end());
    } else if (tag == "!SCECFile") {
      for (const auto& [name, column] : scecParameters()) {
        own.insert(name);
      }
    } else {
      fail(node, "Unknown tag " + tag);
    }
    const auto children = childrenOf(node);
    if (children.empty()) {
      return own;
    }
    std::set<std::string> supplied;
    for (const auto& child : children) {
      const auto sub = suppliedBy(child, own);
      supplied.insert(sub.begin(), sub.end());
    }
    return supplied;
  }

  EasiQuery query_;
  std::set<std::string> rootDimensions_;
  expr::Program program_;
  std::vector<std::string> files_;
  std::vector<std::pair<NodeId, std::string>> checks_;
  /// the text of the files read, by path: where a scalar is a block scalar
  std::map<std::string, std::string> texts_;
  bool groupDeclared_{false};
};

} // namespace

expr::Program easiProgram(const std::string& path, const EasiQuery& query) {
  return Walker(query).translate(path);
}

std::vector<std::string> easiParameters(const std::string& path,
                                        const std::vector<std::string>& dimensions) {
  EasiQuery query;
  query.dimensions = dimensions;
  const auto parameters = Walker(query).parameters(path);
  return {parameters.begin(), parameters.end()};
}

EasiProgramReader::EasiProgramReader(std::string path,
                                     const std::vector<std::string>& dimensions,
                                     ReferenceMaker makeReference)
    : path_(std::move(path)), makeReference_(makeReference) {
  // easi orders the dimensions by name
  const std::set<std::string> sorted(dimensions.begin(), dimensions.end());
  dimensions_.assign(sorted.begin(), sorted.end());
  parameters_ = easiParameters(path_, dimensions_);
}

EasiProgramReader::~EasiProgramReader() = default;

void EasiProgramReader::prepare(const DataTable& table) {
  std::vector<std::string> wanted;
  bool group = false;
  for (const auto& entry : table.dataEntries()) {
    if (entry.direction != Direction::In) {
      wanted.push_back(entry.name);
    } else if (entry.name == "group") {
      group = true;
    }
  }
  if (compiled_ == nullptr || wanted != madeFor_ || group != madeWithGroup_) {
    EasiQuery query;
    query.dimensions = dimensions_;
    query.parameters = wanted;
    query.group = group;
    compiled_ = std::make_unique<CompiledReader>(
        easiProgram(path_, query),
        makeReference_ != nullptr ? makeReference_(path_, dimensions_) : nullptr);
    madeFor_ = std::move(wanted);
    madeWithGroup_ = group;
  }
  compiled_->prepare(table);
}

void EasiProgramReader::call(const DataTable& table) {
  if (compiled_ == nullptr) {
    prepare(table);
  }
  compiled_->call(table);
}

} // namespace seissol::reader::scripting
