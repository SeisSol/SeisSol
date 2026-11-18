// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The Python module of sderiv (`import sderiv`): loads a program -- sderiv in the full language or
// the core format, a Lua model, or the program the easi walker makes of an easi file --, validates
// it, prints it in the normal form of the core format, and evaluates it with the interpreter, the
// reference of every backend, at points given as NumPy arrays. What the command line tool does with
// a program, for scripts and tests: the tests of a frontend, say, which compare the program it
// writes with the model it was written from.

#include "BuildInfo.h"
#include "Expr/Backend.h"
#include "Expr/Core.h"
#include "Expr/Program.h"
#include "Expr/Record.h"
#include "Expr/SderivFrontend.h"
#include "Model/ScriptTypes.h"
#include "Parallel/MPI.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/EasiWalker.h"
#include "Reader/Scripting/LuaFrontend.h"
#include "Reader/Scripting/LuaSyntax.h"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iterator>
#include <memory>
#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

using seissol::expr::Program;

/// An error of a source: what sderiv.Error carries, its file, line and column where known (0
/// where not).
class SourceError : public std::runtime_error {
  public:
  SourceError(const std::string& what, std::string file, int line, int column)
      : std::runtime_error(what), file_(std::move(file)), line_(line), column_(column) {}
  [[nodiscard]] const std::string& file() const { return file_; }
  [[nodiscard]] int line() const { return line_; }
  [[nodiscard]] int column() const { return column_; }

  private:
  std::string file_;
  int line_;
  int column_;
};

/// The grids of the programs, and MPI, which loading a grid needs (AsagiLite reads collectively):
/// started on the first evaluation of a program with grids, unless the process started it (as
/// mpi4py does), and ended when Python ends.
struct Session {
  bool mpiReady{false};
  bool mpiStarted{false};
  std::unique_ptr<seissol::reader::datafield::GridStore> grids;
};

Session& session() {
  static Session state;
  return state;
}

seissol::reader::datafield::GridStore& gridsFor(const Program& program) {
  auto& state = session();
  if (!program.grids().empty() && !state.mpiReady) {
#ifdef USE_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized != 0) {
      seissol::Mpi::mpi.setComm(MPI_COMM_WORLD);
    } else {
      int argc = 0;
      char** argv = nullptr;
      seissol::Mpi::mpi.init(argc, argv);
      state.mpiStarted = true;
    }
#else
    int argc = 0;
    char** argv = nullptr;
    seissol::Mpi::mpi.init(argc, argv);
#endif
    state.mpiReady = true;
  }
  if (state.grids == nullptr) {
    state.grids = std::make_unique<seissol::reader::datafield::GridStore>();
  }
  return *state.grids;
}

void endSession() {
  auto& state = session();
  state.grids.reset();
  if (state.mpiStarted) {
    seissol::Mpi::finalize();
    state.mpiStarted = false;
  }
}

std::string readFile(const std::string& path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    throw SourceError("cannot read " + path, path, 0, 0);
  }
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

bool startsWith(const std::string& text, const std::string& prefix) {
  return text.compare(0, prefix.size(), prefix) == 0;
}

bool endsWith(const std::string& text, const std::string& suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/// The program of the sderiv text `text` of the file `file`: in the core format where it starts
/// with the header, else in the full language.
Program compileText(const std::string& text, const std::string& file) {
  seissol::expr::SderivOptions options;
  options.file = file;
  try {
    return seissol::expr::compileSderivModule(text, options);
  } catch (const seissol::expr::CoreError& error) {
    throw SourceError(error.what(), error.where().file, error.where().line, error.where().column);
  } catch (const seissol::expr::SderivError& error) {
    if (error.position() < 0) {
      throw SourceError((file.empty() ? "" : file + ": ") + error.what(), file, 0, 0);
    }
    const auto at =
        seissol::expr::positionOf(text, static_cast<std::size_t>(error.position()), file);
    throw SourceError((file.empty() ? "" : file + ":") + std::to_string(at.line) + ":" +
                          std::to_string(at.column) + ": " + error.what(),
                      file,
                      at.line,
                      at.column);
  } catch (const std::invalid_argument& error) {
    throw SourceError((file.empty() ? "" : file + ": ") + error.what(), file, 0, 0);
  }
}

/// The program of the Lua model of the file `path` (the Lua frontend).
Program translateLua(const std::string& path) {
  namespace lua = seissol::reader::scripting::lua;
  const lua::LuaOrigin origin{path, 1, ':'};
  try {
    return lua::translateLuaModule(readFile(path), origin);
  } catch (const lua::LuaError& error) {
    const int line = error.line();
    throw SourceError(
        (line > 0 ? origin.at(line) + ": " : path + ": ") + error.message(), path, line, 0);
  }
}

/// The program of the file `path`: sderiv, else a Lua model (`.lua`, `lua:`).
Program load(const std::string& path) {
  if (endsWith(path, ".lua") || startsWith(path, "lua:")) {
    return translateLua(startsWith(path, "lua:") ? path.substr(4) : path);
  }
  if (endsWith(path, ".yaml") || endsWith(path, ".yml") || startsWith(path, "easi:")) {
    throw SourceError(path + ": an easi file is read with load_easi(), which is given the "
                             "parameters to compute",
                      path,
                      0,
                      0);
  }
  const std::string file = startsWith(path, "sderiv:") ? path.substr(7) : path;
  return compileText(readFile(file), file);
}

/// The program the easi walker makes of the easi file `path`, for the parameters `parameters` (all
/// it supplies where none are given) at points of the dimensions `dimensions`, with or without the
/// input `group`.
Program loadEasi(const std::string& path,
                 const std::optional<std::vector<std::string>>& parameters,
                 const std::vector<std::string>& dimensions,
                 bool group) {
  namespace sc = seissol::reader::scripting;
  try {
    sc::EasiQuery query;
    query.dimensions = dimensions;
    query.parameters = parameters.has_value() ? *parameters : sc::easiParameters(path, dimensions);
    query.group = group;
    return sc::easiProgram(path, query);
  } catch (const sc::EasiUnsupported& error) {
    throw SourceError(error.what(), path, 0, 0);
  } catch (const sc::EasiModelError& error) {
    throw SourceError(error.what(), path, 0, 0);
  }
}

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// The values of `names` in `given`, a mapping of names to numbers or arrays of one dimension. A
/// name `given` lacks takes its value of `defaults`, where there are those; else it is an error
/// that names it.
std::vector<Array> arraysOf(const std::vector<std::string>& names,
                            const py::dict& given,
                            const std::vector<double>* defaults,
                            const char* what) {
  std::vector<Array> arrays;
  arrays.reserve(names.size());
  for (std::size_t i = 0; i < names.size(); ++i) {
    const py::str name(names[i]);
    if (!given.contains(name)) {
      if (defaults == nullptr) {
        throw py::key_error(std::string("the program reads the ") + what + " `" + names[i] +
                            "`, which is not given");
      }
      arrays.emplace_back(std::vector<py::ssize_t>{}, &(*defaults)[i]);
      continue;
    }
    auto array = Array::ensure(given[name]);
    if (!array || array.ndim() > 1) {
      throw py::value_error(std::string("the ") + what + " `" + names[i] +
                            "` is a number or an array of one dimension");
    }
    arrays.push_back(std::move(array));
  }
  return arrays;
}

/// The number of points of `arrays` of `names` and the points `count` found so far: that of their
/// arrays of more than one value, which agree; a number, or an array of one value, is at every
/// point.
void countPoints(const std::vector<Array>& arrays,
                 const std::vector<std::string>& names,
                 std::optional<std::size_t>& count) {
  for (std::size_t i = 0; i < arrays.size(); ++i) {
    const auto size = static_cast<std::size_t>(arrays[i].size());
    if (arrays[i].ndim() == 0 || size == 1) {
      continue;
    }
    if (count.has_value() && *count != size) {
      throw py::value_error(
          "the values have different numbers of points: " + std::to_string(*count) + " and " +
          std::to_string(size) + " (`" + names[i] + "`)");
    }
    count = size;
  }
}

/// The values of `arrays` at `points` points, [array][point].
std::vector<double> columnsOf(const std::vector<Array>& arrays, std::size_t points) {
  std::vector<double> columns;
  columns.reserve(arrays.size() * points);
  for (const auto& array : arrays) {
    const double* data = array.data();
    for (std::size_t k = 0; k < points; ++k) {
      columns.push_back(array.size() == 1 ? data[0] : data[k]);
    }
  }
  return columns;
}

/// A mapping of names to values, as a dict: none an empty one.
py::dict dictOf(const py::object& mapping) {
  return mapping.is_none() ? py::dict() : py::dict(mapping);
}

/// The inputs of a call: those of the mapping `inputs`, and those given by name.
py::dict merged(const py::object& inputs, const py::kwargs& named) {
  py::dict all = dictOf(inputs);
  for (const auto& [key, value] : named) {
    all[key] = value;
  }
  return all;
}

/// The values of `program` at `points` points, [output, then next state][point], by the
/// interpreter. The GIL is kept: the grids of the module are shared.
std::vector<double> interpret(const Program& program,
                              const std::vector<double>& inputs,
                              std::size_t points,
                              const std::vector<double>* states) {
  return seissol::expr::interpretValues(
      program, gridsFor(program), inputs, points, nullptr, states, true);
}

/// `count` values from `values`, from `first` on, as an array.
py::array_t<double>
    arrayOf(const std::vector<double>& values, std::size_t first, std::size_t count) {
  py::array_t<double> array(static_cast<py::ssize_t>(count));
  std::copy_n(values.begin() + static_cast<std::ptrdiff_t>(first), count, array.mutable_data());
  return array;
}

std::vector<std::string> namesOf(const std::vector<seissol::expr::VarSpec>& specs) {
  std::vector<std::string> names;
  names.reserve(specs.size());
  for (const auto& spec : specs) {
    names.push_back(spec.name);
  }
  return names;
}

/// The outputs of `values` of `program` at `points` points, by name.
py::dict outputsOf(const Program& program, const std::vector<double>& values, std::size_t points) {
  py::dict outputs;
  for (std::size_t o = 0; o < program.outputs().size(); ++o) {
    outputs[py::str(program.outputs()[o].name)] = arrayOf(values, o * points, points);
  }
  return outputs;
}

py::dict evaluate(const Program& program, const py::object& inputs, const py::kwargs& named) {
  if (!program.state().empty()) {
    throw py::value_error("the program has state, which a call advances: evaluate it with step()");
  }
  const auto names = namesOf(program.inputs());
  const auto arrays = arraysOf(names, merged(inputs, named), nullptr, "input");
  std::optional<std::size_t> count;
  countPoints(arrays, names, count);
  const std::size_t points = count.value_or(1);
  return outputsOf(program, interpret(program, columnsOf(arrays, points), points, nullptr), points);
}

py::tuple step(const Program& program, const py::object& inputs, const py::object& states) {
  const auto inputNames = namesOf(program.inputs());
  const auto inputArrays = arraysOf(inputNames, dictOf(inputs), nullptr, "input");
  std::vector<std::string> stateNames;
  std::vector<double> initial;
  for (const auto& state : program.state()) {
    stateNames.push_back(state.name);
    initial.push_back(state.initial);
  }
  const auto stateArrays = arraysOf(stateNames, dictOf(states), &initial, "state");
  // as many points as the arrays of the inputs and of the states have
  std::optional<std::size_t> count;
  countPoints(inputArrays, inputNames, count);
  countPoints(stateArrays, stateNames, count);
  const std::size_t points = count.value_or(1);
  const auto stateColumns = columnsOf(stateArrays, points);
  const auto values = interpret(program, columnsOf(inputArrays, points), points, &stateColumns);
  py::dict next;
  const std::size_t base = program.outputs().size();
  for (std::size_t s = 0; s < program.state().size(); ++s) {
    next[py::str(program.state()[s].name)] = arrayOf(values, (base + s) * points, points);
  }
  return py::make_tuple(outputsOf(program, values, points), next);
}

std::string computeTypeOf(const Program& program) {
  return program.computeType() == seissol::expr::ComputeType::F32 ? "f32" : "f64";
}

std::string describe(const Program& program) {
  const auto list = [](const std::vector<std::string>& names) {
    std::string text;
    for (const auto& name : names) {
      text += (text.empty() ? "" : ", ") + name;
    }
    return text;
  };
  return "<sderiv.Program " + seissol::expr::fingerprintText(program) + ": " +
         list(namesOf(program.inputs())) + " -> " + list(namesOf(program.outputs())) + ">";
}

} // namespace

PYBIND11_MODULE(sderiv, module) {
  module.doc() = "sderiv, the programs SeisSol evaluates: load a program (sderiv in the full "
                 "language or the core format, a Lua model, an easi file through the easi "
                 "walker), validate it, print it in the normal form of the core format, and "
                 "evaluate it with the interpreter at points given as NumPy arrays.";
  module.attr("__version__") = seissol::BuildInfo::VersionString;

  // the materials of SeisSol as types of the full language, as a run has them
  seissol::model::registerScriptTypes();

  // sderiv.Error, a ValueError with the file, the line and the column of the source where known;
  // sderiv.EvaluationError, a RuntimeError where a point does not meet a check of the program
  PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> errorType;
  errorType.call_once_and_store_result([&]() -> py::object {
    return py::exception<SourceError>(module, "Error", PyExc_ValueError);
  });
  py::register_exception<seissol::expr::EvaluationError>(
      module, "EvaluationError", PyExc_RuntimeError);
  py::register_exception_translator([](std::exception_ptr pointer) {
    try {
      if (pointer) {
        std::rethrow_exception(pointer);
      }
    } catch (const SourceError& error) {
      const py::object& type = errorType.get_stored();
      py::object value = type(error.what());
      value.attr("file") = error.file().empty() ? py::object(py::none()) : py::str(error.file());
      value.attr("line") = error.line() > 0 ? py::object(py::int_(error.line())) : py::none();
      value.attr("column") = error.column() > 0 ? py::object(py::int_(error.column())) : py::none();
      py::set_error(type, value);
    }
  });

  py::class_<Program>(module,
                      "Program",
                      "A program: what SeisSol evaluates of a script, the same whatever "
                      "frontend wrote it.")
      .def_property_readonly(
          "fingerprint",
          [](const Program& program) { return seissol::expr::fingerprintText(program); },
          "The fingerprint of what the program computes, 16 hexadecimal digits: the same for "
          "two programs that compute the same, whatever their names.")
      .def_property_readonly(
          "inputs",
          [](const Program& program) { return namesOf(program.inputs()); },
          "The names of the inputs, in their order.")
      .def_property_readonly(
          "outputs",
          [](const Program& program) { return namesOf(program.outputs()); },
          "The names of the outputs, in their order.")
      .def_property_readonly(
          "states",
          [](const Program& program) {
            py::dict states;
            for (const auto& state : program.state()) {
              states[py::str(state.name)] = state.initial;
            }
            return states;
          },
          "The states, by name, with the values they start from.")
      .def_property_readonly(
          "compute_type", &computeTypeOf, "The type the program computes in: 'f64' or 'f32'.")
      .def_property_readonly(
          "meta",
          [](const Program& program) {
            py::list entries;
            for (const auto& entry : program.meta()) {
              py::list values;
              for (const auto& value : entry.values) {
                if (value.quoted) {
                  values.append(py::str(value.text));
                } else {
                  values.append(py::int_(py::str(value.text)));
                }
              }
              entries.append(py::make_tuple(entry.key, values));
            }
            return entries;
          },
          "The meta lines of the program (its provenance), each a key and its values, strings "
          "and whole numbers.")
      .def_property_readonly(
          "grids",
          [](const Program& program) {
            std::vector<std::string> labels;
            for (const auto& grid : program.grids()) {
              labels.push_back(grid.label());
            }
            return labels;
          },
          "The grids the program samples.")
      .def_property_readonly(
          "checks",
          [](const Program& program) {
            std::vector<std::string> messages;
            for (const auto& check : program.checks()) {
              messages.push_back(check.message);
            }
            return messages;
          },
          "The messages of the checks every point has to meet.")
      .def(
          "validate",
          [](const Program& program) { seissol::expr::validate(program); },
          "Raises ValueError where the program is not valid. One that load() or compile() gave "
          "is.")
      .def(
          "core",
          [](const Program& program) { return seissol::expr::printCore(program); },
          "The program in the normal form of the core format; ValueError for one the format "
          "cannot express.")
      .def("evaluate",
           &evaluate,
           py::arg("inputs") = py::none(),
           "The outputs at points, by name, as arrays: evaluate({'x': x, 'y': y}) or "
           "evaluate(x=x, y=y). The inputs are numbers or arrays of one dimension, a number (or "
           "an array of one value) at every point. Computed by the interpreter, in the type of "
           "the program, the bits of every backend of SeisSol on the CPU. Raises "
           "EvaluationError where a point does not meet a check. For a program without state; "
           "one with state is evaluated with step().")
      .def("step",
           &step,
           py::arg("inputs") = py::none(),
           py::arg("states") = py::none(),
           "One call of a program with state: (outputs, next states), both by name, from the "
           "inputs and the states of the points (by default the values they start from).")
      .def("__str__",
           [](const Program& program) {
             try {
               return seissol::expr::printCore(program);
             } catch (const std::exception&) {
               return describe(program);
             }
           })
      .def("__repr__", &describe);

  module.def("load",
             &load,
             py::arg("path"),
             "The program of a file: sderiv, in the full language or the core format (a file that "
             "starts with #!sderiv), or a Lua model (.lua or lua:). Raises sderiv.Error, with the "
             "file, the line and the column.");
  module.def("compile",
             &compileText,
             py::arg("text"),
             py::arg("file") = "",
             "The program of sderiv text: the core format where it starts with #!sderiv, else "
             "the full language. `file` names it in errors, and its relative grid paths are "
             "relative to its directory.");
  module.def("load_easi",
             &loadEasi,
             py::arg("path"),
             py::arg("parameters") = py::none(),
             py::arg("dimensions") = std::vector<std::string>{"x", "y", "z"},
             py::arg("group") = true,
             "The program the easi walker makes of an easi file, for the parameters asked for "
             "(by default all it supplies) at points of `dimensions`, with the input `group` "
             "unless group=False.");

  // the grids and MPI, when Python ends
  py::module_::import("atexit").attr("register")(py::cpp_function(&endSession));
}
