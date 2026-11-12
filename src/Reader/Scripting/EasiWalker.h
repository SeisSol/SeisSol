// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_EASIWALKER_H_
#define SEISSOL_SRC_READER_SCRIPTING_EASIWALKER_H_

// The easi walker: an easi file (YAML) translated into a program of src/Expr, which the compiled
// backends evaluate, rather than evaluated by easi. It reads the file as easi 1.x does -- the tags,
// the includes and the paths of data files resolved as easi resolves them -- and computes what
// easi computes, in the order easi computes it, so that the bits are easi's: an affine map as
// easi's matrix product (from zero, in the order of the dimensions, the translation last), a
// polynomial by Horner's scheme as easi has it, a grid of ASAGI or of a SCEC file (ScecGrid.h)
// through the grid store (whose interpolation is that of easi 1.7; an earlier easi interpolates
// within a few units in the last place of it), a special component by easi's formulas
// (EasiSpecial.h). A layered model becomes a grid of its nodes (Grid.h), whose step finds the layer
// of a depth: the value of a node, or the slope and the intercept easi computes for the layer.
//
// easi hands a point to the first component of a composite that accepts it, and fails where none
// does; a component in a list is translated as a select between the components, and the points
// no component accepts fail a check of the program, with easi's message -- guarded by what the
// points that reach the component are, so that a component no point reaches fails nothing, as in
// easi. A parameter that a component does not supply fails the same way, for the points that
// reach it.
//
// The Lua of an !LuaMap, and the Lua easi makes of the Impala of a !FunctionMap (ImpalaToLua.h),
// the Lua frontend translates (LuaFrontend.h), as easi runs them. Every component of easi 1.7 is
// translated; what the walker does not translate -- Lua beyond the subset of the frontend, a grid
// of ASAGI whose coordinates decrease -- it refuses with EasiUnsupported, naming the component and
// the line; the caller evaluates the file with easi then (buildReader()).

#include "Expr/Program.h"
#include "Reader/Scripting/DataReader.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

class CompiledReader;

/// A component of an easi file that the walker does not translate, with where it is.
class EasiUnsupported : public std::runtime_error {
  public:
  using std::runtime_error::runtime_error;
};

/// An error of an easi file that easi reports when it reads the file: an unknown tag, an
/// attribute missing or of another form, dimensions of a component that do not match.
class EasiModelError : public std::runtime_error {
  public:
  using std::runtime_error::runtime_error;
};

/// What a consumer asks of an easi model.
struct EasiQuery {
  /// the dimensions of the points (x, y and z, say), as EasiReader is given them; easi orders
  /// them by name
  std::vector<std::string> dimensions;
  /// the parameters asked for, the outputs of the program in this order
  std::vector<std::string> parameters;
  /// whether the points have a group (the input `group`); without, a point is of group 0, as in
  /// easi
  bool group{true};
};

/// The program of the easi model of the file `path` for `query`. Throws EasiUnsupported for a
/// component it does not translate, EasiModelError for a file that easi would not read.
[[nodiscard]] expr::Program easiProgram(const std::string& path, const EasiQuery& query);

/// The parameters the model of `path` supplies, by name (easi's suppliedParameters()), for
/// points of the dimensions `dimensions`. Throws as easiProgram().
[[nodiscard]] std::vector<std::string> easiParameters(const std::string& path,
                                                      const std::vector<std::string>& dimensions);

/// An easi file evaluated through the walker: for the parameters the table asks for (and whether
/// it has a group), the program of easiProgram(), which a CompiledReader evaluates -- with
/// `reference`, made by `makeReference` for the same table, as the oracle it is compared with at
/// its first call and the reader it falls back to, where given.
class EasiProgramReader : public DataReader {
  public:
  using ReferenceMaker = std::unique_ptr<DataReader> (*)(const std::string&,
                                                         const std::vector<std::string>&);

  /// `dimensions` as EasiReader takes them; `makeReference` null for no oracle.
  EasiProgramReader(std::string path,
                    const std::vector<std::string>& dimensions,
                    ReferenceMaker makeReference);
  ~EasiProgramReader() override;

  EasiProgramReader(const EasiProgramReader&) = delete;
  EasiProgramReader& operator=(const EasiProgramReader&) = delete;
  EasiProgramReader(EasiProgramReader&&) = delete;
  EasiProgramReader& operator=(EasiProgramReader&&) = delete;

  const std::vector<std::string>& inputVars() override { return dimensions_; }
  const std::vector<std::string>& outputVars() override { return parameters_; }

  void prepare(const DataTable& table) override;
  void call(const DataTable& table) override;

  /// The reader of the parameters the table prepared last asks for.
  [[nodiscard]] const CompiledReader* compiled() const { return compiled_.get(); }

  private:
  std::string path_;
  std::vector<std::string> dimensions_;
  std::vector<std::string> parameters_;
  ReferenceMaker makeReference_;
  /// what the compiled reader was made for: the parameters asked for, and whether with a group
  std::vector<std::string> madeFor_;
  bool madeWithGroup_{false};
  std::unique_ptr<CompiledReader> compiled_;
};

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_EASIWALKER_H_
