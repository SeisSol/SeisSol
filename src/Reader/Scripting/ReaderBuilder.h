// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_READERBUILDER_H_
#define SEISSOL_SRC_READER_SCRIPTING_READERBUILDER_H_

#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Scripting/DataReader.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

/// The reader of the script `path`: `easi:` or no prefix for an easi file, `lua:` for a Lua model
/// (compiled if it traces), `sderiv:` for an sderiv module.
std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs);

/// A reader of the script `path` that evaluates it point by point rather than through a compiled
/// program, and hence can be made once per thread and called with a new table every time: an
/// easi file or a Lua model; null for an sderiv module, which has no such reader.
std::unique_ptr<DataReader> buildInterpretedReader(const std::string& path,
                                                   const std::vector<std::string>& defaultInArgs);

/// The program of the script `path`, if it compiles to one: an sderiv module (`sderiv:` or a
/// .sderiv file), or a Lua model that traces (`lua:` or a .lua file). Nothing for an easi file, or
/// for a Lua model that does not trace, in which case `reason` (if given) says why. An sderiv
/// module that does not compile is a configuration error. `sderivOptions` are passed on to the
/// sderiv compiler.
std::optional<expr::Program> buildProgram(const std::string& path,
                                          std::string* reason = nullptr,
                                          const expr::SderivOptions& sderivOptions = {});

} // namespace seissol::reader::scripting
#endif // SEISSOL_SRC_READER_SCRIPTING_READERBUILDER_H_
