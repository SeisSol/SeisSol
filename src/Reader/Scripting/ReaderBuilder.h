// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_READERBUILDER_H_
#define SEISSOL_SRC_READER_SCRIPTING_READERBUILDER_H_

#include "Reader/Scripting/DataReader.h"

#include <memory>
namespace seissol::reader::scripting {

std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs);

} // namespace seissol::reader::scripting
#endif // SEISSOL_SRC_READER_SCRIPTING_READERBUILDER_H_
