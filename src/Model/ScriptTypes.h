// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_SCRIPTTYPES_H_
#define SEISSOL_SRC_MODEL_SCRIPTTYPES_H_

#include <string>
#include <vector>

namespace seissol::model {

/// The materials of SeisSol as types of the full language of sderiv: `elastic`, `acoustic`,
/// `anisotropic`, `poroelastic`, `viscoelastic`, `viscoacoustic` and `plasticity`, each the record
/// of the parameters a material file gives (MaterialT::ParameterList), in their order, named as
/// MaterialT::Text names it (without the number of mechanisms). Registered once with
/// expr::registerRecordType(), wherever scripts are read; later calls do nothing.
void registerScriptTypes();

/// The name of the type of a material of the text `text`: "viscoelastic" of "viscoelastic-3".
[[nodiscard]] std::string scriptTypeName(const std::string& text);

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_SCRIPTTYPES_H_
