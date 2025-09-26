// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITLAYOUT_H_

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initLayout(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITLAYOUT_H_
