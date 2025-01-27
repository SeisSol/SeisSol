// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITIO_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITIO_H_

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initIO(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITIO_H_
