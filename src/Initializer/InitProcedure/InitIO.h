// SPDX-FileCopyrightText: 2023 SeisSol Group
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
/**
 * Creates the directory the output prefix lies in, with its parents, if it does not exist yet.
 *
 * Done before anything else of the run, so that what the run records while it initializes has a
 * place to go. Collective; rank 0 creates it.
 */
void createOutputDirectory(const seissol::SeisSol& seissolInstance);

void initIO(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITIO_H_
