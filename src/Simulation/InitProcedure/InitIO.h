// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_INITPROCEDURE_INITIO_H_
#define SEISSOL_SRC_SIMULATION_INITPROCEDURE_INITIO_H_

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
/**
 * Creates the output directory, with its parents, and sets the output prefix inside it: the
 * directory of the outputs (output.outdir), or with output.outsimdir, a directory of the run in
 * there, named by output.outname and the ID of the run.
 *
 * Done before anything else of the run, so that what the run records while it initializes has a
 * place to go. Collective; rank 0 creates it.
 */
void createOutputDirectory(seissol::SeisSol& seissolInstance);

void initIO(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif // SEISSOL_SRC_SIMULATION_INITPROCEDURE_INITIO_H_
