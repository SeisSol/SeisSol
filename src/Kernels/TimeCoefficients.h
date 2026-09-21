// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_TIMECOEFFICIENTS_H_
#define SEISSOL_SRC_KERNELS_TIMECOEFFICIENTS_H_

#include "Common/Constants.h"
#include "Kernels/Precision.h"

#include <array>
#include <vector>

namespace seissol::kernels {

/**
 * The coefficients of one time operation -- an integration over an interval, an
 * evaluation at a point -- for every expansion a solver keeps.
 *
 * A solver that transports only its state keeps one expansion, and the two sets
 * are the same. A solver that transports more keeps the second in another
 * basis: the coefficients of the state fall out of a recursion, the rest have
 * to be won from samples, and the basis that is well conditioned for winning
 * them is not the one the recursion produces.
 *
 * Which set belongs to which basis is therefore a property of the solver, and
 * the two travel together because the choice between coefficient sets -- the
 * whole step or a subinterval of it -- is one choice over a pair.
 */
struct TimeCoefficients {
  std::array<real, ConvergenceOrder> state{};
  std::array<real, ConvergenceOrder> extra{};
};

/**
 * A quadrature over one step in time: where to sample it, with what weight, and
 * the coefficients that read every expansion at each node.
 *
 * A solver whose flux is nonlinear in the state cannot integrate that flux in
 * closed form and samples it instead, so the three belong together -- and
 * pairing a rule with the coefficients of both expansions is the job the time
 * basis already does for a single operation. It is therefore filled there, and
 * the order of these members is what it fills.
 */
struct TimeQuadrature {
  std::vector<double> nodes;
  std::vector<double> weights;
  std::vector<TimeCoefficients> coefficients;
};

/**
 * What a predictor needs of one step in time, beyond how wide the step is.
 *
 * Every predictor needs the coefficients that integrate over the whole step.
 * One whose flux it cannot integrate in closed form also needs the rule it
 * samples that flux with -- and that rule is empty for a solver which does not
 * sample, which is what `Solver::RequiresTimeQuadrature` says.
 *
 * Both are questions about the step and not about a cell, so they are answered
 * once for the step and handed down. A per-cell call that formed a time basis
 * of its own would answer the same question once per cell.
 */
struct TimeStepCoefficients {
  TimeCoefficients integral;
  TimeQuadrature quadrature;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIMECOEFFICIENTS_H_
