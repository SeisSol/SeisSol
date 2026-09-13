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

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIMECOEFFICIENTS_H_
