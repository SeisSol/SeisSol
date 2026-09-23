// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_IMPEDANCE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_IMPEDANCE_H_

#include "Equations/ImpedanceBase.h"
#include "Equations/elastic/Model/Impedance.h"
#include "Equations/viscoelastic/Model/Datastructures.h"

#include <cstddef>

namespace seissol::model {

/**
 * The anelastic functions do not enter the flux: the Riemann problem at a face only sees the
 * elastic part of the Jacobian with the unrelaxed moduli (see getTransposedGodunovState, which
 * forwards to the elastic material). The admittance is therefore the elastic one.
 *
 * Shared by both viscoelastic solvers, like the material itself.
 */
template <std::size_t Mechanisms>
struct ImpedanceCompute<ViscoElasticMaterial<Mechanisms>>
    : public ImpedanceCompute<ElasticMaterial> {};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_IMPEDANCE_H_
