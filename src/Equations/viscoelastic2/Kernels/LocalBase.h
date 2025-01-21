// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_LOCALBASE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_LOCALBASE_H_

#include "Physics/InitialField.h"
#include "generated_code/kernel.h"
#include <memory>

namespace seissol {
namespace kernels {
class LocalBase {
  protected:
  double gravitationalAcceleration;
  kernel::volumeExt m_volumeKernelPrototype;
  kernel::localFluxExt m_localFluxKernelPrototype;
  kernel::local m_localKernelPrototype;
  const std::vector<std::unique_ptr<physics::InitialField>>* initConds;

  public:
  virtual void setInitConds(decltype(initConds) initConds) { this->initConds = initConds; }

  void setGravitationalAcceleration(double g) { gravitationalAcceleration = g; }
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_LOCALBASE_H_
