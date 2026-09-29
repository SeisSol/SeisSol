// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
#define SEISSOL_SRC_PHYSICS_INITIALFIELD_H_

#include "GeneratedCode/init.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

#include <array>
#include <cstddef>

namespace seissol::physics {

class InitialField {
  public:
  virtual ~InitialField() = default;
  virtual void evaluate(double time,
                        const std::array<double, 3>* points,
                        std::size_t count,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2, real, unsigned>& dofsQP) const = 0;
};

class ZeroField : public InitialField {
  public:
  void evaluate(double /*time*/,
                const std::array<double, 3>* /*points*/,
                std::size_t /*count*/,
                const CellMaterialData& /*materialData*/,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override {
    dofsQP.setZero();
  }
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
