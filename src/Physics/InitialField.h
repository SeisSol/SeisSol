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

#include <array>
#include <cstddef>

namespace seissol::physics {

/// A field that is known analytically. It is evaluated in the reals of whoever asks: into views of
/// float and of double alike.
class InitialField {
  public:
  virtual ~InitialField() = default;
  virtual void evaluate(double time,
                        const std::array<double, 3>* points,
                        std::size_t count,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2, float, unsigned>& dofsQP) const = 0;
  virtual void evaluate(double time,
                        const std::array<double, 3>* points,
                        std::size_t count,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2, double, unsigned>& dofsQP) const = 0;
};

/// Implements `evaluate` of `Base`, an `InitialField`, for both types of reals by `evaluateIn` of
/// `Derived`, a template of the type of the reals.
template <typename Derived, typename Base = InitialField>
class InitialFieldOf : public Base {
  public:
  using Base::Base;

  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, float, unsigned>& dofsQP) const override {
    static_cast<const Derived*>(this)->evaluateIn(time, points, count, materialData, dofsQP);
  }
  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, double, unsigned>& dofsQP) const override {
    static_cast<const Derived*>(this)->evaluateIn(time, points, count, materialData, dofsQP);
  }
};

class ZeroField : public InitialFieldOf<ZeroField> {
  public:
  template <typename RealT>
  void evaluateIn(double /*time*/,
                  const std::array<double, 3>* /*points*/,
                  std::size_t /*count*/,
                  const CellMaterialData& /*materialData*/,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const {
    dofsQP.setZero();
  }
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
