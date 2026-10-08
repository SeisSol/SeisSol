// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCRIPTFIELD_H_
#define SEISSOL_SRC_PHYSICS_SCRIPTFIELD_H_

#include "Physics/InitialField.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace seissol::physics {

/// A field given by a script: an initial condition, an analytical solution, and with it the values
/// of an analytic boundary condition, which ask for it at their points and times as they go.
///
/// The script gives the quantities by name. It reads `x`, `y`, `z`, and as it needs them the time
/// `t`, the simulation `sim`, and the material of the cell `rho`, `mu` and `lambda`. A quantity it
/// does not give is zero.
///
/// A script that compiles to a program (an sderiv module, or a Lua model that traces) is evaluated
/// by one kernel per thread, made here, with the points and values of a call passed as bases --
/// the boundary condition asks for a few points per face and time step, which is far too often to
/// bind a table each time. Any other script (easi) is evaluated by one reader per thread, which is
/// correct but slow.
class ScriptField : public InitialFieldOf<ScriptField> {
  public:
  /// A script that reads anything else, keeps state, or gives what is no quantity is an error
  /// once the field is evaluated -- not before, since a script that serves as an initial condition
  /// only is fine with that (it may read the mesh group, for one).
  ///
  /// `hasTime` says whether an easi script reads the time, as for the initial condition.
  ScriptField(const std::string& path,
              std::vector<std::string> quantities,
              std::size_t simulation,
              bool hasTime);
  ~ScriptField() override;

  ScriptField(const ScriptField&) = delete;
  ScriptField& operator=(const ScriptField&) = delete;
  ScriptField(ScriptField&&) = delete;
  ScriptField& operator=(ScriptField&&) = delete;

  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const {
    auto& values = scratch(count);
    evaluateValues(time, points, count, materialData, values.data());
    for (std::size_t j = 0; j < quantityCount(); ++j) {
      for (std::size_t i = 0; i < count; ++i) {
        dofsQP(i, j) = static_cast<RealT>(values[j * count + i]);
      }
    }
  }

  /// Quantity j at point i, into values[j * count + i].
  void evaluateValues(double time,
                      const std::array<double, 3>* points,
                      std::size_t count,
                      const CellMaterialData& materialData,
                      double* values) const;

  [[nodiscard]] std::size_t quantityCount() const;

  /// Whether the script compiled, and calls run a kernel rather than a reader.
  [[nodiscard]] bool compiled() const;

  private:
  [[nodiscard]] std::vector<double>& scratch(std::size_t count) const;

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_SCRIPTFIELD_H_
