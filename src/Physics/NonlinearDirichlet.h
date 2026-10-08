// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_NONLINEARDIRICHLET_H_
#define SEISSOL_SRC_PHYSICS_NONLINEARDIRICHLET_H_

#include "Initializer/Typedefs.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace seissol::expr {
class Program;
} // namespace seissol::expr

namespace seissol::physics {

/// The condition of a nonlinear Dirichlet boundary (FaceType::NonlinearDirichlet): the state
/// beyond a face, the ghost state, as a function of the state inside.
///
/// A script gives it, an sderiv module or a Lua model that traces. It reads the inner state by the
/// names of the quantities of the material (`s_xx`, ..., `v1`, `v2`, `v3` for an elastic one), and
/// x, y, z, t, sim, rho, mu and lambda; it gives the ghost state by the same names. A quantity it
/// does not give is that of the inside. In an sderiv module, a quantity always reads the inner
/// state, also next to the definition of its ghost value: `out def v1 = -v1` mirrors the normal
/// velocity, in the face-aligned basis.
///
/// The boundary evaluates the condition at the nodes of a face, at the times of the quadrature of
/// a time step, with the inner state of the cell at that time (kernels::ApplyNonlinearDirichlet);
/// it is not folded into the flux solver, as the affine condition of a Dirichlet boundary is.
///
/// A script that gives `frame` as 1 states both in the face-aligned basis, as with the Dirichlet
/// boundary: the first axis is the outward normal, the others are the strike and the dip
/// direction. 0, the default, is global coordinates. The frame is the same on all faces, so it
/// may read nothing.
class NonlinearDirichlet {
  public:
  /// The condition of the script at `path`, for a material with the quantities `quantities`. A
  /// script that cannot be used is a configuration error.
  NonlinearDirichlet(const std::string& path, std::vector<std::string> quantities);

  /// The condition of a program, which `name` names in messages; from sderiv, it is compiled with
  /// the quantities as its inputs (expr::SderivOptions::inputs). Throws std::invalid_argument for a
  /// program that cannot be used.
  NonlinearDirichlet(const expr::Program& program,
                     const std::string& name,
                     std::vector<std::string> quantities);

  ~NonlinearDirichlet();

  NonlinearDirichlet(const NonlinearDirichlet&) = delete;
  NonlinearDirichlet& operator=(const NonlinearDirichlet&) = delete;
  NonlinearDirichlet(NonlinearDirichlet&&) = delete;
  NonlinearDirichlet& operator=(NonlinearDirichlet&&) = delete;

  /// Whether the condition is stated in the face-aligned basis (and not in global coordinates).
  [[nodiscard]] bool faceAligned() const;

  [[nodiscard]] std::size_t quantityCount() const;

  /// The ghost state at `count` points at `time`, in the simulation `simulation`, beyond a cell of
  /// the material `materialData`: from the inner state, quantity j at point i in
  /// inner[j * count + i], into ghost likewise, both in the basis of the frame.
  void evaluate(double time,
                std::size_t simulation,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                const double* inner,
                double* ghost) const;

  private:
  void setUp(const expr::Program& program);

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_NONLINEARDIRICHLET_H_
