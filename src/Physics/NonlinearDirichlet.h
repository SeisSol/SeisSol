// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_NONLINEARDIRICHLET_H_
#define SEISSOL_SRC_PHYSICS_NONLINEARDIRICHLET_H_

#include "Physics/GhostState.h"

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
/// names of the quantities of the material (`s_xx`, ..., `v1`, `v2`, `v3` for an elastic one), x,
/// y, z, t, sim, and the parameters of the material (model::scriptParameterNames); it gives the
/// ghost state by the same names. A quantity it does not give is that of the inside. In an sderiv
/// module, a quantity always reads the inner state, also next to the definition of its ghost
/// value: `out def v1 = -v1` mirrors the normal velocity, in the face-aligned basis.
///
/// The boundary evaluates the condition at the nodes of a face, at the times of the quadrature of
/// a time step, with the inner state of the cell at that time, and integrates it over the step:
/// all in one program (integral(), physics::GhostStateIntegral). It is not folded into the flux
/// solver, as the affine condition of a Dirichlet boundary is.
///
/// A script that gives `frame` as 1 states both in the face-aligned basis, as with the Dirichlet
/// boundary: the first axis is the outward normal, the others are the strike and the dip
/// direction. 0, the default, is global coordinates. The frame is the same on all faces, so it
/// may read nothing.
class NonlinearDirichlet {
  public:
  /// The condition of the script at `path`, for a material with the quantities `quantities` and the
  /// parameters `parameters`, integrated over a step with `timePoints` points. A script that cannot
  /// be used is a configuration error.
  NonlinearDirichlet(const std::string& path,
                     std::vector<std::string> quantities,
                     std::vector<std::string> parameters,
                     std::size_t timePoints);

  /// The condition of a program, which `name` names in messages; from sderiv, it is compiled with
  /// the quantities as its inputs (expr::SderivOptions::inputs). Throws std::invalid_argument for a
  /// program that cannot be used.
  NonlinearDirichlet(const expr::Program& program,
                     const std::string& name,
                     std::vector<std::string> quantities,
                     std::vector<std::string> parameters,
                     std::size_t timePoints);

  ~NonlinearDirichlet();

  NonlinearDirichlet(const NonlinearDirichlet&) = delete;
  NonlinearDirichlet& operator=(const NonlinearDirichlet&) = delete;
  NonlinearDirichlet(NonlinearDirichlet&&) = delete;
  NonlinearDirichlet& operator=(NonlinearDirichlet&&) = delete;

  /// Whether the condition is stated in the face-aligned basis (and not in global coordinates).
  [[nodiscard]] bool faceAligned() const;

  [[nodiscard]] std::size_t quantityCount() const;

  /// The ghost state integrated over a step, from the state inside at its time points (all
  /// quantities, in the basis of the frame), into the basis of the frame.
  [[nodiscard]] const GhostStateIntegral& integral() const;

  /// The same for a device, in a device build, which takes a quantity the script does not give
  /// from the integrated state of the cell, and evaluates the rest: the integral of the quantities
  /// the script gives, which reads the state inside only where the script does. Null in a build
  /// without a device.
  [[nodiscard]] const GhostStateIntegral* deviceIntegral() const;

  private:
  void setUp(const expr::Program& program, std::size_t timePoints);

  std::string name_;
  std::vector<std::string> quantities_;
  std::vector<std::string> parameters_;
  bool faceAligned_{false};
  std::unique_ptr<GhostStateIntegral> integral_;
  std::unique_ptr<GhostStateIntegral> deviceIntegral_;
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_NONLINEARDIRICHLET_H_
