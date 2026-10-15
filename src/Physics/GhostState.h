// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_GHOSTSTATE_H_
#define SEISSOL_SRC_PHYSICS_GHOSTSTATE_H_

#include "Initializer/Typedefs.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace seissol::expr {
class Program;
} // namespace seissol::expr

namespace seissol::physics {

class InitialField;

/// The ghost state of a boundary given by a script, integrated over a time step: one compiled
/// program that evaluates the script at the points of the Gauss-Legendre rule of the step and adds
/// up the weighted values, as kernels::AnalyticalBoundary does with a function of the time. Only
/// what depends on the time (and on the state inside) is computed per time point.
///
/// The script gives components of the state by the names `components`; what else it gives is not
/// asked for. It reads x, y, z, t, sim, rho, mu and lambda, and the state inside by the names
/// `inner`, which then are inputs per time point. A component it does not give is zero
/// (Missing::Zero, the analytical boundary) or the state inside (Missing::Inner, the nonlinear
/// Dirichlet boundary, whose inner names are its components).
class GhostStateIntegral {
  public:
  enum class Missing : std::uint8_t { Zero, Inner };

  /// Throws std::invalid_argument for a script that reads anything else or keeps state; `name`
  /// names it in the message.
  GhostStateIntegral(const expr::Program& ghost,
                     const std::string& name,
                     std::vector<std::string> components,
                     std::vector<std::string> inner,
                     Missing missing,
                     std::size_t timePoints);
  ~GhostStateIntegral();

  GhostStateIntegral(const GhostStateIntegral&) = delete;
  GhostStateIntegral& operator=(const GhostStateIntegral&) = delete;
  GhostStateIntegral(GhostStateIntegral&&) = delete;
  GhostStateIntegral& operator=(GhostStateIntegral&&) = delete;

  [[nodiscard]] std::size_t timePointCount() const;

  /// The components of the ghost state, the names `components`.
  [[nodiscard]] std::size_t componentCount() const;

  /// The components of the state inside that the integral reads per time point (the names
  /// `inner`, or none if the script does not read the state inside and leaves nothing to it).
  [[nodiscard]] std::size_t innerCount() const;

  /// The times of the quadrature of the step [start, start + width], one per time point.
  void timePoints(double start, double width, double* times) const;

  /// The integral of the ghost state over the step [start, start + width] at `count` points, in
  /// the simulation `simulation`, beyond a cell of the material `materialData`: component j at
  /// point p into ghost[j * count + p]. The state inside, if read, is inner component k at point p
  /// at time point i in inner[(i * innerCount() + k) * count + p].
  void evaluate(double start,
                double width,
                std::size_t simulation,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                const double* inner,
                double* ghost) const;

  private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/// The ghost state of the analytical boundary (FaceType::Analytical) from an initial condition
/// that is a script that compiles, integrated over a step: per fused simulation, the integral of
/// the script of the field that serves it. Fields of one script share an integral, which then
/// reads the simulation as an input, so that one program serves all simulations.
class AnalyticalGhostState {
  public:
  /// The ghost state of the fields `fields`, of which field s % fields.size() serves simulation s
  /// of `simulations`, as for the initial condition; `components` are the components the boundary
  /// takes. Null unless all fields are scripts that compile (physics::ScriptField): the boundary
  /// then evaluates the fields at one time after the other. A script that cannot serve as a field
  /// is an error.
  static std::unique_ptr<AnalyticalGhostState>
      fromFields(const std::vector<std::unique_ptr<InitialField>>& fields,
                 const std::vector<std::string>& components,
                 std::size_t simulations,
                 std::size_t timePoints);

  /// The integral of simulation `simulation`.
  [[nodiscard]] const GhostStateIntegral& integral(std::size_t simulation) const;

  /// What the script of simulation `simulation` reads as `sim`: the simulation of its field.
  [[nodiscard]] std::size_t scriptSimulation(std::size_t simulation) const;

  private:
  std::vector<std::unique_ptr<GhostStateIntegral>> integrals_;
  // per simulation, its integral and what its script reads as `sim`
  std::vector<std::size_t> integralOf_;
  std::vector<std::size_t> scriptSimulation_;
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_GHOSTSTATE_H_
