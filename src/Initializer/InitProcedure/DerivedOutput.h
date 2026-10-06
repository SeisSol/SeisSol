// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_DERIVEDOUTPUT_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_DERIVEDOUTPUT_H_

// Derived quantities at the output points of the wave field, computed by an expression program.
//
// The program is written pointwise. It reads the quantities by name, and the consumer decides
// where each of them comes from; here, from the coefficients of the cell:
//
//   <q>                 quantity q of the solution, e.g. v1 or s_xx
//   <q>_r0, _r1, _r2    its derivative along the reference coordinates xi, eta, zeta
//   dx_<q>, dy_<q>, dz_<q>
//                       its derivative in space, by the chain rule over the inverse Jacobian
//   int_<q>, int_<q>_rk, d*_int_<q>
//                       the same for the time integral of q
//   <p>, <p>_rk, d*_<p> the plastic strain quantity p, e.g. ep_xx (with plasticity)
//   jinv<k><d>          d xi_k / d x_d, the inverse Jacobian of the cell
//   x, y, z             the output point
//   t, dt               the time, and the time since the previous evaluation of the point
//
// Every value and reference derivative becomes a contraction of the coefficients of its quantity
// against a projection matrix (expr::Kind::Contract), so the program never sees a basis function
// and the coefficients are not projected into a buffer first. Interning makes a contraction one
// node however many outputs read it: the strain and rotation outputs together need the nine
// reference derivatives of the velocities once, not once per output.
//
// The subcells of a refined output are stacked into the rows of one matrix: the points of a cell
// are its subcells times the output points of a subcell, so one program and one call cover them
// all, and a state slot (a maximum over time, say) exists once per output point.
//
// CADENCE. A program without state is evaluated when the output is written, for all written
// points at once, into a buffer the writer copies from. A program with state accumulates: the
// points of a layer are evaluated after every time step of its cluster, and the buffer then
// always holds the running values. A layer whose cluster did not step since the last write (at
// the start, or one that runs on a device) is evaluated when the output is written.

#include "Expr/Program.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "Numerical/Projection.h"
#include "Reader/Scripting/DataTable.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer {

/// A quantity a derived program can read, stored as coefficients per cell.
struct DerivedSource {
  std::string name;
  /// Point values at the nodal set of the volume (plastic strain) rather than modal coefficients.
  bool nodal{false};
};

/// Where the output points lie: the subcells of a cell and the points of a subcell, both in
/// reference coordinates, and how the solution is transferred onto them.
struct DerivedGeometry {
  std::vector<io::instance::geometry::Subcell<3>> subcells{io::instance::geometry::unrefined<3>()};
  std::vector<std::array<double, 3>> dataBase;
  /// The polynomial degree of the output points; only relevant for Target::Project.
  std::size_t dataOrder{0};
  /// The convergence order of the solution.
  std::size_t order{1};
  numerical::projection::Target target{numerical::projection::Target::Interpolate};
  numerical::projection::NodalSet nodalSet{numerical::projection::NodalSet::WarpBlend};
};

/// A derived-output program prepared for the output points of a geometry: its channels resolved
/// against the vocabulary above, the contractions in place, and the projection matrices built.
class DerivedProgram {
  public:
  /// Throws std::invalid_argument for a channel outside the vocabulary, and for an output that is
  /// not a 64-bit floating-point value.
  DerivedProgram(expr::Program program,
                 const std::vector<DerivedSource>& sources,
                 const DerivedGeometry& geometry);

  [[nodiscard]] const expr::Program& program() const { return program_; }

  /// Output points per cell: subcells times points per subcell.
  [[nodiscard]] std::size_t pointsPerCell() const { return pointsPerCell_; }
  [[nodiscard]] std::size_t pointsPerSubcell() const { return pointsPerSubcell_; }

  /// The sources the program reads, as indices into the list it was prepared with, in the order
  /// of Program::blocks() -- i.e. block i of the program is source usedSources()[i].
  [[nodiscard]] const std::vector<std::size_t>& usedSources() const { return usedSources_; }

  [[nodiscard]] bool readsJacobian() const { return readsJacobian_; }
  [[nodiscard]] bool readsCoordinates() const { return readsCoordinates_; }
  [[nodiscard]] bool readsTime() const { return readsTime_; }
  [[nodiscard]] bool readsTimeStep() const { return readsTimeStep_; }

  /// The output points of a cell in its reference coordinates, subcell by subcell.
  [[nodiscard]] const std::vector<std::array<double, 3>>& referencePoints() const {
    return referencePoints_;
  }

  /// Binds the projection matrices, which this object owns, to `table`.
  void bindMatrices(reader::scripting::DataTable& table) const;

  private:
  struct Matrix {
    std::string name;
    std::vector<double> values;
    std::size_t cols{0};
  };

  expr::Program program_;
  std::size_t pointsPerSubcell_{0};
  std::size_t pointsPerCell_{0};
  std::vector<std::size_t> usedSources_;
  std::vector<Matrix> matrices_;
  std::vector<std::array<double, 3>> referencePoints_;
  bool readsJacobian_{false};
  bool readsCoordinates_{false};
  bool readsTime_{false};
  bool readsTimeStep_{false};
};

/// Whether `name` reads the time integral of the solution, in the vocabulary above.
[[nodiscard]] bool readsTimeIntegral(const std::string& name);

/// What the wave field output of a configuration writes, with the quantities of its material.
struct WaveFieldSelection {
  /// The quantities of the material, and where the velocities start among them.
  std::vector<std::string> quantities;
  std::size_t velocityOffset{0};
  /// Per quantity: write its value, and its time integral (`int-<q>`).
  std::vector<bool> outputMask;
  std::vector<bool> integrationMask;
  /// The strain (epsxx ... epsxz, from the time integral of the velocity) and the rotation (rot1
  /// ... rot3, from the velocity).
  bool strain{false};
  bool rotation{false};
  /// The plastic strain quantities and which of them to write; empty without plasticity.
  std::vector<std::string> plasticQuantities;
  std::vector<bool> plasticityMask;
};

/// The built-in outputs of the wave field as one derived program, in the order and with the names
/// of the hand-written outputs they replace, and with their operations in the same order.
[[nodiscard]] expr::Program waveFieldProgram(const WaveFieldSelection& selection);

/// Whether the channel `name` reads the time integral of the solution, in the vocabulary above.
[[nodiscard]] bool readsTimeIntegral(const std::string& name);

/// Loads a derived-output program from a file: `sderiv:` or a .sderiv file is compiled, `lua:`
/// or a .lua file is traced. A failure is a configuration error.
[[nodiscard]] expr::Program loadDerivedProgram(const std::string& path);

/// Whether the configured wave field output reads the time integral of the solution, which then
/// has to be kept.
[[nodiscard]] bool waveFieldOutputReadsIntegrals(seissol::SeisSol& seissolInstance);

/// The derived outputs of the cells of one configuration, evaluated into a buffer at once and
/// handed to the writer from there.
class DerivedOutput {
  public:
  DerivedOutput() = default;
  virtual ~DerivedOutput() = default;
  DerivedOutput(const DerivedOutput&) = delete;
  DerivedOutput& operator=(const DerivedOutput&) = delete;
  DerivedOutput(DerivedOutput&&) = delete;
  DerivedOutput& operator=(DerivedOutput&&) = delete;

  /// Output names, with the simulation suffix for fused simulations.
  [[nodiscard]] virtual const std::vector<std::string>& names() const = 0;

  /// Whether the program keeps state, and hence wants to see every time step (see step()).
  [[nodiscard]] virtual bool accumulates() const = 0;

  /// Brings the values of all written points to time `time`, before they are written.
  virtual void write(double time) = 0;

  /// Evaluates the written points of the layer `layerId` after a time step of its cluster to time
  /// `time`. Only for a program that accumulates.
  virtual void step(std::size_t layerId, double time) = 0;

  /// Copies the values of output `output` at the points of `subcell` of the written cell `index`.
  virtual void
      copy(std::size_t output, double* target, std::size_t index, std::size_t subcell) const = 0;
};

/// The derived outputs of `program` at the cells `cells` (mesh ids, in the order written) of the
/// configuration `Cfg`; cells of other configurations are skipped. Null if none of the cells is
/// of the configuration `Cfg`.
template <typename Cfg>
std::shared_ptr<DerivedOutput>
    makeDerivedVolumeOutput(seissol::SeisSol& seissolInstance,
                            const std::shared_ptr<const std::vector<std::size_t>>& cells,
                            const DerivedGeometry& geometry,
                            const expr::Program& program);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_DERIVEDOUTPUT_H_
