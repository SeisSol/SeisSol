// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_DERIVEDOUTPUT_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_DERIVEDOUTPUT_H_

// Derived quantities at the output points of the wave field and of the free surface, computed by
// an expression program.
//
// The program is written pointwise. It reads the quantities by name, and the consumer decides
// where each of them comes from; here, from the coefficients of the cell:
//
//   <q>                 quantity q of the solution, e.g. v1 or s_xx
//   <q>_r0, _r1, _r2    its derivative along the reference coordinates xi, eta, zeta of the cell
//   dx_<q>, dy_<q>, dz_<q>
//                       its derivative in space, by the chain rule over the inverse Jacobian
//   int_<q>, int_<q>_rk, d*_int_<q>
//                       the same for the time integral of q
//   <p>, <p>_rk, d*_<p> the plastic strain quantity p, e.g. ep_xx (with plasticity)
//   u1, u2, u3          the displacement of a free surface face (surface output only)
//   jinv<k><d>          d xi_k / d x_d, the inverse Jacobian of the cell
//   x, y, z             the output point
//   t, dt               the time, and the time since the previous evaluation of the point
//
// Every value and reference derivative becomes a contraction of the coefficients of its quantity
// against a projection matrix (expr::Kind::Contract), so the program never sees a basis function
// and the coefficients are not projected into a buffer first. Interning makes a contraction one
// node however many outputs read it: the strain and rotation outputs together need the nine
// reference derivatives of the velocities once, not once per output. The coordinates are a
// contraction as well, of the affine map of the element against its reference points, so that a
// program reads nothing that a device kernel could not.
//
// The subcells of a refined output are stacked into the rows of one matrix: the points of an
// element are its subcells times the output points of a subcell, so one program and one call
// cover them all, and a state slot (a maximum over time, say) exists once per output point. The
// faces of a cell come in four variants, one per side, with projections of their own; the points
// of a call are faces of one side, and the matrices are passed per call like the blocks.
//
// CADENCE. A program without state is evaluated when the output is written, for all written
// points at once, into a buffer the writer copies from. A program with state accumulates: the
// points of a layer are evaluated after every time step of its cluster, and the buffer then
// always holds the running values. A layer whose cluster did not step since the last write (at
// the start, or one that runs on a device) is evaluated when the output is written.
//
// STATE. The state of the configured program lives with the elements, in the storage of the
// cells (LTS::DerivedState) or of the faces (SurfaceLTS::DerivedState): per element, a value per
// declared state, output point and fused simulation. It is set to the initial values before a
// checkpoint is loaded, and written to the checkpoints under a name that changes with its layout,
// so that a restart continues it when the layout is the same and starts it over otherwise.

#include "Expr/Program.h"
#include "IO/Instance/Geometry/Refinement.h"
#include "Numerical/Projection.h"
#include "Reader/Scripting/DataTable.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace seissol {
class SeisSol;
namespace io::instance::checkpoint {
class CheckpointManager;
} // namespace io::instance::checkpoint
namespace initializer::parameters {
struct FreeSurfaceOutputParameters;
struct WaveFieldOutputParameters;
} // namespace initializer::parameters
} // namespace seissol

namespace seissol::initializer {

/// How a quantity a derived program reads is stored per element.
enum class Representation : std::uint8_t {
  /// Modal coefficients of the volume basis of the cell.
  Modal,
  /// Point values at the nodal set of the volume (the plastic strain).
  Nodal,
  /// Point values at the nodes of the face basis (the displacement of a free surface face).
  FaceNodal
};

/// A quantity a derived program can read.
struct DerivedSource {
  std::string name;
  Representation representation{Representation::Modal};
};

/// Where the output points of the wave field lie: the subcells of a cell and the points of a
/// subcell, both in reference coordinates, and how the solution is transferred onto them.
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

/// The same for the free surface: the subcells of a face and the points of a subcell, in the
/// reference coordinates of the face.
struct DerivedSurfaceGeometry {
  std::vector<io::instance::geometry::Subcell<2>> subcells{io::instance::geometry::unrefined<2>()};
  std::vector<std::array<double, 2>> dataBase;
  std::size_t dataOrder{0};
  std::size_t order{1};
  numerical::projection::Target target{numerical::projection::Target::Project};
  /// The nodal set of the volume (the plastic strain); the displacement of a face is always on
  /// the warp&blend nodes of the face.
  numerical::projection::NodalSet nodalSet{numerical::projection::NodalSet::WarpBlend};
};

/// The output points of an element (a cell, or a face of a cell) as the projections onto them see
/// them.
class DerivedPoints {
  public:
  DerivedPoints() = default;
  virtual ~DerivedPoints() = default;
  DerivedPoints(const DerivedPoints&) = delete;
  DerivedPoints& operator=(const DerivedPoints&) = delete;
  DerivedPoints(DerivedPoints&&) = delete;
  DerivedPoints& operator=(DerivedPoints&&) = delete;

  /// Output points per element: subcells times points per subcell.
  [[nodiscard]] virtual std::size_t pointsPerElement() const = 0;
  [[nodiscard]] virtual std::size_t pointsPerSubcell() const = 0;

  /// Elements of different variants -- the side of its cell a face lies on -- have projections of
  /// their own.
  [[nodiscard]] virtual std::size_t variants() const = 0;

  /// Coefficients per element of a source of `representation`; 0 if an element cannot read it.
  [[nodiscard]] virtual std::size_t columns(Representation representation) const = 0;

  /// Whether a source of `representation` has reference derivatives at the points.
  [[nodiscard]] virtual bool derivatives(Representation representation) const = 0;

  /// The projection of a source of `representation`, or of its derivative along the reference
  /// direction `derivative` of the cell, onto the points of an element of `variant`: entry
  /// (point, m) at point + m * pointsPerElement().
  [[nodiscard]] virtual std::vector<double> matrix(std::size_t variant,
                                                   Representation representation,
                                                   std::optional<std::size_t> derivative) const = 0;

  /// The output points of an element in its reference coordinates -- three for a cell, two (and
  /// a zero) for a face -- subcell by subcell.
  [[nodiscard]] virtual std::vector<std::array<double, 3>> referencePoints() const = 0;
};

/// The output points in the cells of the wave field output.
class VolumePoints final : public DerivedPoints {
  public:
  explicit VolumePoints(DerivedGeometry geometry);

  [[nodiscard]] std::size_t pointsPerElement() const override;
  [[nodiscard]] std::size_t pointsPerSubcell() const override;
  [[nodiscard]] std::size_t variants() const override { return 1; }
  [[nodiscard]] std::size_t columns(Representation representation) const override;
  [[nodiscard]] bool derivatives(Representation representation) const override;
  [[nodiscard]] std::vector<double> matrix(std::size_t variant,
                                           Representation representation,
                                           std::optional<std::size_t> derivative) const override;
  [[nodiscard]] std::vector<std::array<double, 3>> referencePoints() const override;

  private:
  DerivedGeometry geometry_;
};

/// The output points on the faces of the free surface output; variant f is a face on side f of
/// its cell.
class SurfacePoints final : public DerivedPoints {
  public:
  explicit SurfacePoints(DerivedSurfaceGeometry geometry);

  [[nodiscard]] std::size_t pointsPerElement() const override;
  [[nodiscard]] std::size_t pointsPerSubcell() const override;
  [[nodiscard]] std::size_t variants() const override;
  [[nodiscard]] std::size_t columns(Representation representation) const override;
  [[nodiscard]] bool derivatives(Representation representation) const override;
  [[nodiscard]] std::vector<double> matrix(std::size_t variant,
                                           Representation representation,
                                           std::optional<std::size_t> derivative) const override;
  [[nodiscard]] std::vector<std::array<double, 3>> referencePoints() const override;

  private:
  DerivedSurfaceGeometry geometry_;
};

/// A derived-output program prepared for the output points of an element: its channels resolved
/// against the vocabulary above, the contractions in place, and the projection matrices built.
class DerivedProgram {
  public:
  /// Throws std::invalid_argument for a channel outside the vocabulary, for one the elements
  /// cannot read, and for an output that is not a 64-bit floating-point value.
  DerivedProgram(expr::Program program,
                 const std::vector<DerivedSource>& sources,
                 const DerivedPoints& points);

  /// The same for the cells of the wave field.
  DerivedProgram(expr::Program program,
                 const std::vector<DerivedSource>& sources,
                 const DerivedGeometry& geometry);

  [[nodiscard]] const expr::Program& program() const { return program_; }

  /// Output points per element: subcells times points per subcell.
  [[nodiscard]] std::size_t pointsPerElement() const { return pointsPerElement_; }
  [[nodiscard]] std::size_t pointsPerSubcell() const { return pointsPerSubcell_; }

  /// The sources the program reads, as indices into the list it was prepared with, in the order
  /// of Program::blocks(): block i is source usedSources()[i]. If the program reads the
  /// coordinates, the blocks end with the affine map of the element along x, y and z (see
  /// bindGeometry()).
  [[nodiscard]] const std::vector<std::size_t>& usedSources() const { return usedSources_; }

  [[nodiscard]] bool readsJacobian() const { return readsJacobian_; }
  [[nodiscard]] bool readsCoordinates() const { return readsCoordinates_; }
  [[nodiscard]] bool readsTime() const { return readsTime_; }
  [[nodiscard]] bool readsTimeStep() const { return readsTimeStep_; }

  /// Binds the projection matrices of the first variant, which this object owns, to `table`.
  void bindMatrices(reader::scripting::DataTable& table) const;

  /// The bases of the matrices of `variant`, in the order of Program::matrices(), for a call.
  [[nodiscard]] std::vector<const void*> matrixBases(std::size_t variant) const;

  /// Binds the affine maps of the elements if the program reads the coordinates: 12 per element
  /// in the order of the point set, the origin, then the images of the reference unit vectors less
  /// the origin (zero for the third of a face).
  void bindGeometry(reader::scripting::DataTable& table, const double* transforms) const;

  private:
  struct Matrix {
    std::string name;
    /// One per variant.
    std::vector<std::vector<double>> values;
    std::size_t cols{0};
  };

  expr::Program program_;
  std::size_t pointsPerSubcell_{0};
  std::size_t pointsPerElement_{0};
  std::vector<std::size_t> usedSources_;
  std::vector<Matrix> matrices_;
  bool readsJacobian_{false};
  bool readsCoordinates_{false};
  bool readsTime_{false};
  bool readsTimeStep_{false};
};

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

/// The built-in outputs of the free surface as one derived program: the quantities of
/// `outputMask` from the cell, and the displacement u1, u2, u3 of the face.
[[nodiscard]] expr::Program surfaceProgram(const std::vector<std::string>& quantities,
                                           const std::vector<bool>& outputMask);

/// The output points of the wave field output as configured: the subcells of a cell, the points
/// of a subcell and the projection onto them. The convergence order and the nodal set are the
/// configuration's, and left to the caller.
[[nodiscard]] DerivedGeometry
    waveFieldGeometry(const parameters::WaveFieldOutputParameters& parameters);

/// The same for the free surface output.
[[nodiscard]] DerivedSurfaceGeometry
    surfaceGeometry(const parameters::FreeSurfaceOutputParameters& parameters);

/// The output a derived program serves.
enum class DerivedOutputKind : std::uint8_t { WaveField, Surface };

/// The state the configured program of an output keeps per element (see STATE above): in an
/// element, value (state, point, simulation) at offset(state, point, simulation).
struct DerivedStateLayout {
  /// The declared states, in the order of Program::state(), and their initial values.
  std::vector<std::string> states;
  std::vector<double> initial;
  std::size_t pointsPerElement{0};
  /// The most fused simulations of a configuration of the run.
  std::size_t simulations{0};
  /// Changes with whatever changes the meaning of a value: the states, the output points, the
  /// simulations.
  std::uint64_t fingerprint{0};

  /// Values per element; 0 without state.
  [[nodiscard]] std::size_t size() const { return states.size() * pointsPerElement * simulations; }
  [[nodiscard]] std::size_t
      offset(std::size_t state, std::size_t point, std::size_t simulation) const {
    return (state * pointsPerElement + point) * simulations + simulation;
  }
  /// The name of the dataset in the checkpoints.
  [[nodiscard]] std::string checkpointName() const;
};

/// The state layout of the configured program of output `kind`; empty if the output is off or its
/// program keeps no state.
[[nodiscard]] DerivedStateLayout derivedStateLayout(seissol::SeisSol& seissolInstance,
                                                    DerivedOutputKind kind);

/// Sets the state of the derived outputs to its initial values, in the storage of the cells and
/// of the faces. Before a checkpoint is loaded.
void initializeDerivedState(seissol::SeisSol& seissolInstance);

/// Registers the state of the derived outputs with the checkpoints, after the trees of the cells
/// and of the faces. A checkpoint without it (or with another layout) leaves the initial values.
void registerDerivedStateCheckpoints(io::instance::checkpoint::CheckpointManager& checkpoint,
                                     seissol::SeisSol& seissolInstance);

/// Whether the channel `name` reads the time integral of the solution, in the vocabulary above.
[[nodiscard]] bool readsTimeIntegral(const std::string& name);

/// Loads a derived-output program from a file: `sderiv:` or a .sderiv file is compiled, `lua:`
/// or a .lua file is traced. A failure is a configuration error.
[[nodiscard]] expr::Program loadDerivedProgram(const std::string& path);

/// Whether the configured outputs of the wave field and of the free surface read the time
/// integral of the solution, which then has to be kept.
[[nodiscard]] bool derivedOutputsReadIntegrals(seissol::SeisSol& seissolInstance);

/// The derived outputs of the elements of one configuration, evaluated into a buffer at once and
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

  /// Copies the values of output `output` at the points of `subcell` of the written element
  /// `index`.
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

/// The derived outputs of `program` at the written faces of the free surface (in the order of the
/// free surface integrator) whose cells are of the configuration `Cfg`. Null if there is none.
template <typename Cfg>
std::shared_ptr<DerivedOutput> makeDerivedSurfaceOutput(seissol::SeisSol& seissolInstance,
                                                        const DerivedSurfaceGeometry& geometry,
                                                        const expr::Program& program);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_DERIVEDOUTPUT_H_
