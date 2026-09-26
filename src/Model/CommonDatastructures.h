// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_MODEL_COMMONDATASTRUCTURES_H_
#define SEISSOL_SRC_MODEL_COMMONDATASTRUCTURES_H_

#include "Initializer/Parameters/ModelParameters.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace seissol::model {

/// Where a scalar coefficient of the operator gets its value, which is what
/// decides whether a cell has to carry it.
enum class CoefficientOrigin : std::uint8_t {
  /// Read off the material. It varies from cell to cell, and within a cell
  /// wherever the material is sampled at the nodal points.
  Material,
  /// One number for the whole domain, fixed once the run is set up. The
  /// relaxation frequencies are the case: they follow the frequency band and
  /// nothing else, so no material parameter moves them.
  Global,
};

/// One entry of a transposed coefficient matrix, written as a scalar
/// coefficient of the material times a constant factor. A material that
/// declares these gives
///
///     A_dim(row, column) = sum_a coefficients[a] * factor
///
/// over the entries carrying that dim, row and column, so that the geometry
/// and the material can be folded into the matrix separately.
struct CoefficientEntry {
  std::size_t coefficient;
  std::size_t dim;
  std::size_t row;
  std::size_t column;
  double factor;
};

/// The same for a source term, which has no direction. One relaxation
/// mechanism contributes the whole table, so the mechanism is not an index
/// here either: a caller walks the table once per mechanism with that
/// mechanism's coefficients.
/// One entry of a relaxation mechanism's coupling block, with the column
/// given relative to that mechanism's block. It carries no coefficient index:
/// the whole block is weighted by one scalar, and which one that is -- the
/// relaxation frequency, or nothing at all -- is the solver's decision.
struct AnelasticCoefficientEntry {
  std::size_t dim;
  std::size_t row;
  std::size_t columnOffset;
  double factor;
};

struct SourceCoefficientEntry {
  std::size_t coefficient;
  std::size_t row;
  std::size_t column;
  double factor;
};
enum class MaterialType {
  Solid,
  Acoustic,
  Elastic,
  Viscoelastic,
  Viscoacoustic,
  Anisotropic,
  Poroelastic
};

// the local solvers. CK is the default for elastic, acoustic etc.
// viscoelastic uses CauchyKovalevskiAnelastic (maybe all other materials may be extended to use
// that one as well) poroelastic uses SpaceTimePredictorPoroelastic (someone may generalize that
// one, but so long I(David) had decided to put poroelastic in its name) the solver Unknown is a
// dummy to let all other implementations fail
enum class LocalSolver {
  Unknown,
  CauchyKovalevski,
  CauchyKovalevskiAnelastic,
  SpaceTimePredictorPoroelastic
};

/**
 * A source row stiff enough that a space-time predictor has to factorise it
 * separately: it carries a damping term on the diagonal and feeds one other
 * quantity through an off-diagonal entry.
 */
struct StiffSourceRow {
  std::size_t quantity;
  std::size_t target;
};

struct Material {
  static constexpr std::size_t NumQuantities = 0;      // ?
  static constexpr std::size_t NumberPerMechanism = 0; // ?
  /// Materials whose source term is not stiff declare none.
  static constexpr std::array<StiffSourceRow, 0> StiffSourceRows{};

  static constexpr std::size_t VelocityOffset = 0;
  static constexpr std::size_t TractionComponents = 0;
  static constexpr std::size_t Mechanisms = 0;                // ?
  static constexpr MaterialType Type = MaterialType::Solid;   // ?
  static constexpr LocalSolver Solver = LocalSolver::Unknown; // ?
  static inline const std::string Text = "material";
  static inline const std::array<std::string, NumQuantities> Quantities = {};
  static constexpr std::size_t Parameters = 1; // rho

  virtual ~Material() = default;

  double rho{};
  Material() = default;
  explicit Material(double rho) : rho(rho) {}
  explicit Material(const std::vector<double>& data) : rho(data.at(0)) {}

  virtual void initialize(const initializer::parameters::ModelParameters& parameters) {}

  [[nodiscard]] virtual double getMaxWaveSpeed() const = 0;
  [[nodiscard]] virtual double getPWaveSpeed() const = 0;
  [[nodiscard]] virtual double getSWaveSpeed() const = 0;
  [[nodiscard]] virtual double getMuBar() const = 0;
  [[nodiscard]] virtual double getLambdaBar() const = 0;
  [[nodiscard]] virtual double getDensity() const { return rho; }
  [[nodiscard]] virtual double maximumTimestep() const {
    return std::numeric_limits<double>::infinity();
  }
  virtual void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const = 0;
  [[nodiscard]] virtual MaterialType getMaterialType() const = 0;

  virtual void setLameParameters(double mu, double lambda) {}
  virtual void setDensity(double rho) { this->rho = rho; }
};

struct Plasticity {
  static const inline std::string Text = "plasticity";
  double bulkFriction;
  double plastCo;
  double sXX;
  double sYY;
  double sZZ;
  double sXY;
  double sYZ;
  double sXZ;

  static const std::unordered_map<std::string, double Plasticity::*> ParameterMap;
};

inline const std::unordered_map<std::string, double Plasticity::*> Plasticity::ParameterMap{
    {"bulkFriction", &Plasticity::bulkFriction},
    {"plastCo", &Plasticity::plastCo},
    {"s_xx", &Plasticity::sXX},
    {"s_yy", &Plasticity::sYY},
    {"s_zz", &Plasticity::sZZ},
    {"s_xy", &Plasticity::sXY},
    {"s_yz", &Plasticity::sYZ},
    {"s_xz", &Plasticity::sXZ},
};

struct IsotropicWaveSpeeds {
  double density;
  double pWaveVelocity;
  double sWaveVelocity;
};
} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_COMMONDATASTRUCTURES_H_
