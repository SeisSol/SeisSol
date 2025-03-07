// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_MODEL_COMMONDATASTRUCTURES_H_
#define SEISSOL_SRC_MODEL_COMMONDATASTRUCTURES_H_

#include <Initializer/Parameters/ModelParameters.h>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace seissol::model {
enum class MaterialType { Solid, Acoustic, Elastic, Viscoelastic, Anisotropic, Poroelastic };

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

struct Material {
  static constexpr std::size_t NumQuantities = 0;             // ?
  static constexpr std::size_t NumberPerMechanism = 0;        // ?
  static constexpr std::size_t TractionQuantities = 0;        // ?
  static constexpr std::size_t Mechanisms = 0;                // ?
  static constexpr MaterialType Type = MaterialType::Solid;   // ?
  static constexpr LocalSolver Solver = LocalSolver::Unknown; // ?
  static inline const std::string Text = "material";
  static inline const std::array<std::string, NumQuantities> Quantities = {};
  static constexpr std::size_t Parameters = 1; // rho

  static constexpr bool VaryingWavespeeds = false;
  static constexpr bool SupportsDR = false;
  static constexpr bool SupportsLTS = false;

  virtual ~Material() = default;

  double rho;
  Material() = default;
  Material(const std::vector<double>& data) : rho(data.at(0)) {}

  virtual void initialize(const initializer::parameters::ModelParameters& parameters) {}

  [[nodiscard]] virtual double getMaxWaveSpeed() const = 0;
  [[nodiscard]] virtual double getPWaveSpeed() const = 0;
  [[nodiscard]] virtual double getSWaveSpeed() const = 0;
  [[nodiscard]] virtual double getRhoBar() const { return rho; }
  [[nodiscard]] virtual double getMuBar() const = 0;
  [[nodiscard]] virtual double getLambdaBar() const = 0;
  [[nodiscard]] virtual double maximumTimestep() const {
    return std::numeric_limits<double>::infinity();
  }
  virtual void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const = 0;
  [[nodiscard]] virtual MaterialType getMaterialType() const = 0;

  // very temporary; if someone has a better solution, remove that one here
  virtual void setMuLambda(double mu, double lambda) {}
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
};
} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_COMMONDATASTRUCTURES_H_
