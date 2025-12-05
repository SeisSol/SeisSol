// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_DATASTRUCTURES_H_

#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/PreProcessorMacros.h"
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Model/CommonDatastructures.h"
#include "Physics/Attenuation.h"

#include <Common/Typedefs.h>
#include <Equations/viscoelastic/Model/IntegrationData.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Kernels/LinearCK/Solver.h>
#include <Kernels/LinearCKAnelastic/Solver.h>
#include <Physics/Attenuation.h>
#include <array>
#include <cstddef>
#include <limits>
#include <string>

namespace seissol::model {
template <typename>
struct ViscoElasticQELocalData;
template <typename>
struct ViscoElasticATLocalData;
template <typename>
struct ViscoElasticATNeighborData;

template <ViscoImplementation Implementation>
struct ViscoSolver {
  using Type = kernels::solver::linearck::Solver;

  template <typename Cfg>
  using LocalData = ViscoElasticQELocalData<Cfg>;

  template <typename Cfg>
  using NeighborData = std::monostate;
};

template <>
struct ViscoSolver<ViscoImplementation::QuantityExtension> {
  using Type = kernels::solver::linearck::Solver;

  template <typename Cfg>
  using LocalData = ViscoElasticQELocalData<Cfg>;

  template <typename Cfg>
  using NeighborData = std::monostate;
};

template <>
struct ViscoSolver<ViscoImplementation::AnelasticTensor> {
  using Type = kernels::solver::linearckanelastic::Solver;

  template <typename Cfg>
  using LocalData = ViscoElasticATLocalData<Cfg>;

  template <typename Cfg>
  using NeighborData = ViscoElasticATNeighborData<Cfg>;
};

template <std::size_t MechanismsP>
struct ViscoElasticMaterialParametrized : public ElasticMaterial {
  static constexpr std::size_t NumberPerMechanism = 6;
  static constexpr std::size_t NumElasticQuantities = 9;
  static constexpr std::size_t NumQuantities =
      NumElasticQuantities + MechanismsP * NumberPerMechanism;
  static constexpr std::size_t TractionQuantities = 6;
  static constexpr std::size_t Mechanisms = MechanismsP;
  static constexpr MaterialType Type = MaterialType::Viscoelastic;
  static inline const std::string Text = "viscoelastic-" + std::to_string(MechanismsP);
  static inline const std::array<std::string, NumElasticQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};
  static constexpr std::size_t Parameters = ElasticMaterial::Parameters + 4 * Mechanisms;

  static constexpr bool SupportsDR = true;
  static constexpr bool SupportsLTS = true;

  template <typename Config>
  using LocalSpecificData = typename ViscoSolver<Config::ViscoMode>::template LocalData<Config>;

  template <typename Config>
  using NeighborSpecificData =
      typename ViscoSolver<Config::ViscoMode>::template NeighborData<Config>;

  template <typename Config>
  using Solver = typename ViscoSolver<Config::ViscoMode>::Type;

  //! Relaxation frequencies
  double omega[zeroLengthArrayHandler(Mechanisms)]{};
  /** Entries of the source matrix (E)
   * theta[0] = -(lambda * Y_lambda + 2.0 * mu * Y_mu)
   * theta[1] = -lambda * Y_lambda
   * theta[2] = -2.0 * mu * Y_mu
   **/
  double theta[zeroLengthArrayHandler(Mechanisms)][3]{};
  double qp{};
  double qs{};

  static const std::unordered_map<std::string, double ViscoElasticMaterialParametrized::*>
      ParameterMap;

  ViscoElasticMaterialParametrized() = default;
  explicit ViscoElasticMaterialParametrized(const std::vector<double>& materialValues)
      : ElasticMaterial(materialValues) {
    for (int mech = 0; mech < Mechanisms; ++mech) {
      this->omega[mech] = materialValues.at(3 + 4 * mech);
      for (unsigned i = 1; i < 4; ++i) {
        this->theta[mech][i - 1] = materialValues.at(3 + 4 * mech + i);
      }
    }
    // This constructor is used to initialize a ViscoElasticMaterial
    // from the values in Fortran. Qp and Qs are not part of the
    // material in Fortran, so we set these to NaN.
    qp = std::numeric_limits<double>::signaling_NaN();
    qs = std::numeric_limits<double>::signaling_NaN();
  }

  ~ViscoElasticMaterialParametrized() override = default;

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }

  // for now, keep per material entry
  double maxTimestep{std::numeric_limits<double>::infinity()};

  void initialize(const initializer::parameters::ModelParameters& parameters) override {
    maxTimestep = 0.25 / (parameters.freqCentral * std::sqrt(parameters.freqRatio));

    physics::fitAttenuation<Mechanisms>(*this, parameters.freqCentral, parameters.freqRatio);
  }

  [[nodiscard]] double maximumTimestep() const override { return maxTimestep; }
};

template <std::size_t N>
inline const std::unordered_map<std::string, double ViscoElasticMaterialParametrized<N>::*>
    ViscoElasticMaterialParametrized<N>::ParameterMap{
        {"rho", &ViscoElasticMaterialParametrized<N>::rho},
        {"lambda", &ViscoElasticMaterialParametrized<N>::lambda},
        {"mu", &ViscoElasticMaterialParametrized<N>::mu},
        {"Qp", &ViscoElasticMaterialParametrized<N>::qp},
        {"Qs", &ViscoElasticMaterialParametrized<N>::qs}};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_DATASTRUCTURES_H_
