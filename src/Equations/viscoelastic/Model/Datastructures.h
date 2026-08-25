// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_DATASTRUCTURES_H_

#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/viscoelastic/Model/Attenuation.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/PreProcessorMacros.h"
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Model/CommonDatastructures.h"

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <string>
#include <utils/logger.h>
#include <variant>

namespace seissol::model {

template <ViscoImplementation Implementation>
struct ViscoSolver {
  using Type = kernels::solver::linearck::Solver;
};

template <>
struct ViscoSolver<ViscoImplementation::QuantityExtension> {
  using Type = kernels::solver::linearck::Solver;
};

template <>
struct ViscoSolver<ViscoImplementation::AnelasticTensor> {
  using Type = kernels::solver::linearckanelastic::Solver;
};

template <std::size_t MechanismsP>
struct ViscoElasticMaterial : public ElasticMaterial {
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
  static constexpr bool SupportsEnergy = true;

  static constexpr ViscoImplementation ViscoMode = Config::ViscoMode;

  using LocalSpecificData = std::monostate;
  using NeighborSpecificData = std::monostate;

  using Solver = ViscoSolver<ViscoMode>::Type;

  using EnergyData = std::monostate;

  //! Relaxation frequencies
  double omega[zeroGuard(Mechanisms)]{};
  /** Entries of the source matrix (E)
   * theta[0] = -(lambda * Y_lambda + 2.0 * mu * Y_mu)
   * theta[1] = -lambda * Y_lambda
   * theta[2] = -2.0 * mu * Y_mu
   **/
  double theta[zeroGuard(Mechanisms)][3]{};
  double qp{};
  double qs{};

  static const std::unordered_map<std::string, double ViscoElasticMaterial::*> ParameterMap;

  ViscoElasticMaterial() = default;
  explicit ViscoElasticMaterial(const std::vector<double>& materialValues)
      : ElasticMaterial(materialValues) {
    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      this->omega[mech] = materialValues.at(3 + 4 * mech);
      for (std::size_t i = 0; i < 3; ++i) {
        this->theta[mech][i] = materialValues.at(4 + i + 4 * mech);
      }
    }
    // This constructor is used to initialize a ViscoElasticMaterial
    // from the values in Fortran. Qp and Qs are not part of the
    // material in Fortran, so we set these to NaN.
    qp = std::numeric_limits<double>::signaling_NaN();
    qs = std::numeric_limits<double>::signaling_NaN();
  }

  explicit ViscoElasticMaterial(const ElasticMaterial& elastic) : ElasticMaterial(elastic) {}
  explicit ViscoElasticMaterial(const AcousticMaterial& acoustic) : ElasticMaterial(acoustic) {}

  ~ViscoElasticMaterial() override = default;

  /**
   * Modulus defect of relaxation mechanism `mech`, i.e. the stiffness of the
   * spring in that Maxwell branch.
   *
   * fitAttenuation keeps varP/varMu/alpha/beta as locals and stores only theta,
   * but the defects are exactly recoverable:
   *   theta[0] = -(dLambda + 2 dMu),  theta[1] = -dLambda,  theta[2] = -2 dMu
   * These are accessors rather than stored fields on purpose -- the material has
   * two construction paths (easi + fitAttenuation, and the legacy vector
   * constructor that reads theta directly), and both populate theta.
   */
  [[nodiscard]] double getDeltaLambda(std::size_t mech) const { return -theta[mech][1]; }
  [[nodiscard]] double getDeltaMu(std::size_t mech) const { return -0.5 * theta[mech][2]; }
  /// Bulk modulus defect, dLambda + 2/3 dMu.
  [[nodiscard]] double getDeltaBulk(std::size_t mech) const {
    return -theta[mech][1] - theta[mech][2] / 3.0;
  }

  /// lambda/mu are the *unrelaxed* moduli once fitAttenuation has run.
  [[nodiscard]] double getLambdaUnrelaxed() const { return lambda; }
  [[nodiscard]] double getMuUnrelaxed() const { return mu; }

  [[nodiscard]] double getLambdaRelaxed() const {
    double result = lambda;
    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      result -= getDeltaLambda(mech);
    }
    return result;
  }
  [[nodiscard]] double getMuRelaxed() const {
    double result = mu;
    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      result -= getDeltaMu(mech);
    }
    return result;
  }
  [[nodiscard]] double getBulkRelaxed() const {
    return getLambdaRelaxed() + 2.0 / 3.0 * getMuRelaxed();
  }

  /**
   * Whether the fitted mechanisms form a physically admissible generalized
   * Maxwell body: every branch stiffness positive definite, and the relaxed
   * moduli still positive.
   *
   * alpha/beta in fitAttenuation come out of a least-squares solve and are not
   * sign-constrained, so a poor Q/frequency-band combination can produce a fit
   * that reproduces the target Q but is thermodynamically inconsistent. The
   * energy output would then report a negative viscous dissipation.
   */
  [[nodiscard]] bool attenuationWellPosed() const {
    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      if (getDeltaMu(mech) <= 0 || getDeltaBulk(mech) <= 0 || omega[mech] <= 0) {
        return false;
      }
      // theta[0] must be the sum of the other two; a violation means theta was
      // not produced by the documented parametrization at all
      if (std::abs(theta[mech][0] - (theta[mech][1] + theta[mech][2])) >
          1e-10 * std::max(1.0, std::abs(theta[mech][0]))) {
        return false;
      }
    }
    return getMuRelaxed() > 0 && getBulkRelaxed() > 0;
  }

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }

  void initialize(const initializer::parameters::ModelParameters& parameters) override {
    physics::fitAttenuation<Mechanisms>(*this, parameters.freqCentral, parameters.freqRatio);

    if constexpr (Mechanisms > 0) {
      if (!attenuationWellPosed()) {
        // initialize() runs per cell, and from an OpenMP loop -- warn once
        static std::atomic<bool> warned{false};
        if (!warned.exchange(true)) {
          logWarning() << "The attenuation fit produced at least one mechanism that is not"
                       << "positive definite, or non-positive relaxed moduli. The fit may still"
                       << "reproduce the target Q, but it is thermodynamically inconsistent:"
                       << "the reported viscous dissipation can become negative. Consider a"
                       << "different central frequency or frequency ratio. Further occurrences"
                       << "are not reported.";
        }
      }
    }
  }
};

template <std::size_t N>
inline const std::unordered_map<std::string, double ViscoElasticMaterial<N>::*>
    ViscoElasticMaterial<N>::ParameterMap{{"rho", &ViscoElasticMaterial<N>::rho},
                                          {"lambda", &ViscoElasticMaterial<N>::lambda},
                                          {"mu", &ViscoElasticMaterial<N>::mu},
                                          {"Qp", &ViscoElasticMaterial<N>::qp},
                                          {"Qs", &ViscoElasticMaterial<N>::qs}};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_DATASTRUCTURES_H_
