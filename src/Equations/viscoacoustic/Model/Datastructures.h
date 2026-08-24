// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_DATASTRUCTURES_H_

#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/acoustic/Model/Datastructures.h"
#include "Equations/viscoacoustic/Model/Attenuation.h"
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
struct ViscoAcousticLocalDataSplit;
struct ViscoAcousticNeighborDataSplit;
struct ViscoAcousticLocalDataExtend;
struct ViscoAcousticNeighborDataExtend;

// TODO: unify by solver class
template <ViscoImplementation Implementation>
struct ViscoSolverAcoustic {
  using Type = kernels::solver::linearck::Solver;
  using LocalData = std::monostate;
  using NeighborData = std::monostate;
};

template <>
struct ViscoSolverAcoustic<ViscoImplementation::QuantityExtension> {
  using Type = kernels::solver::linearck::Solver;
  using LocalData = ViscoAcousticLocalDataExtend;
  using NeighborData = ViscoAcousticNeighborDataExtend;
};

template <>
struct ViscoSolverAcoustic<ViscoImplementation::AnelasticTensor> {
  using Type = kernels::solver::linearckanelastic::Solver;
  using LocalData = ViscoAcousticLocalDataSplit;
  using NeighborData = ViscoAcousticNeighborDataSplit;
};

template <std::size_t MechanismsP>
struct ViscoAcousticMaterial : public AcousticMaterial {
  static constexpr std::size_t NumberPerMechanism = 1;
  static constexpr std::size_t NumElasticQuantities = 4;
  static constexpr std::size_t NumQuantities =
      NumElasticQuantities + MechanismsP * NumberPerMechanism;
  static constexpr std::size_t TractionQuantities = 1;
  static constexpr std::size_t Mechanisms = MechanismsP;
  static constexpr MaterialType Type = MaterialType::Viscoacoustic;
  static inline const std::string Text = "viscoacoustic-" + std::to_string(MechanismsP);
  static inline const std::array<std::string, NumElasticQuantities> Quantities{
      "pprime", "v1", "v2", "v3"};
  static constexpr std::size_t Parameters = AcousticMaterial::Parameters + 2 * Mechanisms;

  static constexpr bool SupportsDR = false;
  static constexpr bool SupportsLTS = true;
  static constexpr bool SupportsEnergy = true;

  using LocalSpecificData = ViscoSolverAcoustic<Config::ViscoMode>::LocalData;
  using NeighborSpecificData = ViscoSolverAcoustic<Config::ViscoMode>::NeighborData;

  using Solver = ViscoSolverAcoustic<Config::ViscoMode>::Type;

  using EnergyData = std::monostate;

  //! Relaxation frequencies
  double omega[zeroLengthArrayHandler(Mechanisms)]{};
  /** Entries of the source matrix (E)
   * theta[0] = -lambda * Y_lambda
   **/
  double theta[zeroLengthArrayHandler(Mechanisms)][1]{};
  double qp{};

  static const std::unordered_map<std::string, double ViscoAcousticMaterial::*> ParameterMap;

  ViscoAcousticMaterial() = default;
  explicit ViscoAcousticMaterial(const std::vector<double>& materialValues)
      : AcousticMaterial(materialValues) {
    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      this->omega[mech] = materialValues.at(2 + 2 * mech);
      for (std::size_t i = 0; i < 3; ++i) {
        this->theta[mech][i] = materialValues.at(3 + i + 2 * mech);
      }
    }
    // This constructor is used to initialize a ViscoAcousticMaterial
    // from the values in Fortran. Qp and Qs are not part of the
    // material in Fortran, so we set these to NaN.
    qp = std::numeric_limits<double>::signaling_NaN();
  }

  explicit ViscoAcousticMaterial(const AcousticMaterial& acoustic) : AcousticMaterial(acoustic) {}

  ~ViscoAcousticMaterial() override = default;

  /**
   * Modulus defect of relaxation mechanism `mech`, i.e. the stiffness of the
   * spring in that Maxwell branch.
   *
   * fitAttenuation keeps varP/varMu/alpha/beta as locals and stores only theta,
   * but the defects are exactly recoverable:
   *   theta[0] = -dLambda
   * These are accessors rather than stored fields on purpose -- the material has
   * two construction paths (easi + fitAttenuation, and the legacy vector
   * constructor that reads theta directly), and both populate theta.
   */
  [[nodiscard]] double getDeltaLambda(std::size_t mech) const { return -theta[mech][0]; }
  /// Bulk modulus defect, dLambda + 2/3 dMu.
  [[nodiscard]] double getDeltaBulk(std::size_t mech) const { return -theta[mech][0] / 3.0; }

  /// lambda/mu are the *unrelaxed* moduli once fitAttenuation has run.
  [[nodiscard]] double getLambdaUnrelaxed() const { return lambda; }

  [[nodiscard]] double getLambdaRelaxed() const {
    double result = lambda;
    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      result -= getDeltaLambda(mech);
    }
    return result;
  }
  [[nodiscard]] double getBulkRelaxed() const { return getLambdaRelaxed(); }

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
      if (getDeltaBulk(mech) <= 0 || omega[mech] <= 0) {
        return false;
      }
    }
    return getBulkRelaxed() > 0;
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
inline const std::unordered_map<std::string, double ViscoAcousticMaterial<N>::*>
    ViscoAcousticMaterial<N>::ParameterMap{{"rho", &ViscoAcousticMaterial<N>::rho},
                                           {"lambda", &ViscoAcousticMaterial<N>::lambda},
                                           {"Qp", &ViscoAcousticMaterial<N>::qp}};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_DATASTRUCTURES_H_
