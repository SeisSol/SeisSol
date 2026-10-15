// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Physics/Scenario/PlaneWave.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Equations/Datastructures.h"
#include "Initializer/Typedefs.h"
#include "Model/Common.h"
#include "Model/MaterialType.h"
#include "Numerical/Eigenvalues.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/ScriptedScenario.h"

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <initializer_list>
#include <math.h>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

// the specializations of the plane wave operator for the equations, which it needs to be correct
// IWYU pragma: no_include "Kernels/LinearCKAnelastic/Setup.h"
#include "Equations/Setup.h" // IWYU pragma: keep

namespace seissol::physics {

namespace {

using scenario::concat;
using scenario::literal;

/// theta of a mode: Re(omega) t - k.(x - o) + phi, in the order in which a complex argument of
/// exp is formed, so that both round alike.
std::string argumentOf(const PlaneWaveMode& mode, const std::optional<Eigen::Vector3d>& origin) {
  const std::array<const char*, 3> coordinates{"x", "y", "z"};
  std::string argument = concat({literal(mode.omega.real()), " * t"});
  for (std::size_t d = 0; d < 3; ++d) {
    const std::string position =
        origin.has_value() ? concat({"(", coordinates[d], " - ", literal((*origin)[d]), ")"})
                           : std::string(coordinates[d]);
    argument.append(concat({" - ", literal(mode.k[d]), " * ", position}));
  }
  return argument.append(" + phase");
}

/// The fields of the module of `sum`, one per fused simulation if `perSimulation`, else one.
std::vector<std::unique_ptr<InitialField>> fieldsOf(const PlaneWaveSum& sum,
                                                    ConfigId config,
                                                    const std::string& name,
                                                    bool perSimulation) {
  return scenario::scriptedFields(
      planeWaveModule(sum, scenario::stateNamesOf(config)), config, name, perSimulation);
}

/// The phase of fused simulation `sim` out of the simulations of `config`: 2 pi sim / simulations.
std::string phasePerSimulation(ConfigId config) {
  return concat({"2.0 * pi * sim / ", std::to_string(configValue(config).numSimulations)});
}

} // namespace

std::vector<PlaneWaveMode> planeWaveModes(const CellMaterialData& materialData,
                                          ConfigId config,
                                          const Eigen::Vector3d& k,
                                          const std::vector<int>& modes,
                                          const std::vector<std::complex<double>>& amplitudes) {
  std::size_t numQuantities = 0;
  std::vector<std::complex<double>> values;
  std::vector<std::complex<double>> vectors;
  dispatchConfig(config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using MaterialT = model::MaterialOf<Cfg>;
    constexpr auto NumQuantities = MaterialT::NumQuantities;

    std::array<std::complex<double>, NumQuantities * NumQuantities> planeWaveOperator{};
    model::getPlaneWaveOperator<Cfg>(
        *dynamic_cast<MaterialT*>(materialData.local), k.data(), planeWaveOperator.data());
    eigenvalues::Eigenpair<std::complex<double>, NumQuantities> eigendecomposition;
    eigenvalues::computeEigenvalues<model::EigenvalueBackend<MaterialT>>(planeWaveOperator,
                                                                         eigendecomposition);
    numQuantities = NumQuantities;
    values.assign(eigendecomposition.values.begin(), eigendecomposition.values.end());
    vectors.assign(eigendecomposition.vectors.begin(), eigendecomposition.vectors.end());
  });

  std::vector<PlaneWaveMode> result;
  for (std::size_t v = 0; v < modes.size(); ++v) {
    const auto mode = static_cast<std::size_t>(modes[v]);
    PlaneWaveMode planeWaveMode{k, values.at(mode), {}};
    for (std::size_t j = 0; j < numQuantities; ++j) {
      // the eigenvectors are the columns of the matrix, which is column-major
      planeWaveMode.coefficients.push_back(vectors.at(j + numQuantities * mode) * amplitudes[v]);
    }
    result.push_back(std::move(planeWaveMode));
  }
  return result;
}

std::pair<std::vector<int>, std::vector<std::complex<double>>>
    defaultPlaneWaveModes(const CellMaterialData& materialData, ConfigId config) {
  const auto materialType = configValue(config).materialType;
  if (materialType == model::MaterialType::Acoustic ||
      materialType == model::MaterialType::Viscoacoustic) {
    // the wave modes of an acoustic material are P, N, N, -P: the P mode
    return {{0}, {1.0}};
  }
  if (materialType == model::MaterialType::Poroelastic) {
    // -P, -S2, -S1, -Ps, N, N, N, N, N, Ps, S1, S2, P: -S1, -Ps and P
    return {{2, 3, 12}, {1.0, 1.0, 1.0}};
  }
  if (materialData.local->getMuBar() <= 1e-15) {
    // an elastic material without shear: -P, N, N, N, N, N, N, N, P: the P mode
    return {{8}, {1.0}};
  }
  // -P, -S2, -S1, N, N, N, S1, S2, P: -S2 and P
  return {{1, 8}, {1.0, 1.0}};
}

std::string planeWaveModule(const PlaneWaveSum& sum, const std::vector<std::string>& names) {
  std::string source = "# a plane wave, written by SeisSol from the eigenpairs of the plane wave "
                       "operator of the material\n";
  source.append(concat({"def phase = ", sum.phase, "\n"}));

  // per component of the state, the sum of each wave
  std::vector<std::vector<std::string>> terms(names.size());
  std::size_t index = 0;
  for (const auto& wave : sum.waves) {
    std::vector<std::vector<std::string>> waveTerms(names.size());
    for (const auto& mode : wave) {
      ++index;
      const auto theta = "th" + std::to_string(index);
      source.append(concat({"def ", theta, " = ", argumentOf(mode, sum.origin), "\n"}));
      std::string cosine = concat({"cos(", theta, ")"});
      std::string sine = concat({"sin(", theta, ")"});
      if (mode.omega.imag() != 0.0) {
        // the damping of an attenuated mode
        const auto damping = "d" + std::to_string(index);
        source.append(concat({"def ", damping, " = exp(", literal(-mode.omega.imag()), " * t)\n"}));
        cosine = concat({"(", damping, " * ", cosine, ")"});
        sine = concat({"(", damping, " * ", sine, ")"});
      }
      for (std::size_t j = 0; j < names.size(); ++j) {
        const auto coefficient = mode.coefficients.at(j);
        std::string term;
        if (coefficient.real() != 0.0) {
          term.append(concat({literal(coefficient.real()), " * ", cosine}));
        }
        if (coefficient.imag() != 0.0) {
          term.append(term.empty() ? literal(-coefficient.imag())
                                   : concat({" - ", literal(coefficient.imag())}));
          term.append(concat({" * ", sine}));
        }
        if (term.empty()) {
          continue;
        }
        if (sum.origin.has_value()) {
          // one period of the mode
          term = concat({"select(land(lt(",
                         literal(-0.5 * M_PI),
                         ", ",
                         theta,
                         "), lt(",
                         theta,
                         ", ",
                         literal(1.5 * M_PI),
                         ")), ",
                         term,
                         ", 0.0)"});
        }
        waveTerms[j].push_back(std::move(term));
      }
    }
    for (std::size_t j = 0; j < names.size(); ++j) {
      if (waveTerms[j].empty()) {
        continue;
      }
      std::string waveSum;
      for (const auto& term : waveTerms[j]) {
        waveSum.append(waveSum.empty() ? "" : " + ").append(term);
      }
      terms[j].push_back(sum.waves.size() > 1 ? concat({"(", waveSum, ")"}) : waveSum);
    }
  }

  for (std::size_t j = 0; j < names.size(); ++j) {
    if (terms[j].empty()) {
      continue;
    }
    std::string value;
    for (const auto& term : terms[j]) {
      value.append(value.empty() ? "" : " + ").append(term);
    }
    source.append(concat({"out def ", names[j], " = ", value, "\n"}));
  }
  return source;
}

std::vector<std::unique_ptr<InitialField>> planarwave(const CellMaterialData& materialData,
                                                      ConfigId config) {
  const auto [modes, amplitudes] = defaultPlaneWaveModes(materialData, config);
  PlaneWaveSum sum;
  sum.waves.push_back(
      planeWaveModes(materialData, config, Eigen::Vector3d(M_PI, M_PI, M_PI), modes, amplitudes));
  sum.phase = phasePerSimulation(config);
  return fieldsOf(sum, config, "the planar wave", true);
}

std::vector<std::unique_ptr<InitialField>>
    superimposedPlanarwave(const CellMaterialData& materialData, ConfigId config) {
  const auto [modes, amplitudes] = defaultPlaneWaveModes(materialData, config);
  PlaneWaveSum sum;
  for (const auto& k : {Eigen::Vector3d(M_PI, 0.0, 0.0),
                        Eigen::Vector3d(0.0, M_PI, 0.0),
                        Eigen::Vector3d(0.0, 0.0, M_PI)}) {
    sum.waves.push_back(planeWaveModes(materialData, config, k, modes, amplitudes));
  }
  sum.phase = phasePerSimulation(config);
  return fieldsOf(sum, config, "the superimposed planar wave", true);
}

std::vector<std::unique_ptr<InitialField>>
    travellingWave(const CellMaterialData& materialData,
                   ConfigId config,
                   const TravellingWaveParameters& parameters) {
  logInfo() << "Impose a travelling wave as initial condition";
  logInfo() << "Origin = (" << parameters.origin[0] << ", " << parameters.origin[1] << ", "
            << parameters.origin[2] << ")";
  logInfo() << "kVec = (" << parameters.kVec[0] << ", " << parameters.kVec[1] << ", "
            << parameters.kVec[2] << ")";
  logInfo() << "Combine following wave modes";
  for (std::size_t i = 0; i < parameters.ampField.size(); ++i) {
    logInfo() << "(" << parameters.varField[i] << ": " << parameters.ampField[i] << ")";
  }

  PlaneWaveSum sum;
  sum.waves.push_back(planeWaveModes(
      materialData, config, parameters.kVec, parameters.varField, parameters.ampField));
  // a zero at the origin; the wave travels in the direction of kVec, its wave length is
  // 2 pi / |kVec|
  sum.phase = literal(0.5 * M_PI);
  sum.origin = parameters.origin;
  return fieldsOf(sum, config, "the travelling wave", false);
}

} // namespace seissol::physics
