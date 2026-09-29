// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_ENERGYQUADRATURE_H_
#define SEISSOL_SRC_RESULTWRITER_ENERGYQUADRATURE_H_

#include "Alignment.h"
#include "Common/Constants.h"
#include "Equations/Energy.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Numerical/Quadrature.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace seissol::writer {

namespace detail {
/// Whether a material carries a relaxation, whose theta it derives from its
/// parameters.
template <typename MaterialT, typename = void>
struct HasRelaxation : std::false_type {};
template <typename MaterialT>
struct HasRelaxation<MaterialT, std::void_t<decltype(std::declval<MaterialT&>().theta)>>
    : std::true_type {};
} // namespace detail

/// The volume energies of a cell whose material varies inside it.
///
/// They are integrated point by point, over the volume quadrature a modal field
/// is evaluated at, each with the material at that point. That quadrature is
/// exact for the product of two modal fields, so for a material that does not
/// vary it gives what the cell moments give.
///
/// The material at a point interpolates every parameter the samples were
/// queried for, and the theta of an attenuating material, which is derived
/// from those and varies with them. What does not vary, like the relaxation
/// frequencies, comes from a sample as it is. The weights are those of the
/// samples at the points, which are the samples themselves where they are that
/// set, and otherwise the values the operator is formed from where it is formed
/// at those points.
template <typename MaterialT>
class EnergyQuadrature {
  public:
  using Energy = model::EnergyCompute<MaterialT>;
  static constexpr std::size_t Points = tensor::materialToQuadrature::Shape[0];
  static constexpr std::size_t Samples = tensor::materialToQuadrature::Shape[1];
  /// per simulation, as the energy output sums them up
  static constexpr std::size_t Count = Energy::EnergyCount * multisim::NumSimulations;

  static_assert(Samples == LTS::MaterialNodes,
                "The interpolation to the quadrature starts from other samples than a cell has.");
  static_assert(tensor::dofsQP::Shape[multisim::BasisFunctionDimension] == Points,
                "The quantities are evaluated at other points than the material.");

  EnergyQuadrature() {
    constexpr auto PerDirection = ConvergenceOrder + 1;
    static_assert(PerDirection * PerDirection * PerDirection == Points,
                  "The material is interpolated to another quadrature than the one used here.");
    std::array<std::array<double, 3>, Points> points{};
    seissol::quadrature::TetrahedronQuadrature(points.data(), weights_.data(), PerDirection);

    const auto interpolation =
        init::materialToQuadrature::view::create(init::materialToQuadrature::Values);
    for (std::size_t point = 0; point < Points; ++point) {
      for (std::size_t sample = 0; sample < Samples; ++sample) {
        if (interpolation.isInRange(point, sample) && interpolation(point, sample) != 0) {
          interpolation_[point].emplace_back(sample, interpolation(point, sample));
        }
      }
    }
    for (const auto& parameter : MaterialT::ParameterMap) {
      members_.push_back(parameter.second);
    }
  }

  [[nodiscard]] MaterialT materialAt(const std::array<MaterialT, Samples>& samples,
                                     std::size_t point) const {
    const auto interpolate = [&](const auto& read) {
      double value = 0;
      for (const auto& [sample, weight] : interpolation_[point]) {
        value += weight * read(samples[sample]);
      }
      return value;
    };

    MaterialT material = samples[0];
    for (const auto member : members_) {
      material.*member = interpolate([member](const MaterialT& sample) { return sample.*member; });
    }
    if constexpr (detail::HasRelaxation<MaterialT>::value) {
      using Theta = std::remove_reference_t<decltype(material.theta)>;
      for (std::size_t mech = 0; mech < std::extent_v<Theta, 0>; ++mech) {
        for (std::size_t i = 0; i < std::extent_v<Theta, 1>; ++i) {
          material.theta[mech][i] =
              interpolate([mech, i](const MaterialT& sample) { return sample.theta[mech][i]; });
        }
      }
    }
    return material;
  }

  /// The energies over the reference cell, indexed as simulation times
  /// EnergyCount plus energy. The shear modulus of each point goes to
  /// `shearModulus`, for what else is integrated there.
  std::array<double, Count> energies(const std::array<MaterialT, Samples>& samples,
                                     const real* dofs,
                                     const real* dofsAne,
                                     const seissol::Pool& pool,
                                     std::array<double, Points>& shearModulus) const {
    std::array<double, Count> result{};

    alignas(Alignment) real dofsAtPointsData[tensor::dofsQP::size()];
    kernel::evalAtQP evalKrnl;
    evalKrnl.bindGlobals(pool);
    evalKrnl.Q = dofs;
    evalKrnl.dofsQP = dofsAtPointsData;
    evalKrnl.execute();
    const auto dofsAtPoints = init::dofsQP::view::create(dofsAtPointsData);

    const auto anelastic = Energy::evaluateAnelastic(dofsAne, pool);

    for (std::size_t point = 0; point < Points; ++point) {
      const auto material = materialAt(samples, point);
      const auto data = Energy::initEnergyData(material);
      const auto moments = Energy::pointMoments(dofsAtPointsData, anelastic, point);
      shearModulus[point] = material.getMuBar();

      for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
        const auto values = multisim::simtensor(dofsAtPoints, sim);
        // the value at the point, and the products there, stand in for the
        // moments
        const auto linear = [&](std::size_t /*row*/, std::size_t quantity) -> double {
          return values(point, quantity);
        };
        const auto quadratic = [&](std::size_t first, std::size_t second) -> double {
          return static_cast<double>(values(point, first)) * values(point, second);
        };

        const auto local = Energy::computeEnergies(material, data, linear, quadratic, moments, sim);
        for (std::size_t i = 0; i < local.size(); ++i) {
          result[local.size() * sim + i] += weights_[point] * local[i];
        }
      }
    }
    return result;
  }

  private:
  std::array<double, Points> weights_{};
  std::array<std::vector<std::pair<std::size_t, double>>, Points> interpolation_;
  std::vector<double MaterialT::*> members_;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_ENERGYQUADRATURE_H_
