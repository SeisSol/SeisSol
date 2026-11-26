// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WavePropagation/Solver/PlaneWaveOperator.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "WavePropagation/Physics/CommonDatastructures.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Solver/MaterialOperators.h"

#include <Eigen/Dense>
#include <complex>
#include <vector>

namespace seissol::solver {

std::vector<std::complex<double>>
    planeWaveOperator(ConfigId config, const model::Material& material, const Eigen::Vector3d& k) {
  return dispatchConfig(config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using MaterialT = model::MaterialOf<Cfg>;
    constexpr auto NumQuantities = MaterialT::NumQuantities;

    std::vector<std::complex<double>> matrix(NumQuantities * NumQuantities);
    model::getPlaneWaveOperator<Cfg>(
        dynamic_cast<const MaterialT&>(material), k.data(), matrix.data());
    return matrix;
  });
}

} // namespace seissol::solver
