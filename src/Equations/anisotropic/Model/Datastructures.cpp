// SPDX-FileCopyrightText: 2019-2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Wolf

#include "Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace seissol::model {
double AnisotropicMaterial::getLambdaBar() const {
  return (c11 + c22 + c33) / 3.0 - 2.0 * getMuBar();
}

double AnisotropicMaterial::getMuBar() const { return (c44 + c55 + c66) / 3.0; }

AnisotropicMaterial::AnisotropicMaterial() = default;

AnisotropicMaterial::AnisotropicMaterial(const ElasticMaterial& m)
    : c11(m.lambda + 2 * m.mu), c12(m.lambda), c13(m.lambda), c14(0), c15(0), c16(0),
      c22(m.lambda + 2 * m.mu), c23(m.lambda), c24(0), c25(0), c26(0), c33(m.lambda + 2 * m.mu),
      c34(0), c35(0), c36(0), c44(m.mu), c45(0), c46(0), c55(m.mu), c56(0), c66(m.mu) {
  rho = m.rho;
}

AnisotropicMaterial::AnisotropicMaterial(const std::vector<double>& materialValues)
    : Material(materialValues), c11(materialValues.at(1)), c12(materialValues.at(2)),
      c13(materialValues.at(3)), c14(materialValues.at(4)), c15(materialValues.at(5)),
      c16(materialValues.at(6)), c22(materialValues.at(7)), c23(materialValues.at(8)),
      c24(materialValues.at(9)), c25(materialValues.at(10)), c26(materialValues.at(11)),
      c33(materialValues.at(12)), c34(materialValues.at(13)), c35(materialValues.at(14)),
      c36(materialValues.at(15)), c44(materialValues.at(16)), c45(materialValues.at(17)),
      c46(materialValues.at(18)), c55(materialValues.at(19)), c56(materialValues.at(20)),
      c66(materialValues.at(21)) {}

AnisotropicMaterial::~AnisotropicMaterial() = default;

void AnisotropicMaterial::getFullStiffnessTensor(std::array<double, 81>& fullTensor) const {
  auto stiffnessTensorView =
      seissol_general::init::stiffnessTensor::view::create(fullTensor.data());
  stiffnessTensorView.setZero();
  stiffnessTensorView(0, 0, 0, 0) = c11;
  stiffnessTensorView(0, 0, 0, 1) = c16;
  stiffnessTensorView(0, 0, 0, 2) = c15;
  stiffnessTensorView(0, 0, 1, 0) = c16;
  stiffnessTensorView(0, 0, 1, 1) = c12;
  stiffnessTensorView(0, 0, 1, 2) = c14;
  stiffnessTensorView(0, 0, 2, 0) = c15;
  stiffnessTensorView(0, 0, 2, 1) = c14;
  stiffnessTensorView(0, 0, 2, 2) = c13;
  stiffnessTensorView(0, 1, 0, 0) = c16;
  stiffnessTensorView(0, 1, 0, 1) = c66;
  stiffnessTensorView(0, 1, 0, 2) = c56;
  stiffnessTensorView(0, 1, 1, 0) = c66;
  stiffnessTensorView(0, 1, 1, 1) = c26;
  stiffnessTensorView(0, 1, 1, 2) = c46;
  stiffnessTensorView(0, 1, 2, 0) = c56;
  stiffnessTensorView(0, 1, 2, 1) = c46;
  stiffnessTensorView(0, 1, 2, 2) = c36;
  stiffnessTensorView(0, 2, 0, 0) = c15;
  stiffnessTensorView(0, 2, 0, 1) = c56;
  stiffnessTensorView(0, 2, 0, 2) = c55;
  stiffnessTensorView(0, 2, 1, 0) = c56;
  stiffnessTensorView(0, 2, 1, 1) = c25;
  stiffnessTensorView(0, 2, 1, 2) = c45;
  stiffnessTensorView(0, 2, 2, 0) = c55;
  stiffnessTensorView(0, 2, 2, 1) = c45;
  stiffnessTensorView(0, 2, 2, 2) = c35;
  stiffnessTensorView(1, 0, 0, 0) = c16;
  stiffnessTensorView(1, 0, 0, 1) = c66;
  stiffnessTensorView(1, 0, 0, 2) = c56;
  stiffnessTensorView(1, 0, 1, 0) = c66;
  stiffnessTensorView(1, 0, 1, 1) = c26;
  stiffnessTensorView(1, 0, 1, 2) = c46;
  stiffnessTensorView(1, 0, 2, 0) = c56;
  stiffnessTensorView(1, 0, 2, 1) = c46;
  stiffnessTensorView(1, 0, 2, 2) = c36;
  stiffnessTensorView(1, 1, 0, 0) = c12;
  stiffnessTensorView(1, 1, 0, 1) = c26;
  stiffnessTensorView(1, 1, 0, 2) = c25;
  stiffnessTensorView(1, 1, 1, 0) = c26;
  stiffnessTensorView(1, 1, 1, 1) = c22;
  stiffnessTensorView(1, 1, 1, 2) = c24;
  stiffnessTensorView(1, 1, 2, 0) = c25;
  stiffnessTensorView(1, 1, 2, 1) = c24;
  stiffnessTensorView(1, 1, 2, 2) = c23;
  stiffnessTensorView(1, 2, 0, 0) = c14;
  stiffnessTensorView(1, 2, 0, 1) = c46;
  stiffnessTensorView(1, 2, 0, 2) = c45;
  stiffnessTensorView(1, 2, 1, 0) = c46;
  stiffnessTensorView(1, 2, 1, 1) = c24;
  stiffnessTensorView(1, 2, 1, 2) = c44;
  stiffnessTensorView(1, 2, 2, 0) = c45;
  stiffnessTensorView(1, 2, 2, 1) = c44;
  stiffnessTensorView(1, 2, 2, 2) = c34;
  stiffnessTensorView(2, 0, 0, 0) = c15;
  stiffnessTensorView(2, 0, 0, 1) = c56;
  stiffnessTensorView(2, 0, 0, 2) = c55;
  stiffnessTensorView(2, 0, 1, 0) = c56;
  stiffnessTensorView(2, 0, 1, 1) = c25;
  stiffnessTensorView(2, 0, 1, 2) = c45;
  stiffnessTensorView(2, 0, 2, 0) = c55;
  stiffnessTensorView(2, 0, 2, 1) = c45;
  stiffnessTensorView(2, 0, 2, 2) = c35;
  stiffnessTensorView(2, 1, 0, 0) = c14;
  stiffnessTensorView(2, 1, 0, 1) = c46;
  stiffnessTensorView(2, 1, 0, 2) = c45;
  stiffnessTensorView(2, 1, 1, 0) = c46;
  stiffnessTensorView(2, 1, 1, 1) = c24;
  stiffnessTensorView(2, 1, 1, 2) = c44;
  stiffnessTensorView(2, 1, 2, 0) = c45;
  stiffnessTensorView(2, 1, 2, 1) = c44;
  stiffnessTensorView(2, 1, 2, 2) = c34;
  stiffnessTensorView(2, 2, 0, 0) = c13;
  stiffnessTensorView(2, 2, 0, 1) = c36;
  stiffnessTensorView(2, 2, 0, 2) = c35;
  stiffnessTensorView(2, 2, 1, 0) = c36;
  stiffnessTensorView(2, 2, 1, 1) = c23;
  stiffnessTensorView(2, 2, 1, 2) = c34;
  stiffnessTensorView(2, 2, 2, 0) = c35;
  stiffnessTensorView(2, 2, 2, 1) = c34;
  stiffnessTensorView(2, 2, 2, 2) = c33;
}

// calculate maximal wave speed
double AnisotropicMaterial::getMaxWaveSpeed() const {
  // Wavespeeds for anisotropic materials depend on the direction of propagation.
  // An analytic solution for the maximal wave speed is hard to obtain.
  // Instead of solving an optimization problem we sample the velocitiy for
  // different directions and take the maximum.
  auto samplingDirections = seissol_general::init::samplingDirections::view::create(
      const_cast<double*>(seissol_general::init::samplingDirections::Values));

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> saes;

  double maxEv = 0;

  std::array<double, 81> fullTensor;
  getFullStiffnessTensor(fullTensor);
  seissol_general::kernel::computeChristoffel computeChristoffel;
  computeChristoffel.stiffnessTensor = fullTensor.data();

  for (unsigned j = 0; j < 200; ++j) {
    double n[3] = {samplingDirections(j, 0), samplingDirections(j, 1), samplingDirections(j, 2)};
    double m[9];
    computeChristoffel.direction = n;
    computeChristoffel.christoffel = m;
    computeChristoffel.execute();

    saes.compute(Eigen::Matrix<double, 3, 3>(m).cast<double>());
    auto eigenvalues = saes.eigenvalues();
    for (unsigned i = 0; i < 3; ++i) {
      maxEv = std::max(eigenvalues(i), maxEv);
    }
  }
  return std::sqrt(maxEv / rho);
}

// calculate P-wave speed based on averaged material parameters
double AnisotropicMaterial::getPWaveSpeed() const {
  const double muBar = (c44 + c55 + c66) / 3.0;
  const double lambdaBar = (c11 + c22 + c33) / 3.0 - 2.0 * muBar;
  return std::sqrt((lambdaBar + 2 * muBar) / rho);
}

// calculate S-wave speed based on averaged material parameters
double AnisotropicMaterial::getSWaveSpeed() const {
  const double muBar = (c44 + c55 + c66) / 3.0;
  return std::sqrt(muBar / rho);
}

MaterialType AnisotropicMaterial::getMaterialType() const { return MaterialType::Anisotropic; }
} // namespace seissol::model
