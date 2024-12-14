/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 * @section LICENSE
 * Copyright (c) 2019 - 2020, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#ifndef MODEL_ANISOTROPIC_DATASTRUCTURES_H_
#define MODEL_ANISOTROPIC_DATASTRUCTURES_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <cstddef>
#include <string>

namespace seissol::model {
struct AnisotropicMaterial : Material {
  static constexpr std::size_t NumQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Anisotropic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "anisotropic";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

  double c11;
  double c12;
  double c13;
  double c14;
  double c15;
  double c16;
  double c22;
  double c23;
  double c24;
  double c25;
  double c26;
  double c33;
  double c34;
  double c35;
  double c36;
  double c44;
  double c45;
  double c46;
  double c55;
  double c56;
  double c66;

  double getLambdaBar() const override { return (c11 + c22 + c33) / 3.0 - 2.0 * getMuBar(); }

  double getMuBar() const override { return (c44 + c55 + c66) / 3.0; }

  AnisotropicMaterial() = default;

  explicit AnisotropicMaterial(ElasticMaterial m) {
    rho = m.rho;
    c11 = m.lambda + 2 * m.mu;
    c12 = m.lambda;
    c13 = m.lambda;
    c14 = 0;
    c15 = 0;
    c16 = 0;
    c22 = m.lambda + 2 * m.mu;
    c23 = m.lambda;
    c24 = 0;
    c25 = 0;
    c26 = 0;
    c33 = m.lambda + 2 * m.mu;
    c34 = 0;
    c35 = 0;
    c36 = 0;
    c44 = m.mu;
    c45 = 0;
    c46 = 0;
    c55 = m.mu;
    c56 = 0;
    c66 = m.mu;
  }

  AnisotropicMaterial(const double* materialValues, int numMaterialValues) {
    assert(numMaterialValues == 22);

    this->rho = materialValues[0];
    this->c11 = materialValues[1];
    this->c12 = materialValues[2];
    this->c13 = materialValues[3];
    this->c14 = materialValues[4];
    this->c15 = materialValues[5];
    this->c16 = materialValues[6];
    this->c22 = materialValues[7];
    this->c23 = materialValues[8];
    this->c24 = materialValues[9];
    this->c25 = materialValues[10];
    this->c26 = materialValues[11];
    this->c33 = materialValues[12];
    this->c34 = materialValues[13];
    this->c35 = materialValues[14];
    this->c36 = materialValues[15];
    this->c44 = materialValues[16];
    this->c45 = materialValues[17];
    this->c46 = materialValues[18];
    this->c55 = materialValues[19];
    this->c56 = materialValues[20];
    this->c66 = materialValues[21];
  }

  ~AnisotropicMaterial() override = default;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const final {
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
  // Wavespeeds for anisotropic materials depend on the direction of propagation.
  // An analytic solution for the maximal wave speed is hard to obtain.
  // Instead of solving an optimization problem we sample the velocitiy for
  // different directions and take the maximum.
  double getMaxWaveSpeed() const override {
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
    return sqrt(maxEv / rho);
  }

  // calculate P-wave speed based on averaged material parameters
  double getPWaveSpeed() const override {
    double muBar = (c44 + c55 + c66) / 3.0;
    double lambdaBar = (c11 + c22 + c33) / 3.0 - 2.0 * muBar;
    return std::sqrt((lambdaBar + 2 * muBar) / rho);
  }

  // calculate S-wave speed based on averaged material parameters
  double getSWaveSpeed() const override {
    double muBar = (c44 + c55 + c66) / 3.0;
    return std::sqrt(muBar / rho);
  }

  MaterialType getMaterialType() const override { return MaterialType::Anisotropic; }
};
} // namespace seissol::model

#endif
