/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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

#ifndef MODEL_ELASTIC_DATASTRUCTURES_H_
#define MODEL_ELASTIC_DATASTRUCTURES_H_

#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <string>

namespace seissol::model {
struct ElasticMaterial : Material {
  static constexpr std::size_t NumQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Elastic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "elastic";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

  double lambda;
  double mu;

  double getLambdaBar() const override { return lambda; }

  double getMuBar() const override { return mu; }

  ElasticMaterial() = default;
  ElasticMaterial(const double* materialValues, int numMaterialValues) {
    assert(numMaterialValues == 3);

    this->rho = materialValues[0];
    this->mu = materialValues[1];
    this->lambda = materialValues[2];
  }

  ~ElasticMaterial() override = default;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {

    auto stiffnessTensorView =
        seissol_general::init::stiffnessTensor::view::create(fullTensor.data());
    stiffnessTensorView.setZero();
    stiffnessTensorView(0, 0, 0, 0) = lambda + 2 * mu;
    stiffnessTensorView(0, 0, 1, 1) = lambda;
    stiffnessTensorView(0, 0, 2, 2) = lambda;
    stiffnessTensorView(0, 1, 0, 1) = mu;
    stiffnessTensorView(0, 1, 1, 0) = mu;
    stiffnessTensorView(0, 2, 0, 2) = mu;
    stiffnessTensorView(0, 2, 2, 0) = mu;
    stiffnessTensorView(1, 0, 0, 1) = mu;
    stiffnessTensorView(1, 0, 1, 0) = mu;
    stiffnessTensorView(1, 1, 0, 0) = lambda;
    stiffnessTensorView(1, 1, 1, 1) = lambda + 2 * mu;
    stiffnessTensorView(1, 1, 2, 2) = lambda;
    stiffnessTensorView(1, 2, 1, 2) = mu;
    stiffnessTensorView(1, 2, 2, 1) = mu;
    stiffnessTensorView(2, 0, 0, 2) = mu;
    stiffnessTensorView(2, 0, 2, 0) = mu;
    stiffnessTensorView(2, 1, 2, 1) = mu;
    stiffnessTensorView(2, 2, 0, 0) = lambda;
    stiffnessTensorView(2, 2, 1, 1) = lambda;
    stiffnessTensorView(2, 2, 2, 2) = lambda + 2 * mu;
  }

  double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  double getPWaveSpeed() const override { return std::sqrt((lambda + 2 * mu) / rho); }

  double getSWaveSpeed() const override { return std::sqrt(mu / rho); }

  MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif
