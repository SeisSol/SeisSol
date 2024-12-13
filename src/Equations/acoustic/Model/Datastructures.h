/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 * @author Jinwen Pan (jinwen.pan AT tum.de)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2024, SeisSol Group
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
 * This header file defines the AcousticMaterial struct, which represents the acoustic material
 * properties used in simulations. It inherits from the Material class and includes:
 * 
 * 1. Static constants for the number of quantities, material type, and solver type.
 * 
 * 2. The material's properties such as density (rho) and lambda, along with methods for 
 *    calculating wave speeds (P-wave and S-wave) and the stiffness tensor.
 * 
 * 3. A constructor to initialize the material from provided values, and a method to return 
 *    the maximum wave speed for the acoustic material.
 **/

#ifndef MODEL_ACOUSTIC_DATASTRUCTURES_H_
#define MODEL_ACOUSTIC_DATASTRUCTURES_H_

#include "Model/CommonDatastructures.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <string>

namespace seissol::model {
struct AcousticMaterial : Material {
  static constexpr std::size_t NumQuantities = 4;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Acoustic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "acoustic";
  // The stress-velocity formulation of the elastic model is reused. 
  // By definition, the normal stress and pressure are negatives of each other.
  static inline const std::array<std::string, NumQuantities> Quantities = {
      "-p", "v1", "v2", "v3"};

  double lambda;

  double getLambdaBar() const override { return lambda; }

  double getMuBar() const override { return 0.0; }

  AcousticMaterial() = default;
  AcousticMaterial(const double* materialValues, int numMaterialValues) {
    assert(numMaterialValues == 2);

    this->rho = materialValues[0];
    this->lambda = materialValues[1];
  }

  ~AcousticMaterial() override = default;

  // The stiffness tensor of the elastic model is reused. 
  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {

    auto stiffnessTensorView =
        seissol_general::init::stiffnessTensor::view::create(fullTensor.data());
    stiffnessTensorView.setZero();
    stiffnessTensorView(0, 0, 0, 0) = lambda;
    stiffnessTensorView(0, 0, 1, 1) = lambda;
    stiffnessTensorView(0, 0, 2, 2) = lambda;
    stiffnessTensorView(1, 1, 0, 0) = lambda;
    stiffnessTensorView(1, 1, 1, 1) = lambda;
    stiffnessTensorView(1, 1, 2, 2) = lambda;
    stiffnessTensorView(2, 2, 0, 0) = lambda;
    stiffnessTensorView(2, 2, 1, 1) = lambda;
    stiffnessTensorView(2, 2, 2, 2) = lambda;
  }

  double getMaxWaveSpeed() const override { return getPWaveSpeed(); }

  double getPWaveSpeed() const override { return std::sqrt(lambda / rho); }

  double getSWaveSpeed() const override { return 0.0; }

  MaterialType getMaterialType() const override { return Type; }
};
} // namespace seissol::model

#endif
