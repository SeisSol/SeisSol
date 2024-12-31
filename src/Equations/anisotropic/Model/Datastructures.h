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
#include <array>
#include <cstddef>
#include <string>

namespace seissol::model {
class AnisotropicLocalData;
class AnisotropicNeighborData;

struct AnisotropicMaterial : public Material {
  static constexpr std::size_t NumQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Anisotropic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevski;
  static inline const std::string Text = "anisotropic";
  static inline const std::array<std::string, NumQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

  using LocalSpecificData = AnisotropicLocalData;
  using NeighborSpecificData = AnisotropicNeighborData;

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

  [[nodiscard]] double getLambdaBar() const override;

  [[nodiscard]] double getMuBar() const override;

  AnisotropicMaterial();

  explicit AnisotropicMaterial(const ElasticMaterial& m);

  AnisotropicMaterial(const std::vector<double>& materialValues);

  ~AnisotropicMaterial() override;

  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override;

  // calculate maximal wave speed
  [[nodiscard]] double getMaxWaveSpeed() const override;

  // calculate P-wave speed based on averaged material parameters
  [[nodiscard]] double getPWaveSpeed() const override;

  // calculate S-wave speed based on averaged material parameters
  [[nodiscard]] double getSWaveSpeed() const override;

  [[nodiscard]] MaterialType getMaterialType() const override;
};
} // namespace seissol::model

#endif
