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

#ifndef MODEL_COMMONDATASTRUCTURES_HPP_
#define MODEL_COMMONDATASTRUCTURES_HPP_

#include "Kernels/precision.hpp"
#include <array>
#include <string>

namespace seissol::model {
enum class MaterialType { solid, acoustic, elastic, viscoelastic, anisotropic, poroelastic };

// the local solvers. CK is the default for elastic, acoustic etc.
// viscoelastic uses CauchyKovalevskiAnelastic (maybe all other materials may be extended to use
// that one as well) poroelastic uses SpaceTimePredictorPoroelastic (someone may generalize that
// one, but so long I(David) had decided to put poroelastic in its name) the solver Unknown is a
// dummy to let all other implementations fail
enum class LocalSolver {
  Unknown,
  CauchyKovalevski,
  CauchyKovalevskiAnelastic,
  SpaceTimePredictorPoroelastic
};

struct Material {
  static constexpr std::size_t NumberOfQuantities = 0;        // ?
  static constexpr std::size_t NumberPerMechanism = 0;        // ?
  static constexpr std::size_t Mechanisms = 0;                // ?
  static constexpr MaterialType Type = MaterialType::solid;   // ?
  static constexpr LocalSolver Solver = LocalSolver::Unknown; // ?
  static inline const std::string Text = "material";
  static inline const std::array<std::string, NumberOfQuantities> Quantities = {};

  double rho;
  virtual double getMaxWaveSpeed() const = 0;
  virtual double getPWaveSpeed() const = 0;
  virtual double getSWaveSpeed() const = 0;
  virtual double getMu() const = 0;
  virtual void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const = 0;
  virtual MaterialType getMaterialType() const = 0;
};

struct Plasticity {
  static const inline std::string Text = "plasticity";
  double bulkFriction;
  double plastCo;
  double s_xx;
  double s_yy;
  double s_zz;
  double s_xy;
  double s_yz;
  double s_xz;
};

struct IsotropicWaveSpeeds {
  double density;
  double pWaveVelocity;
  double sWaveVelocity;
};
} // namespace seissol::model

#endif
