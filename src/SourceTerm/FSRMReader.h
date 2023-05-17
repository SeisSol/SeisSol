/**
 * @file
 * This file is part of SeisSol.
 *
 * @author David Schneller
 *
 * @section LICENSE
 * Copyright (c) 2023, SeisSol Group
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
 * Point source datastructure.
 **/

#ifndef READER_FSRMREADER_H_
#define READER_FSRMREADER_H_

#include <vector>
#include <array>
#include <cstddef>
#include <string>
#include <Eigen/Dense>

#include "Initializer/BasicTypedefs.hpp"

namespace seissol {
namespace sourceterm {

// TODO: when refactoring, replace raw array types
struct FSRMSource {
  real momentTensor[3][3];
  real solidVelocityComponent[3];
  real pressureComponent;
  real fluidVelocityComponent[3];
  size_t numberOfSources;
  std::vector<Eigen::Vector3d> centers;
  std::vector<real> strikes;
  std::vector<real> dips;
  std::vector<real> rakes;
  std::vector<real> onsets;
  std::vector<real> areas;
  real timestep;
  size_t numberOfSamples;
  std::vector<std::vector<real>> timeHistories;

  void read(const std::string& filename);
};
} // namespace sourceterm
} // namespace seissol

#endif
