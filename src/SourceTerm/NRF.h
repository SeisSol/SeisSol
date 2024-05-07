/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * Copyright (c) 2023, Intel Corporation
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
 * Point source computation.
 **/

#ifndef SOURCETERM_NRF_H_
#define SOURCETERM_NRF_H_

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <vector>

namespace seissol::sourceterm {
typedef struct Subfault_units {
  char* tinit;
  char* timestep;
  char* mu;
  char* area;
  char* tan1;
  char* tan2;
  char* normal;
} Subfault_units;

typedef struct Subfault {
  double tinit;
  double timestep;
  double mu;
  double area;
  Eigen::Vector3d tan1;
  Eigen::Vector3d tan2;
  Eigen::Vector3d normal;
} Subfault;

using Offsets = std::array<unsigned, 3u>;

struct NRF {
  std::vector<Eigen::Vector3d> centres;
  std::vector<Subfault> subfaults;
  std::vector<Offsets> sroffsets;
  std::array<std::vector<double>, 3u> sliprates;
  inline std::size_t size() { return centres.size(); }
};
} // namespace seissol::sourceterm

#endif
