/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Vikas Kurapati (vikas.kurapati AT cit.tum.de, https://www.cs.cit.tum.de/sccs/personen/vikas-kurapati/)
 *
 * @section LICENSE
 * Copyright (c) 2024, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Filter Parameters Initialization.
 **/

#ifndef SEISSOL_FILTER_PARAMETERS_H
#define SEISSOL_FILTER_PARAMETERS_H

#include <Initializer/InputParameters.hpp>
#include <Initializer/Parameters/ParameterReader.h>
#include <Kernels/precision.hpp>

namespace seissol::initializer::parameters{
enum class FilterTypes { Identity, Exponential };

struct FilterParameters {
  FilterTypes type;
  real alpha;
  unsigned int order = 32;
  unsigned int cutoff = 0;
};

const static auto validFilters = std::unordered_map<std::string, FilterTypes>{
    {"identity", FilterTypes::Identity},
    {"exponential", FilterTypes::Exponential},
};

// Compare this with Hesthaven Nodal DG: Alpha is set such that it reduces the highest mode to
// epsilon
const static real defaultFilterAlpha = -std::log(std::numeric_limits<real>::epsilon());
const static unsigned int defaultFilterOrder = 32;
const static unsigned int defaultFilterCutoff = 0;

FilterParameters readFilter(ParameterReader* baseReader);

} // namespace seissol::initializer::parameters
#endif // SEISSOL_FILTER_PARAMETERS_H