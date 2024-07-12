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

#include "FilterParameters.h"

namespace seissol::initializer::parameters {
    FilterParameters readFilter(ParameterReader *baseReader){
        // TODO Remove duplication
  // Note: Filter parsing is currently duplicated in DR initialization.
  // If you adjust this code, also adjust in DynamicRupture/Parameters.h
    auto reader = baseReader->readSubNode("discretization");

  const auto type = reader->readWithDefaultStringEnum("filtertype", "identity", validFilters);

  // Compare this with Hesthaven Nodal DG: Alpha is set such that it reduces the highest mode to
  // epsilon
  const auto alpha = reader->readWithDefault("filteralpha", defaultFilterAlpha);

  const auto order = reader->readWithDefault("filterorder", defaultFilterOrder);
  const auto cutoff = reader->readWithDefault("filtercutoff", defaultFilterCutoff);

  if (type == FilterTypes::Exponential) {
    logInfo() << "Using a filter with order" << order << "cutoff" << cutoff
              << "and alpha" << alpha;
  }
  return FilterParameters{type, alpha, order, cutoff};
    }
} // namespace seissol::initializer::parameters
