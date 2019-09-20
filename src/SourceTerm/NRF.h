/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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

#include <cstddef>
#include <glm/vec3.hpp>

namespace seissol {
  namespace sourceterm {
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
        glm::dvec3 tan1;
        glm::dvec3 tan2;
        glm::dvec3 normal;
    } Subfault;
    
    typedef unsigned Offsets[3];

    struct NRF {
      glm::dvec3* centres;
      Subfault* subfaults;
      Offsets* sroffsets;
      double* sliprates[3];
      size_t source;
      NRF() : centres(NULL), subfaults(NULL), sroffsets(NULL), source(0) {
        sliprates[0] = NULL;
        sliprates[1] = NULL;
        sliprates[2] = NULL;
      }
      ~NRF() {
        delete[] centres;
        delete[] subfaults;
        delete[] sroffsets;
        source = 0;
        delete[] sliprates[0];
        delete[] sliprates[1];
        delete[] sliprates[2];
      }
    };
  }
}

#endif
