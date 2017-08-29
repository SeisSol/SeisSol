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
 **/

#include "GMSH.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <unordered_map>

namespace
{
  static std::unordered_map<std::string, unsigned> boundary_codes {
    {"UNSPECIFIED", 0},
    {"AXIS", 1},
    {"CONJUGATE", 2},
    {"CONVECTION", 3},
    {"CYCLIC", 4},
    {"DEAD", 5},
    {"ELEMENT_SID", 6},
    {"ESPECIES", 7},
    {"EXHAUST_FAN", 8},
    {"FAN", 9},
    {"FREE_SURFACE", 10},
    {"GAP", 11},
    {"INFLOW", 12},
    {"INLET", 13},
    {"INLET_VENT", 14},
    {"INTAKE_FAN", 15},
    {"INTERFACE", 16},
    {"INTERIOR", 17},
    {"INTERNAL", 18},
    {"LIVE", 19},
    {"MASS_FLOW_INLET", 20},
    {"MELT", 21},
    {"MELT_INTERFACE", 22},
    {"MOVING_BOUNDARY", 23},
    {"NODE", 24},
    {"OUTFLOW", 25},
    {"OUTLET", 26},
    {"OUTLET_VENT", 27},
    {"PERIODIC", 28},
    {"PLOT", 29},
    {"POROUS", 30},
    {"POROUS_JUMP", 31},
    {"PRESSURE", 32},
    {"PRESSURE_FAR_FIELD", 33},
    {"PRESSURE_INFLOW", 34},
    {"PRESSURE_INLET", 35},
    {"PRESSURE_OUTFLOW", 36},
    {"PRESSURE_OUTLET", 37},
    {"RADIATION", 38},
    {"RADIATOR", 39},
    {"RECIRCULATION_INLET", 40},
    {"RECIRCULATION_OUTLET", 41},
    {"SLIP", 42},
    {"SREACTION", 43},
    {"SURFACE", 44},
    {"SYMMETRY", 45},
    {"TRACTION", 46},
    {"TRAJECTORY", 47},
    {"VELOCITY", 48},
    {"VELOCITY_INLET", 49},
    {"VENT", 50},
    {"WALL", 51},
    {"SPRING", 52},
  };
}

void error(std::string const& errMessage)
{
  std::cerr << "MSH parse error: " << errMessage << std::endl;
  exit(-1);
}

unsigned boundary_code(std::string const& regionName)
{
  auto name = regionName.substr(1, regionName.size()-2);
  std::transform(name.begin(), name.end(),name.begin(), ::toupper);
  return boundary_codes[name];
}

// Uses GAMBIT neu conventions. See GAMBIT NEUTRAL FILE FORMAT Appendix C.2.
template<> unsigned const Simplex<2>::Face2Nodes[][2] = {{0, 1}, {1, 2}, {2, 0}};
template<> unsigned const Simplex<3>::Face2Nodes[][3] = {{1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
