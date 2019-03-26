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

#include "Map.h"
#include <iostream>
#include <math.h>

#ifndef noproj
void printLastError()
{
  if (pj_errno != 0) {
    std::cout << "Proj.4 error: " << pj_strerrno(pj_errno) << std::endl;
  }
}


Map::Map(std::string const& targetCoordSystem)
{
  if (!(pj_lonlat = pj_init_plus("+proj=lonlat +datum=WGS84 +units=km +no_defs"))) { 
    printLastError();
  }
  if (!(pj_mesh = pj_init_plus(targetCoordSystem.c_str()))) { 
    printLastError();
  }
  
  if (pj_mesh->to_meter != 1.0 || pj_mesh->vto_meter != 1.0) {
    std::cout << "Warning: The MCS does not use meter as length unit. Note that area and sliprate will be saved in m^2 and m/s so prepare for trouble." << std::endl;
  }
}

Map::~Map()
{
  pj_free(pj_lonlat);
  pj_free(pj_mesh);  
}

void Map::map(double longitude, double latitude, double depth, double* x, double* y, double* z) const
{
  *x = longitude * DEG_TO_RAD;
  *y = latitude * DEG_TO_RAD;
  *z = -depth;
  pj_transform(pj_lonlat, pj_mesh, 1, 1, x, y, z);
  printLastError();
}

// transform from ned (srf coord system) to pj_mesh->axis
void Map::adjustAxes(double* x, double* y, double* z) const
{
  double p[3] = { *x, *y, *z }; // x = n, y = e, z = d
  double* v[3] = { x, y, z };
  
  for (unsigned i = 0; i < 3; ++i) {
    switch (pj_mesh->axis[i]) {
      case 'e':
        *v[i] = p[1];
        break;
      case 'w':
        *v[i] = -p[1];
        break;
      case 'n':
        *v[i] = p[0];
        break;
      case 's':
        *v[i] = -p[0];
        break;
      case 'd':
        *v[i] = p[2];
        break;
      case 'u':
        *v[i] = -p[2];
        break;
      default:
        std::cerr << "Weird axis definition: " << pj_mesh->axis[i] << std::endl;
        break;
    }
  }  
}
#else

Map::Map(std::string const& targetCoordSystem)
{}

Map::~Map()
{}
#define DEG_TO_RAD   .017453292519943296
#endif

void Map::toMCS(double strike, double dip, double rake, double u1, double u2, double u3, double* x, double* y, double* z) const
{
  strike *= DEG_TO_RAD;
  dip *= DEG_TO_RAD;
  rake *= DEG_TO_RAD;  
  *x = u1*(sin(rake)*sin(strike)*cos(dip) + cos(rake)*cos(strike)) + u2*(-sin(rake)*cos(strike) + sin(strike)*cos(dip)*cos(rake)) - u3*sin(dip)*sin(strike);
  *y = u1*(-sin(rake)*cos(dip)*cos(strike) + sin(strike)*cos(rake)) + u2*(-sin(rake)*sin(strike) - cos(dip)*cos(rake)*cos(strike)) + u3*sin(dip)*cos(strike);
  *z = -u1*sin(dip)*sin(rake) - u2*sin(dip)*cos(rake) - u3*cos(dip);
#ifndef noproj
  adjustAxes(x, y, z);
#endif
}
