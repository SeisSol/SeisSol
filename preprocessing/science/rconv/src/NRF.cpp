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

#include "NRF.h"
#include <netcdf.h>
#include <iostream>
#include <cstring>
#include <limits>
#include <cmath>

void printError(int code)
{
  if (code != 0) {
    std::cout << "Netcdf error: " << nc_strerror(code) << std::endl;
  }
}

void writeNRF(char const* filename, std::vector<SRFPointSource> const& sources, Map const& map)
{  
  unsigned numSources = sources.size();  
  double minTinit = std::numeric_limits<double>::max();
  unsigned numSamples[3] = { 0, 0, 0 };
  
  for (std::vector<SRFPointSource>::const_iterator source = sources.begin(); source != sources.end(); ++source) {
    minTinit = std::min(minTinit, source->tinit);
    for (unsigned sr = 0; sr < 3; ++sr) {
      numSamples[sr] += source->slipRate[sr].size();
    }
  }  
  std::cout << "Minimal tinit: " << minTinit << " -> 0" << std::endl;
  
  Offsets* offsets = new Offsets[numSources+1];
  memset(offsets, 0, (numSources+1) * sizeof(Offsets));

  double** sliprates = new double*[3];
  for (unsigned sr = 0; sr < 3; ++sr) {
    sliprates[sr] = new double[numSamples[sr]];
    memset(sliprates[sr], 0, numSamples[sr] * sizeof(double));
  }

  double* basicparams = new double[numSources * NumberOfBasicParams];
  for (unsigned i = 0; i < numSources; ++i) {
    SRFPointSource const& source = sources[i];
    double* param = &basicparams[i * NumberOfBasicParams];
    param[Tinit] = source.tinit - minTinit;
    param[Timestep] = source.dt;
    
    map.map(source.longitude, source.latitude, source.depth, &param[CentreX], &param[CentreY], &param[CentreZ]);
    map.toMCS(source.strike, source.dip, source.rake, 1.0, 0.0, 0.0, &param[SlipX], &param[SlipY], &param[SlipZ]);
    map.toMCS(source.strike, source.dip, source.rake, 0.0, 1.0, 0.0, &param[TanSlipX], &param[TanSlipY], &param[TanSlipZ]);
    map.toMCS(source.strike, source.dip, source.rake, 0.0, 0.0, 1.0, &param[NormalX], &param[NormalY], &param[NormalZ]);    
    
    // cm^2 -> m^2
    param[Area] = source.area * 1.0e-4;
    for (unsigned sr = 0; sr < 3; ++sr) {
      unsigned offset = offsets[i][sr];
      std::copy(source.slipRate[sr].begin(), source.slipRate[sr].end(), &sliprates[sr][ offset ]);
      offsets[i+1][sr] = offset + source.slipRate[sr].size();
    }
  }
  
  // convert cm / s -> m / s
  for (unsigned sr = 0; sr < 3; ++sr) {
    for (unsigned sample = 0; sample < numSamples[sr]; ++sample) {
      sliprates[sr][sample] *= 1.0e-2;
    }
  }
  
  int err, ncid, basicparamsid, offsetsid, sourcedimid, offsetdimid, dimdimid, basicparamdimid;
  int slipratesid[3];
  int sampledimid[3];
  int dimids[2];
  printError( nc_create(filename, NC_NETCDF4, &ncid) );
  
  printError( nc_def_dim(ncid, "source", numSources, &sourcedimid) );
  printError( nc_def_dim(ncid, "offset", numSources+1, &offsetdimid) );
  printError( nc_def_dim(ncid, "dim", 3, &dimdimid) );
  printError( nc_def_dim(ncid, "basicparam", NumberOfBasicParams, &basicparamdimid) );
  printError( nc_def_dim(ncid, "sample1", numSamples[0], &sampledimid[0]) );
  printError( nc_def_dim(ncid, "sample2", numSamples[1], &sampledimid[1]) );
  printError( nc_def_dim(ncid, "sample3", numSamples[2], &sampledimid[2]) );
  
  dimids[0] = sourcedimid;
  dimids[1] = basicparamdimid;
  printError( nc_def_var(ncid, "basicparams", NC_DOUBLE, 2, dimids, &basicparamsid) );
  
  dimids[0] = offsetdimid;
  dimids[1] = dimdimid;
  printError( nc_def_var(ncid, "offsets", NC_UINT, 2, dimids, &offsetsid) );

  printError( nc_def_var(ncid, "sliprates1", NC_DOUBLE, 1, &sampledimid[0], &slipratesid[0]) );
  printError( nc_def_var(ncid, "sliprates2", NC_DOUBLE, 1, &sampledimid[1], &slipratesid[1]) );
  printError( nc_def_var(ncid, "sliprates3", NC_DOUBLE, 1, &sampledimid[2], &slipratesid[2]) );
  
  // Turn on compression
  printError( nc_def_var_deflate(ncid, basicparamsid, 0, 1, 1) );
  printError( nc_def_var_deflate(ncid, offsetsid, 0, 1, 1) );
  printError( nc_def_var_deflate(ncid, slipratesid[0], 0, 1, 1) );
  printError( nc_def_var_deflate(ncid, slipratesid[1], 0, 1, 1) );
  printError( nc_def_var_deflate(ncid, slipratesid[2], 0, 1, 1) );
  
  printError( nc_enddef(ncid) );
  
  printError( nc_put_var_double(ncid, basicparamsid, &basicparams[0]) );
  printError( nc_put_var_uint(ncid, offsetsid, &offsets[0][0]) );
  for (unsigned sr = 0; sr < 3; ++sr) {
    printError( nc_put_var_double(ncid, slipratesid[sr], sliprates[sr]) );
  }
  
  printError( nc_close(ncid) );
  
  delete[] basicparams;
  for (unsigned sr = 0; sr < 3; ++sr) {
    delete[] sliprates[sr];
  }
  delete[] sliprates;
  delete[] offsets;
}

