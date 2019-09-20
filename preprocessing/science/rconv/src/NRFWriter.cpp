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

#include "NRFWriter.h"
#include <netcdf.h>
#include <iostream>
#include <cstring>
#include <limits>
#include <cmath>

using seissol::sourceterm::Subfault;
using seissol::sourceterm::Subfault_units;
using seissol::sourceterm::Offsets;

void check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
        (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
        exit(1);
    }
}

void writeNRF(char const* filename, std::vector<SRFPointSource> const& sources, Map const& map, bool normalizeOnset)
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
  if (normalizeOnset) {
    std::cout << "Minimal tinit: " << minTinit << " -> 0" << std::endl;
  } else {
    minTinit = 0.0;
  }
  
  Offsets* offsets = new Offsets[numSources+1];
  memset(offsets, 0, (numSources+1) * sizeof(Offsets));

  double** sliprates = new double*[3];
  for (unsigned sr = 0; sr < 3; ++sr) {
    sliprates[sr] = new double[numSamples[sr]];
    memset(sliprates[sr], 0, numSamples[sr] * sizeof(double));
  }

  glm::dvec3* centres = new glm::dvec3[numSources];
  Subfault* subfaults = new Subfault[numSources];
  for (unsigned i = 0; i < numSources; ++i) {
    SRFPointSource const& source = sources[i];
    glm::dvec3& centre = centres[i];
    Subfault& sf = subfaults[i];
    sf.tinit = source.tinit - minTinit;
    sf.timestep = source.dt;
    
#ifdef noproj
    centre.x = source.longitude;
    centre.y= source.latitude;
    centre.z=source.depth;
#else
    map.map(source.longitude, source.latitude, source.depth, &centre.x, &centre.y, &centre.z);
#endif

    map.toMCS(source.strike, source.dip, source.rake, 1.0, 0.0, 0.0, &sf.tan1.x, &sf.tan1.y, &sf.tan1.z);
    map.toMCS(source.strike, source.dip, source.rake, 0.0, 1.0, 0.0, &sf.tan2.x, &sf.tan2.y, &sf.tan2.z);
    map.toMCS(source.strike, source.dip, source.rake, 0.0, 0.0, 1.0, &sf.normal.x, &sf.normal.y, &sf.normal.z);
    
    // g / (cm s^2) -> Pa (= kg / (m s^2))
    sf.mu = source.shearModulus * 1.0e-1; 
    // cm^2 -> m^2
    sf.area = source.area * 1.0e-4;
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

  int  stat;  /* return status */
  int  ncid;  /* netCDF id */

  /* group ids */
  int nrf_grp;

  /* type ids */
  int Vector3_typ;
  int Subfault_units_typ;
  int Subfault_typ;

  /* dimension ids */
  int source_dim;
  int sroffset_dim;
  int direction_dim;
  int sample1_dim;
  int sample2_dim;
  int sample3_dim;

  /* dimension lengths */
  size_t source_len = numSources;
  size_t sroffset_len = numSources + 1;
  size_t direction_len = 3;
  size_t sample1_len = numSamples[0];
  size_t sample2_len = numSamples[1];
  size_t sample3_len = numSamples[2];

  /* variable ids */
  int centres_id;
  int subfaults_id;
  int sroffsets_id;
  int sliprates1_id;
  int sliprates2_id;
  int sliprates3_id;

  /* rank (number of dimensions) for each variable */
#   define RANK_centres 1
#   define RANK_subfaults 1
#   define RANK_sroffsets 2
#   define RANK_sliprates1 1
#   define RANK_sliprates2 1
#   define RANK_sliprates3 1

  /* variable shapes */
  int centres_dims[RANK_centres];
  int subfaults_dims[RANK_subfaults];
  int sroffsets_dims[RANK_sroffsets];
  int sliprates1_dims[RANK_sliprates1];
  int sliprates2_dims[RANK_sliprates2];
  int sliprates3_dims[RANK_sliprates3];

  /* enter define mode */
  stat = nc_create(filename, NC_CLOBBER|NC_NETCDF4, &ncid);
  check_err(stat,__LINE__,__FILE__);
  nrf_grp = ncid;

  stat = nc_def_compound(nrf_grp, sizeof(glm::dvec3), "Vector3", &Vector3_typ);    check_err(stat,__LINE__,__FILE__);
  {
  stat = nc_insert_compound(nrf_grp, Vector3_typ, "x", NC_COMPOUND_OFFSET(glm::dvec3,x), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Vector3_typ, "y", NC_COMPOUND_OFFSET(glm::dvec3,y), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Vector3_typ, "z", NC_COMPOUND_OFFSET(glm::dvec3,z), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  }

  stat = nc_def_compound(nrf_grp, sizeof(Subfault_units), "Subfault_units", &Subfault_units_typ);    check_err(stat,__LINE__,__FILE__);
  {
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "tinit", NC_COMPOUND_OFFSET(Subfault_units,tinit), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "timestep", NC_COMPOUND_OFFSET(Subfault_units,timestep), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "mu", NC_COMPOUND_OFFSET(Subfault_units,mu), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "area", NC_COMPOUND_OFFSET(Subfault_units,area), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "tan1", NC_COMPOUND_OFFSET(Subfault_units,tan1), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "tan2", NC_COMPOUND_OFFSET(Subfault_units,tan2), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_units_typ, "normal", NC_COMPOUND_OFFSET(Subfault_units,normal), NC_STRING);    check_err(stat,__LINE__,__FILE__);
  }

  stat = nc_def_compound(nrf_grp, sizeof(Subfault), "Subfault", &Subfault_typ);    check_err(stat,__LINE__,__FILE__);
  {
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "tinit", NC_COMPOUND_OFFSET(Subfault,tinit), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "timestep", NC_COMPOUND_OFFSET(Subfault,timestep), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "mu", NC_COMPOUND_OFFSET(Subfault,mu), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "area", NC_COMPOUND_OFFSET(Subfault,area), NC_DOUBLE);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "tan1", NC_COMPOUND_OFFSET(Subfault,tan1), Vector3_typ);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "tan2", NC_COMPOUND_OFFSET(Subfault,tan2), Vector3_typ);    check_err(stat,__LINE__,__FILE__);
  stat = nc_insert_compound(nrf_grp, Subfault_typ, "normal", NC_COMPOUND_OFFSET(Subfault,normal), Vector3_typ);    check_err(stat,__LINE__,__FILE__);
  }


  /* define dimensions */
  stat = nc_def_dim(nrf_grp, "source", source_len, &source_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(nrf_grp, "sroffset", sroffset_len, &sroffset_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(nrf_grp, "direction", direction_len, &direction_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(nrf_grp, "sample1", sample1_len, &sample1_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(nrf_grp, "sample2", sample2_len, &sample2_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(nrf_grp, "sample3", sample3_len, &sample3_dim);
  check_err(stat,__LINE__,__FILE__);

  /* define variables */

  centres_dims[0] = source_dim;
  stat = nc_def_var(nrf_grp, "centres", Vector3_typ, RANK_centres, centres_dims, &centres_id);
  check_err(stat,__LINE__,__FILE__);

  subfaults_dims[0] = source_dim;
  stat = nc_def_var(nrf_grp, "subfaults", Subfault_typ, RANK_subfaults, subfaults_dims, &subfaults_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_var_deflate(nrf_grp, subfaults_id, NC_NOSHUFFLE, 1, 1);
  check_err(stat,__LINE__,__FILE__);

  sroffsets_dims[0] = sroffset_dim;
  sroffsets_dims[1] = direction_dim;
  stat = nc_def_var(nrf_grp, "sroffsets", NC_UINT, RANK_sroffsets, sroffsets_dims, &sroffsets_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_var_deflate(nrf_grp, sroffsets_id, NC_NOSHUFFLE, 1, 1);
  check_err(stat,__LINE__,__FILE__);

  sliprates1_dims[0] = sample1_dim;
  stat = nc_def_var(nrf_grp, "sliprates1", NC_DOUBLE, RANK_sliprates1, sliprates1_dims, &sliprates1_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_var_deflate(nrf_grp, sliprates1_id, NC_NOSHUFFLE, 1, 1);
  check_err(stat,__LINE__,__FILE__);

  sliprates2_dims[0] = sample2_dim;
  stat = nc_def_var(nrf_grp, "sliprates2", NC_DOUBLE, RANK_sliprates2, sliprates2_dims, &sliprates2_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_var_deflate(nrf_grp, sliprates2_id, NC_NOSHUFFLE, 1, 1);
  check_err(stat,__LINE__,__FILE__);

  sliprates3_dims[0] = sample3_dim;
  stat = nc_def_var(nrf_grp, "sliprates3", NC_DOUBLE, RANK_sliprates3, sliprates3_dims, &sliprates3_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_var_deflate(nrf_grp, sliprates3_id, NC_NOSHUFFLE, 1, 1);
  check_err(stat,__LINE__,__FILE__);

  /* assign per-variable attributes */

  {
  stat = nc_put_att_text(nrf_grp, centres_id, "units", 1, "m");
  check_err(stat,__LINE__,__FILE__);
  }

  {
  static const Subfault_units units_att[1] = {{"s", "s", "pascal", "m^2", "m", "m", "m"}} ;
  stat = nc_put_att(nrf_grp, subfaults_id, "units", Subfault_units_typ, 1, units_att);
  check_err(stat,__LINE__,__FILE__);
  }

  {
  stat = nc_put_att_text(nrf_grp, sliprates1_id, "units", 3, "m/s");
  check_err(stat,__LINE__,__FILE__);
  }

  {
  stat = nc_put_att_text(nrf_grp, sliprates2_id, "units", 3, "m/s");
  check_err(stat,__LINE__,__FILE__);
  }

  {
  stat = nc_put_att_text(nrf_grp, sliprates3_id, "units", 3, "m/s");
  check_err(stat,__LINE__,__FILE__);
  }


  /* leave define mode */
  stat = nc_enddef (nrf_grp);
  check_err(stat,__LINE__,__FILE__);

  /* assign variable data */
  
  stat = nc_put_var(ncid, centres_id, &centres[0]);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_put_var(ncid, subfaults_id, &subfaults[0]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_put_var_uint(ncid, sroffsets_id, &offsets[0][0]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_put_var_double(ncid, sliprates1_id, sliprates[0]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_put_var_double(ncid, sliprates2_id, sliprates[1]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_put_var_double(ncid, sliprates3_id, sliprates[2]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_close(nrf_grp);
  check_err(stat,__LINE__,__FILE__);

  delete[] centres;
  delete[] subfaults;
  for (unsigned sr = 0; sr < 3; ++sr) {
    delete[] sliprates[sr];
  }
  delete[] sliprates;
  delete[] offsets;
}

