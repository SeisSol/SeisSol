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

#include "SRF.h"
#include <iostream>
#include <fstream>

void parseSamples(std::ifstream& in, std::vector<double>& samples, unsigned numSamples)
{
  samples.resize(numSamples);
  for (unsigned sampleNo = 0; sampleNo < numSamples; ++sampleNo) {
    in >> samples[sampleNo];
  }
}

void parsePOINTS_1_0(std::ifstream& in, SRFPointSource& ps)
{
  double dummy;
  unsigned nt[3];
  
  // Point source, line 1
  in >> ps.longitude >> ps.latitude >> ps.depth >> ps.strike >> ps.dip >> ps.area >> ps.tinit >> ps.dt;
  // Point source, line 2
  in >> ps.rake >> dummy >> nt[0] >> dummy >> nt[1] >> dummy >> nt[2];
  
  ps.shearModulus = 0.0; // = unknown
  
  for (unsigned i = 0; i < 3; ++i) {
    parseSamples(in, ps.slipRate[i], nt[i]);
  }
}

void parsePOINTS_2_0(std::ifstream& in, SRFPointSource& ps)
{
  double vs, den, dummy;
  unsigned nt[3];
  
  // Point source, line 1
  in >> ps.longitude >> ps.latitude >> ps.depth >> ps.strike >> ps.dip >> ps.area >> ps.tinit >> ps.dt >> vs >> den;
  // Point source, line 2
  in >> ps.rake >> dummy >> nt[0] >> dummy >> nt[1] >> dummy >> nt[2];
  
  ps.shearModulus = (vs > 0.0 && den > 0.0) ? vs * vs * den : 0.0;
  
  for (unsigned i = 0; i < 3; ++i) {
    parseSamples(in, ps.slipRate[i], nt[i]);
  }
}

std::vector<SRFPointSource> parseSRF(char const* filename)
{
  std::vector<SRFPointSource> srf;
  
  std::ifstream in(filename, std::ifstream::in);
  
  double version;
  in >> version;
  void (*parsePoints)(std::ifstream&, SRFPointSource&);
  
  if (version == 2.0) {
    parsePoints = parsePOINTS_2_0;
  } else if (version == 1.0) {
    parsePoints = parsePOINTS_1_0;
  } else {
    std::cerr << "Unsupported SRF version " << version << " in " << filename << "." << std::endl;
    return srf;
  }

  while (in.good()) {
    std::string block;
    unsigned num;
    in >> block;
    if (block.compare("POINTS") == 0) {
      in >> num;
      unsigned numSources = srf.size();
      unsigned newSize = numSources + num;
      srf.resize(newSize);
      for (unsigned point = numSources; point < newSize; ++point) {
        parsePoints(in, srf[point]);
      }
    } else if (block.compare("PLANE") == 0) {
      in >> num;
      // skip plane segments
      for (unsigned l = 0; l < 2*num && in.good(); ++ l) {
        std::string line;
        getline(in, line);
      }
    } else if (block.compare(0, 1, "#") == 0) { // ignore comments
      std::string line;
      getline(in, line); // skip line
    } else if (block.find_first_not_of(" \t\n\v\f\r") == std::string::npos) { // ignore whitespace line
      continue;
    } else {
      std::cerr << "Unsupported block type " << block << " in " << filename << "." << std::endl;
      return srf;
    }
  }
  
  in.close();
  
  // Erasing sources with zero samples
  unsigned numRemoved = 0;
  std::vector<SRFPointSource>::iterator source = srf.begin();
  while (source != srf.end()) {
    if (source->slipRate[0].size() != 0 || source->slipRate[1].size() != 0 || source->slipRate[2].size() != 0) {
      ++source;
    } else {
      ++numRemoved;
      source = srf.erase(source);
    }
  }
  if (numRemoved > 0) {
    std::cout << "Removed " << numRemoved << " point sources that did not contain samples." << std::endl;
  }

  return srf;
}

