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
#include "NRFWriter.h"
#include "XMFWriter.h"

#include <utils/args.h>
#include <iostream>


int main(int argc, char** argv)
{
  utils::Args args;
  args.addOption("input", 'i', "Input file (.srf)");
	args.addOption("mcs", 'm', "Proj.4 string that describes the mesh coordinate system (e.g. \"+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs\").", utils::Args::Required, false);
  args.addOption("output", 'o', "Output file (.nrf)", utils::Args::Required, false);
  args.addOption("normalize-onset", 'n', "Subtract the minimum onset time from all onsets.", utils::Args::No, false);
	args.addOption("vcs", 'v', "Proj.4 string that describes the coordinate system for visualisation (defaults to geocentric if mcs not given, i.e. \"+proj=geocent +datum=WGS84 +units=m +no_def\").", utils::Args::Required, false);
  args.addOption("xdmf", 'x', "Output for visualisation (.xmf)", utils::Args::Required, false);
  
  args.setCustomHelpMessage("\nWith rconv you may either convert a SRF file to a NRF file, which you can use as input in SeisSol.\n"
                            "In this case, give the options -i, -m, -o, and optionally -n.\n\n"
                            "You may also write a file which may be loaded in Paraview for visualisation of the SRF file.\n"
                            "In this case, give the options -i, -x, and optionally -v.\n\n"
                            "You may write both files simultaneously by giving all options.\n");

	if (args.parse(argc, argv) == utils::Args::Success) {
    std::string in = args.getArgument<std::string>("input");
    std::string mcs = args.getArgument<std::string>("mcs", "");
    std::string out = args.getArgument<std::string>("output", "");
    bool normalizeOnset = args.isSet("normalize-onset");
    std::string vcs = args.getArgument<std::string>("vcs", "");
    std::string xdmf = args.getArgument<std::string>("xdmf", "");
    
    if (mcs.empty() && vcs.empty()) {
      vcs = "+proj=geocent +datum=WGS84 +units=m +no_defs";
    } else if (vcs.empty()) {
      vcs = mcs;
    }
    
    std::cout << "Reading SRF..." << std::flush;
    std::vector<SRFPointSource> srf = parseSRF(in.c_str());
    std::cout << "finished." << std::endl;

    if (mcs.empty() != out.empty()) {
      std::cerr << "Error: -o and -m may only be given simultaneously." << std::endl;
      return -1;
    } else if (!mcs.empty()) {
      Map map(mcs);    
      std::cout << "Writing NRF..." << std::flush;
      writeNRF(out.c_str(), srf, map, normalizeOnset);
      std::cout << "finished." << std::endl;
    }
    
    if (!xdmf.empty()) {
      Map mapGeocent(vcs);
      std::cout << "Writing XDMF..." << std::flush;
      writeXMF(xdmf.c_str(), srf, mapGeocent);
      std::cout << "finished." << std::endl;
    }
	} else {
		return -1;
	}
  
  return 0;
}
