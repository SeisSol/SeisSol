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
#include "GAMBITWriter.h"

#include <utils/args.h>
#include <iostream>

template<unsigned DIM>
void convert(std::string const& in, std::string const& out) {
  std::cout << "Reading MSH..." << std::flush;
  GMSH<DIM> msh(in.c_str());
  std::cout << "finished." << std::endl;
  std::cout << "Writing NEU..." << std::flush;
  writeNEU(out.c_str(), msh);
  std::cout << "finished." << std::endl;
}

int main(int argc, char** argv)
{
  utils::Args args;
  args.addOption("dim", 'd', "Dimension (Defaults to 3, put 2 here for SeisSol 2D)", utils::Args::Required, false);
  args.addOption("input", 'i', "Input file (.msh)");
  args.addOption("output", 'o', "Output file (.neu)");

	if (args.parse(argc, argv) == utils::Args::Success) {
    unsigned dim = args.getArgument<unsigned>("dim", 3);
    std::string in = args.getArgument<std::string>("input");
    std::string out = args.getArgument<std::string>("output");
    
    std::cout << "Reading MSH..." << std::flush;
    switch (dim) {
      case 3:
        convert<3>(in, out);
        break;
      case 2:
        std::cerr << "Warning: The z-coordinate will be ignored. If your mesh is not parallel to the z=0 plane, you will be in trouble." << std::endl;
        convert<2>(in, out);
        break;
      default:
        std::cerr << "Error: Only dim=2 and dim=3 are supported." << std::endl;
        break;
    }    
	} else {
		return -1;
	}
  
  return 0;
}
