/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 * Wave field Fortran interface
 */

#include "Parallel/MPI.h"

#include "utils/logger.h"

#include "SeisSol.h"
#include "WaveFieldWriter.h"
#include "Geometry/MeshReader.h"

/** Map from cells to dofs (really required for GENERATEDKERNELS */
static unsigned int* cellMap;

extern "C"
{

void wavefield_hdf_init(int rank, const char* outputPrefix,
		const double* dofs, const double* pstrain,
		int numVars, int order, int numBasisFuncs,
        int refinement,  int timestep, int* outputMask, double* outputRegionBounds)
{
	seissol::SeisSol::main.waveFieldWriter().enable();
	seissol::SeisSol::main.waveFieldWriter().setFilename(outputPrefix);

	// Create the map (really required for clustered lts)
	MeshReader& meshReader = seissol::SeisSol::main.meshReader();
	cellMap = new unsigned int[meshReader.getElements().size()];
#ifdef _OPENMP
	#pragma omp parallel for
#endif //_OPENMP
	for (unsigned int i = 0; i < meshReader.getElements().size(); i++)
		cellMap[i] = i;

	seissol::SeisSol::main.waveFieldWriter().init(numVars, order, numBasisFuncs,
			meshReader,	dofs, pstrain, 0L, cellMap, refinement, timestep, outputMask,
			outputRegionBounds, 0);

	// I/O is currently the last initialization that requires the mesh reader
	seissol::SeisSol::main.freeMeshReader();
}

void wavefield_hdf_close()
{
	seissol::SeisSol::main.waveFieldWriter().close();
	delete [] cellMap;
}

/**
 */
void wavefield_hdf_write_step(double time)
{
	seissol::SeisSol::main.waveFieldWriter().write(time);
}

}
