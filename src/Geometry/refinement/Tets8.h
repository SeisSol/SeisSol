/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
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
 * Refines each tet into 8 subtets arcording to this paper:
 * http://www.ams.org/journals/mcom/1996-65-215/S0025-5718-96-00748-X/S0025-5718-96-00748-X.pdf
 **/

#ifndef REFINEMENT_TETS_8_H
#define REFINEMENT_TETS_8_H

#include <set>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include "Refinement.h"
#include "Geometry/MeshReader.h"

namespace refinement
{

class Tets8 : public Refinement
{
private:
	/** Number of variables */
	const unsigned int m_nVariables;

	/** Number of basis functions */
	const unsigned int m_nBasisFunctions;

public:
	Tets8(const MeshReader &meshReader, const unsigned int* cellMap,
			unsigned int numVariables, unsigned int numBasisFunctions)
		  : Refinement(meshReader.getElements().size() * 8),
			m_nVariables(numVariables), m_nBasisFunctions(numBasisFunctions)
	{
		/*
		int numTasks = 1;
#ifdef _OPENMP
		numTasks = omp_get_max_threads();
#endif // _OPENMP

		const std::set<Vertex> *vertexMaps = new std::set<Vertex>[numTasks];
		// we will generate 8 times as many cells (with 4 vertices per cell)
		m_cells = new unsigned int[numCells * 8 * 4];

		int taskId;
#ifdef _OPENMP
		#pragma omp parallel private(taskId)
		{
			taskId = omp_get_thread_num();
			#pragma omp for
#endif // _OPENMP
			// Create the new cells and a vertex map for each task
			for (unsigned int i = 0; i < numCells; i++) {

			}
#ifdef _OPENMP

			// Combine the vertexMaps to one (only required for OpenMP)
		}
#endif // _OPENMP

		delete vertexMaps;
		*/
	}
};

}

#endif // REFINEMENT_TETS_8_H
