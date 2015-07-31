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
 * "No-Refined" implementation (uses only first order data)
 **/

#ifndef REFINEMENT_NO_REFINEMENT_H_
#define REFINEMENT_NO_REFINEMENT_H_

#include "Refinement.h"
#include "Geometry/MeshReader.h"
#include "Numerical_aux/BasicFunction.h"

namespace refinement
{

//------------------------------------------------------------------------------

template<class T>
class NoRefinement : public Refinement<T>
{
private:
	/** Number of basis functions */
	const unsigned int m_nBasisFunctions;

public:
	NoRefinement(
            const MeshReader &meshReader,
            unsigned int order,
			unsigned int numVariables) :
            Refinement<T>(numVariables),
            m_nBasisFunctions(BasicFunction::getBasicFunctionsPerOrder(order))
	{
        this->m_vertices.resize(meshReader.getVertices().size()*3);
        this->m_cells.resize(meshReader.getElements().size()*4);

		const std::vector<Vertex> &vertices = meshReader.getVertices();
		for (unsigned int i = 0; i < vertices.size(); i++)
			this->setVertex(i, vertices[i].coords);

		const std::vector<Element> &elements = meshReader.getElements();
		for (unsigned int i = 0; i < elements.size(); i++) {
			this->setCell(i, elements[i].vertices);
		}
	}

	virtual void get(const double* idata,  const unsigned int* cellMap,
			int variable, double* odata) const
	{
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < this->nCells(); i++)
			odata[i] = idata[(cellMap[i] * this->m_numVariables + variable) * m_nBasisFunctions];
	}
};

//------------------------------------------------------------------------------

} // namespace

#endif // REFINEMENT_NO_REFINEMENT_H_
