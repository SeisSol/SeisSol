/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
 */

#ifndef VARIABLE_SUBSAMPLER_H_
#define VARIABLE_SUBSAMPLER_H_

#include <cassert>
#include <algorithm>

#include <Eigen/Dense>

#include "Geometry/MeshReader.h"
#include "Numerical_aux/BasisFunction.h"
#include "RefinerUtils.h"

namespace seissol
{
namespace refinement
{

//------------------------------------------------------------------------------

template<class T>
class VariableSubsampler
{
private:
    std::vector<basisFunction::SampledBasisFunctions<T> > m_BasisFunctions;

    /** The original number of cells (without refinement) */
    const unsigned int m_numCells;

    const unsigned int kSubCellsPerCell;
    const unsigned int kNumVariables;
    const unsigned int kNumAlignedDOF;


    std::size_t getInVarOffset(unsigned int cell, unsigned int variable,
    		const unsigned int* cellMap) const
    {
        return (cellMap[cell]*kNumVariables + variable) * kNumAlignedDOF;
    }

    std::size_t getOutVarOffset(unsigned cell, unsigned int subcell) const
    {
        return kSubCellsPerCell * cell + subcell;
    }

public:
    VariableSubsampler(
    		unsigned int numCells,
            const TetrahedronRefiner<T>& tetRefiner,
            unsigned int order,
            unsigned int numVariables,
            unsigned int numAlignedDOF
            );

    void get(const real* inData, const unsigned int* cellMap,
            int variable, real* outData) const;
};

//------------------------------------------------------------------------------

template<typename T>
VariableSubsampler<T>::VariableSubsampler(
		unsigned int numCells,
        const TetrahedronRefiner<T>& tetRefiner,
        unsigned int order,
        unsigned int numVariables,
        unsigned int numAlignedDOF)
		: m_numCells(numCells),
		  kSubCellsPerCell(tetRefiner.getDivisionCount()),
    kNumVariables(numVariables), kNumAlignedDOF(numAlignedDOF)
{
    // Generate cell centerpoints in the reference or unit tetrahedron.
	Tetrahedron<T>* subCells = new Tetrahedron<T>[kSubCellsPerCell];
        Eigen::Matrix<T, 3, 1>* additionalVertices = new Eigen::Matrix<T, 3, 1>[tetRefiner.additionalVerticesPerCell()];

    tetRefiner.refine(Tetrahedron<T>::unitTetrahedron(), 0,
    		subCells, additionalVertices);

    // Generate sampled basicfunctions
    for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
        const Eigen::Matrix<T, 3, 1> pnt = subCells[i].center();
        m_BasisFunctions.push_back(
                basisFunction::SampledBasisFunctions<T>(
                    order, pnt(0), pnt(1), pnt(2)));
    }

    delete [] subCells;
    delete [] additionalVertices;
}

//------------------------------------------------------------------------------

template<typename T>
void VariableSubsampler<T>::get(const real* inData,  const unsigned int* cellMap,
        int variable, real* outData) const
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    // Iterate over original Cells
    for (unsigned int c = 0; c < m_numCells; ++c) {
        for (unsigned int sc = 0; sc < kSubCellsPerCell; ++sc) {
            outData[getOutVarOffset(c, sc)] =
            		m_BasisFunctions[sc].evalWithCoeffs(&inData[getInVarOffset(c, variable, cellMap)]);
        }
    }
}

//------------------------------------------------------------------------------

} // namespace
}

#endif // VARIABLE_SUBSAMPLER_H_
