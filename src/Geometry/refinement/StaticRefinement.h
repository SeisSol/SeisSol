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
 * Refines each tet into 4 subtets arcording to this paper:
 * http://www.ams.org/journals/mcom/1996-65-215/S0025-5718-96-00748-X/S0025-5718-96-00748-X.pdf
 **/

#ifndef REFINEMENT_STATIC_REFINEMENT_H_
#define REFINEMENT_STATIC_REFINEMENT_H_

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include "Refinement.h"
#include "Geometry/MeshReader.h"
#include "Numerical_aux/BasicFunction.h"
#include "RefinerUtils.h"

namespace refinement
{

//------------------------------------------------------------------------------

template<class T>
class StaticRefinement : public Refinement<T>
{
private:
    std::vector<BasicFunction::SampledBasicFunctions<T> > m_BasicFunctions;

    std::vector<Vec3D<T> > m_CenterPoints;

    static const unsigned int kCoordsPerCell = 4 * 3;
    static const unsigned int kIndicesPerCell = 4;

    // Set and determined by the constructor
    const unsigned int kSubCellsPerCell;
    const unsigned int kOrder;

public:
    StaticRefinement(
            const MeshReader& meshReader,
            const TetrahedronRefiner<T>& tetRefiner,
            unsigned int order,
            unsigned int numVariables

            ) : Refinement<T>(numVariables),
        kSubCellsPerCell(tetRefiner.getDivisionCount()),
        kOrder(order)
    {
        // Generate cell centerpoints in the reference tetrahedron.
        std::vector<Tetrahedron<T> > subCells(kSubCellsPerCell);
        tetRefiner(Tetrahedron<T>::unitTetrahedron(), subCells.data());
        m_CenterPoints.resize(kSubCellsPerCell);
        for (unsigned int t = 0; t < kSubCellsPerCell; t++)
            m_CenterPoints[t] = subCells[t].center();

        // Generate sampled basicfunctions
        for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
            const Vec3D<T>& pnt = m_CenterPoints[i];
            m_BasicFunctions.push_back(
                    BasicFunction::SampledBasicFunctions<T>::sampleAt(
                        kOrder, pnt.x, pnt.y, pnt.z));
        }

        const std::size_t kInCellCount = meshReader.getElements().size();
        const std::size_t kSubCellCount = kInCellCount * kSubCellsPerCell;

        // Allocate memory for cells and vertices
        this->m_cells.resize(kSubCellCount * kIndicesPerCell);
        this->m_vertices.resize(kSubCellCount * kCoordsPerCell);

        const std::vector<Vertex>& kVertices = meshReader.getVertices();
        const std::vector<Element>& kElements = meshReader.getElements();

        // Start the acctual cell-refinement
#ifdef OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (unsigned int c = 0; c < kInCellCount; c++)
        {
            const unsigned int kCellIndexInVertexMap = c * kSubCellsPerCell * kCoordsPerCell;

            // Build a Terahedron containing the coordinates of the vertices.
            Tetrahedron<T> inTet = Tetrahedron<T>(
                        kVertices[kElements[c].vertices[0]].coords,
                        kVertices[kElements[c].vertices[1]].coords,
                        kVertices[kElements[c].vertices[2]].coords,
                        kVertices[kElements[c].vertices[3]].coords);
            // This works as long as the compiler does not add padding to the
            // Tetrahedron<double> class.
            Tetrahedron<T>* outTets = reinterpret_cast<Tetrahedron<T>*>(this->m_vertices.data()+kCellIndexInVertexMap);
            tetRefiner(inTet, outTets);
        }

        // Set the vertex indices (vid) of each tetrahedron.
        // Each tetrahedron has 4 indices and produces exactly 4 vertices.
        // Bijectivity is given and the order is well known.
#ifdef OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (unsigned int vid  = 0; vid < this->m_cells.size(); vid++) {
            this->m_cells[vid] = vid*3;
        }
    }

    virtual void get(const double* idata,  const unsigned int* cellMap,
            int variable, double* odata) const
    {
#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        // Iterate over original Cells c and subsampled cells i
        for (unsigned int i = 0; i < this->nCells()/kSubCellsPerCell; i++)
            // Iterate over refined subcells which are in a fixed order
            for (unsigned int sc = 0; sc < kSubCellsPerCell; sc++) {
                std::size_t offset = variable * m_NumAlignedBasicFunctions;
                odata[i*kSubCellsPerCell+sc] = m_BasicFunctions[sc].evalWithCoefs(
                        &idata[cellMap[i] * this->m_numVariables + offset]
                        );
            }
    }
};

//------------------------------------------------------------------------------

} // namespace

#endif // REFINEMENT_STATIC_REFINEMENT_H_
