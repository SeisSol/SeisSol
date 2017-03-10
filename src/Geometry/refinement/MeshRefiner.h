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

#ifndef MESH_REFINER_H_
#define MESH_REFINER_H_

#include <cstring>

#include "Geometry/MeshReader.h"
#include "RefinerUtils.h"

namespace seissol
{
namespace refinement
{

//------------------------------------------------------------------------------

template<typename T>
class MeshRefiner
{
private:
    // m_cells contains the indices of the cells
    unsigned int* m_cells;
    T* m_vertices;

    size_t m_numSubCells;
    size_t m_numVertices;

    static const unsigned int kIndicesPerCell = 4;

    const unsigned int kSubCellsPerCell;

public:
    MeshRefiner(const MeshReader& meshReader,
            const TetrahedronRefiner<T>& tetRefiner);

    MeshRefiner(const std::vector<const Element *>& subElements,
            const std::vector<const Vertex *>& subVertices,
            const std::map<int, int>& oldToNewVertexMap,
            const TetrahedronRefiner<T>& tetRefiner);

    ~MeshRefiner();

    const unsigned int* getCellData() const;
    const T* getVertexData() const;

    std::size_t getNumCells() const;
    std::size_t getNumVertices() const;
};

//------------------------------------------------------------------------------

template<typename T>
MeshRefiner<T>::MeshRefiner(
        	const MeshReader& meshReader,
        	const TetrahedronRefiner<T>& tetRefiner)
		: kSubCellsPerCell(tetRefiner.getDivisionCount())

{
    using std::size_t;

    const size_t kInVertexCount = meshReader.getVertices().size();
    const size_t kInCellCount = meshReader.getElements().size();
    m_numSubCells = kInCellCount * kSubCellsPerCell;

    const unsigned int additionalVertices = tetRefiner.additionalVerticesPerCell();
    m_numVertices = kInVertexCount + kInCellCount * additionalVertices;

    m_cells = new unsigned int[m_numSubCells * kIndicesPerCell];
    m_vertices = new T[m_numVertices * 3];

    const std::vector<Vertex>& kVertices = meshReader.getVertices();
    const std::vector<Element>& kElements = meshReader.getElements();

    // Copy original vertices
#ifdef _OPENMP
	#pragma omp parallel for
#endif // _OPENMP
    for (unsigned int i = 0; i < kInVertexCount; i++) {
    	memcpy(&m_vertices[i*3], kVertices[i].coords, sizeof(double)*3);
    }

    // The pointer to the new vertices
    T* newVertices = &m_vertices[kInVertexCount*3];

    // Start the actual cell-refinement
#ifdef _OPENMP
    #pragma omp parallel
    {
#endif // _OPENMPI
    	glm::tvec3<T>* newVerticesTmp = new glm::tvec3<T>[additionalVertices];
    	Tetrahedron<T>* newTetsTmp = new Tetrahedron<T>[kSubCellsPerCell];

#ifdef _OPENMP
        #pragma omp for schedule(static) nowait
#endif // _OPENMP
        for (size_t c = 0; c < kInCellCount; ++c)
        {
            // Build a Terahedron containing the coordinates of the vertices.
            Tetrahedron<T> inTet = Tetrahedron<T>(
                    kVertices[kElements[c].vertices[0]].coords,
                    kVertices[kElements[c].vertices[1]].coords,
                    kVertices[kElements[c].vertices[2]].coords,
                    kVertices[kElements[c].vertices[3]].coords,
					kElements[c].vertices[0],
					kElements[c].vertices[1],
					kElements[c].vertices[2],
					kElements[c].vertices[3]);

            // Generate the tets
            tetRefiner.refine(inTet,
            		kInVertexCount + c*additionalVertices,
					newTetsTmp, newVerticesTmp);

            // Copy new vertices
            for (unsigned int i = 0; i < additionalVertices; i++) {
            	memcpy(&newVertices[(c*additionalVertices + i) * 3],
            			glm::value_ptr(newVerticesTmp[i]), sizeof(T)*3);
            }

            // Copy tets
            for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
            	m_cells[(c*kSubCellsPerCell + i) * 4] = newTetsTmp[i].i;
            	m_cells[(c*kSubCellsPerCell + i) * 4 + 1] = newTetsTmp[i].j;
            	m_cells[(c*kSubCellsPerCell + i) * 4 + 2] = newTetsTmp[i].k;
            	m_cells[(c*kSubCellsPerCell + i) * 4 + 3] = newTetsTmp[i].l;
            }
        }

        delete [] newVerticesTmp;
        delete [] newTetsTmp;
#ifdef _OPENMP
    }
#endif
};

template<typename T>
MeshRefiner<T>::MeshRefiner(
            const std::vector<const Element *>& subElements,
            const std::vector<const Vertex *>& subVertices,
            const std::map<int, int>& oldToNewVertexMap,
        	const TetrahedronRefiner<T>& tetRefiner)
		: kSubCellsPerCell(tetRefiner.getDivisionCount())

{
    using std::size_t;

    const size_t kInVertexCount = subVertices.size();
    const size_t kInCellCount = subElements.size();
    m_numSubCells = kInCellCount * kSubCellsPerCell;

    const unsigned int additionalVertices = tetRefiner.additionalVerticesPerCell();
    m_numVertices = kInVertexCount + kInCellCount * additionalVertices;

    m_cells = new unsigned int[m_numSubCells * kIndicesPerCell];
    m_vertices = new T[m_numVertices * 3];

    const std::vector<const Vertex*>& kVertices = subVertices;
    const std::vector<const Element*>& kElements = subElements;

    // Copy original vertices
#ifdef _OPENMP
	#pragma omp parallel for
#endif // _OPENMP
    for (unsigned int i = 0; i < kInVertexCount; i++) {
    	memcpy(&m_vertices[i*3], kVertices[i]->coords, sizeof(double)*3);
    }

    // The pointer to the new vertices
    T* newVertices = &m_vertices[kInVertexCount*3];

    // Start the actual cell-refinement
#ifdef _OPENMP
    #pragma omp parallel shared(oldToNewVertexMap)
    {
#endif // _OPENMPI
    	glm::tvec3<T>* newVerticesTmp = new glm::tvec3<T>[additionalVertices];
    	Tetrahedron<T>* newTetsTmp = new Tetrahedron<T>[kSubCellsPerCell];

#ifdef _OPENMP
        #pragma omp for schedule(static) nowait
#endif // _OPENMP
        for (size_t c = 0; c < kInCellCount; ++c)
        {
            // Build a Terahedron containing the coordinates of the vertices.
            Tetrahedron<T> inTet = Tetrahedron<T>(
                    kVertices[oldToNewVertexMap.at(kElements[c]->vertices[0])]->coords,
                    kVertices[oldToNewVertexMap.at(kElements[c]->vertices[1])]->coords,
                    kVertices[oldToNewVertexMap.at(kElements[c]->vertices[2])]->coords,
                    kVertices[oldToNewVertexMap.at(kElements[c]->vertices[3])]->coords,
    				oldToNewVertexMap.at(kElements[c]->vertices[0]),
    				oldToNewVertexMap.at(kElements[c]->vertices[1]),
    				oldToNewVertexMap.at(kElements[c]->vertices[2]),
    				oldToNewVertexMap.at(kElements[c]->vertices[3]));

            // Generate the tets
            tetRefiner.refine(inTet,
            		kInVertexCount + c*additionalVertices,
					newTetsTmp, newVerticesTmp);

            // Copy new vertices
            for (unsigned int i = 0; i < additionalVertices; i++) {
            	memcpy(&newVertices[(c*additionalVertices + i) * 3],
            			glm::value_ptr(newVerticesTmp[i]), sizeof(T)*3);
            }

            // Copy tets
            for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
            	m_cells[(c*kSubCellsPerCell + i) * 4] = newTetsTmp[i].i;
            	m_cells[(c*kSubCellsPerCell + i) * 4 + 1] = newTetsTmp[i].j;
            	m_cells[(c*kSubCellsPerCell + i) * 4 + 2] = newTetsTmp[i].k;
            	m_cells[(c*kSubCellsPerCell + i) * 4 + 3] = newTetsTmp[i].l;
            }
        }

        delete [] newVerticesTmp;
        delete [] newTetsTmp;
#ifdef _OPENMP
    }
#endif
};

template<typename T>
MeshRefiner<T>::~MeshRefiner()
{
	delete [] m_cells;
	delete [] m_vertices;
}

//------------------------------------------------------------------------------

template<typename T>
const unsigned int* MeshRefiner<T>::getCellData() const {
    return &m_cells[0];
}

//------------------------------------------------------------------------------

template<typename T>
const T* MeshRefiner<T>::getVertexData() const {
    return &m_vertices[0];
}

//------------------------------------------------------------------------------

template<typename T>
std::size_t MeshRefiner<T>::getNumCells() const {
    return m_numSubCells;
}

//------------------------------------------------------------------------------

template<typename T>
std::size_t MeshRefiner<T>::getNumVertices() const {
    return m_numVertices;
}

//------------------------------------------------------------------------------

} // namespace
}

#endif // MESH_REFINER_H_
