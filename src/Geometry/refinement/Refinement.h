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
 * Base class for triangle/tetrahedron refinement
 **/

#ifndef REFINEMENT_REFINEMENT_H
#define REFINEMENT_REFINEMENT_H

#include <cassert>
#include <cstring>
#include <vector>
#include "RefinerUtils.h"

namespace refinement
{

//------------------------------------------------------------------------------

/**
 * Base class for triangle/tetrahedron refinement
 *
 * @todo free mesh resources after initialization
 */
template<class T>
class Refinement
{
protected:
    const unsigned int m_numVariables;

	/** Vertices of the mesh */
    std::vector<T> m_vertices;

	/** Cells of the original mesh */
    std::vector<unsigned int> m_cells;

	/**
	 * Set a vertex
	 *
	 * @param n The id of the vertex
	 * @param coords The coordinates of the vertex
	 */
	void setVertex(unsigned int n, const T coords[3])
	{
        setVertex(n, coords[0], coords[1], coords[2]);
	}

	/**
	 * Set a vertex
	 *
	 * @param n The id of the vertex
	 * @param x The x coordinate of the vertex
	 * @param y The y coordinate of the vertex
	 * @param z The z coordinate of the vertex
	 */
    void setVertex(unsigned int n, T x, T y, T z)
    {
		assert(n < m_vertices.size());
        m_vertices[n*3]   = x;
        m_vertices[n*3+1] = y;
        m_vertices[n*3+2] = z;
    }

	/**
	 * Set a cell
	 *
	 * @param n The id of the cell
	 * @param vertices The vertex ids
	 */
	void setCell(unsigned int n, const int vertices[4])
	{
        setCell(n, vertices[0], vertices[1], vertices[2], vertices[3]);
	}

	/**
	 * Set a cell
	 *
	 * @param n The id of the cell
	 * @param a The id of the vertex a
	 * @param b The id of the vertex b
	 * @param c The id of the vertex c
	 * @param d The id of the vertex d
	 */
	void setCell(unsigned int n, int a, int b, int c, int d)
	{
		assert(n*4 < m_cells.size());
		m_cells[n*4]   = a;
		m_cells[n*4+1] = b;
		m_cells[n*4+2] = c;
		m_cells[n*4+3] = d;
	}

public:
    Refinement(unsigned int numVariables) : m_numVariables(numVariables)
    {}

    virtual ~Refinement() {};

	/**
	 * @return Number of vertices in the refined mesh
	 */
	unsigned int nVertices() const
	{
		return m_vertices.size() / 3;
	}

	/**
	 * @return Number of cells in the refined mesh
	 */
	unsigned int nCells() const
	{
		return m_cells.size() / 4;
	}

	/**
	 * @return Vertices of the refined mesh
	 */
	const T* vertices() const
	{
		return m_vertices.data();
	}

	/**
	 * @return Cells of the refined mesh
	 */
	const unsigned int* cells() const
	{
		return m_cells.data();
	}

	/**
	 * Extract (and interpolate) one variable from the higher order data set
	 *
	 * @param idata Pointer to the (higher order) data from the original mesh
	 * @param cellMap Mapping from cell index to dof index
	 * @param variable The variable that should be extracted
	 * @param odata Pointer to a buffer where the flat data the refined mesh should be stored
	 */
	virtual void get(const double* idata, const unsigned int* cellMap,
			int variable, double* odata) const = 0;

	// TODO add getAll to extract all variables if required
    
};

//------------------------------------------------------------------------------

} // namespace

#endif // REFINEMENT_REFINEMENT_H
