/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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

#ifndef RESULTWRITER_ASYNCCELLIDS_H
#define RESULTWRITER_ASYNCCELLIDS_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "SeisSol.h"

namespace seissol
{

/**
 * This class can fix cells (vertex ids) in asynchronous mode.
 *
 * Sicne cells assume local vertex ids, we have to add an additional
 * when using the asynchronous MPI mode.
 *
 * @tparam CellVertices Number of vertices per cell
 */
template<int CellVertices>
class AsyncCellIDs
{
private:
	/** Null, if MPI is not enabled */
	unsigned int* m_cells;

	const unsigned int* m_constCells;

public:
	AsyncCellIDs(unsigned int nCells, unsigned int nVertices, const unsigned int* cells)
		: m_cells(0L)
	{
#ifdef USE_MPI
		// Add the offset to the cells
		MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().groupComm();
		unsigned int offset = nVertices;
		MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
		offset -= nVertices;

		// Add the offset to all cells
		m_cells = new unsigned int[nCells * CellVertices];
		for (unsigned int i = 0; i < nCells * CellVertices; i++)
			m_cells[i] = cells[i] + offset;
		m_constCells = m_cells;
#else // USE_MPI
		m_constCells = cells;
#endif // USE_MPI
	}

	~AsyncCellIDs()
	{
		delete [] m_cells;
	}

	const unsigned int* cells() const
	{
		return m_constCells;
	}
};

}

#endif // RESULTWRITER_ASYNCCELLIDS_H