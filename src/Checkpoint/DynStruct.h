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
 * Dynamic structure
 */

#ifndef DYN_STRUCT_H
#define DYN_STRUCT_H

#include "Parallel/MPI.h"

#include <cstdlib>
#include <cstring>
#include <vector>

#include "utils/logger.h"

namespace seissol
{

/**
 * @todo Set MPI and HDF5 type
 */
class DynStruct
{
private:
	/**
	 * Component description
	 */
	struct Component
	{
		size_t offset;
	};

	std::vector<Component> m_components;

	/** Total size of the struct */
	size_t m_size;

	void* m_buffer;

public:
	DynStruct()
		: m_size(0), m_buffer(0L)
	{ }

	~DynStruct()
	{
		free(m_buffer);
	}

	/**
	 * Add a component to the struct
	 *
	 * @return The id of the component
	 */
	template<typename T>
	size_t add()
	{
		if (m_buffer)
			logError() << "The dyn struct was already initialized. Could not add an additional component.";

		Component comp;
		comp.offset = m_size;
		m_components.push_back(comp);

		m_size += sizeof(T);

		return m_components.size()-1;
	}

	void alloc(size_t alignment = 0)
	{
		if (m_buffer)
			logError() << "The dyn struct was already initialized.";

		if (alignment) {
			m_size = (m_size + alignment - 1) / alignment;
			m_size *= alignment;

			if (posix_memalign(&m_buffer, alignment, m_size) != 0)
				logError() << "Could not allocate buffer";
		} else {
			m_buffer = malloc(m_size);
		}
	}

	template<typename T>
	T& value(size_t i)
	{
		return *static_cast<T*>(m_buffer + m_components[i].offset);
	}

	template<typename T>
	const T& value(size_t i) const
	{
		return *static_cast<T*>(m_buffer + m_components[i].offset);
	}

	size_t size() const
	{
		return m_size;
	}

	void* data()
	{
		return m_buffer;
	}

	const void* data() const
	{
		return m_buffer;
	}

	/**
	 * Set all contents to 0
	 */
	void clear()
	{
		if (m_buffer)
			memset(m_buffer, 0, m_size);
		else
			logWarning(seissol::MPI::mpi.rank()) << "Trying to clear an unallocated struct.";
	}
};

}

#endif // DYN_STRUCT_H