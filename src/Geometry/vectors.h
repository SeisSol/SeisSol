/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef UTILS_VECTORS_H_
#define UTILS_VECTORS_H_

#include <algorithm>
#include <vector>

/**
 * Find the first common element in all three vectors. If now element is found, defaultValue is returned.
 */
template<typename T>
T intersection(const std::vector<T> &v1, const std::vector<T> &v2, const std::vector<T> &v3, T defaultValue = T())
{
	for (typename std::vector<T>::const_iterator i = v1.begin(); i != v1.end(); i++) {
		typename std::vector<T>::const_iterator j = std::find(v2.begin(), v2.end(), *i);

		if (j != v2.end()) {
			typename std::vector<T>::const_iterator k = std::find(v3.begin(), v3.end(), *i);

			if (k != v3.end())
				return *i;
		}
	}

	return defaultValue;
}

/**
 * Find the first common element in all three vectors. The first vector is given as begin and end
 * iterator position.
 *
 * If a common element was found, v1begin will point to this element, otherwise v1begin == v1end.
 */
template<typename T>
void intersection(typename std::vector<T>::const_iterator &v1begin, const typename std::vector<T>::const_iterator &v1end,
		const std::vector<T> &v2,
		const std::vector<T> &v3)
{
	for (; v1begin != v1end; v1begin++) {
		typename std::vector<T>::const_iterator j = std::find(v2.begin(), v2.end(), *v1begin);

		if (j != v2.end()) {
			typename std::vector<T>::const_iterator k = std::find(v3.begin(), v3.end(), *v1begin);

			if (k != v3.end())
				return;
		}
	}
}

#endif /* UTILS_VECTOR_H_ */
