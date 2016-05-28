/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
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
 * Stopwatch originally developed by A. Heinecke
 */

#ifndef STOPWATCH_H
#define STOPWATCH_H

#ifdef _WIN32
#include <windows.h>
#else // _WIN32
#include <sys/time.h>
#endif // _WIN32

/**
 * OS-independent (per Preprocessor) version of a stopwatch
 *
 * Part of SeisSol, so you can easily calculate the needed time of SeisSol computations with a high precision
 */
class Stopwatch
{
private:
#ifdef _WIN32
	LARGE_INTEGER m_ticksPerSecond;
	LARGE_INTEGER m_begin;
#else // _WIN32
	timeval m_begin;
#endif // _WIN32

public:
	/**
	 * Constructor
	 *
	 * resets the Stopwatch
	 */
	Stopwatch()
	{
#ifdef _WIN32
		QueryPerformanceFrequency(&m_ticksPerSecond);
#endif // _WIN32
	}

	/**
	 * Destructor
	 */
	~Stopwatch()
	{}

	/**
	 * starts the time measuring
	 */
	void start()
	{
#ifdef _WIN32
		QueryPerformanceCounter(&m_begin);
#else // _WIN32
		gettimeofday(&m_begin,(struct timezone *)0);
#endif // _WIN32
	}

	/**
	 * stops time measuring
	 *
	 * @return measured time in seconds
	 */
	double stop()
	{
#ifdef _WIN32
		LARGE_INTEGER end;
		QueryPerformanceCounter(&end);

		double ret, ticksps;

		end.QuadPart -= m_begin.QuadPart;
		ret = (double)(end.QuadPart);
		ticksps = (double)(m_ticksPerSecond.QuadPart);
		ret /= ticksps;

		return ret;
#else // _WIN32
		timeval end;

		gettimeofday(&end, (struct timezone *) 0);
		double seconds, useconds;
		double ret, tmp;

		if (end.tv_usec >= m_begin.tv_usec) {
			seconds = (double) end.tv_sec - (double) m_begin.tv_sec;
			useconds = (double) end.tv_usec - (double) m_begin.tv_usec;
		} else {
			seconds = (double) end.tv_sec - (double) m_begin.tv_sec;
			seconds -= 1;                                   // Correction
			useconds = (double) end.tv_usec - (double) m_begin.tv_usec;
			useconds += 1000000;                    // Correction
		}

		// get time in seconds
		tmp = (double) useconds;
		ret = (double) seconds;
		tmp /= 1000000;
		ret += tmp;

		return ret;
#endif // _WIN32
	}
};

#endif // STOPWATCH_H
