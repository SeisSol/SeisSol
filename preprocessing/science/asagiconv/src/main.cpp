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
 **/

#include <cassert>
#include <fstream>
#include <limits>
#include <string>

#include <glm/common.hpp>
#include <glm/vec3.hpp>

#include <netcdf.h>

#include "utils/args.h"
#include "utils/logger.h"
#include "utils/stringutils.h"

/**
 * Check netCDF return value for errors
 */
static void checkNcError(int error)
{
	if (error != NC_NOERR)
		logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

/**
 * Computes the index for a coordinate
 */
template<typename T>
static glm::uvec3 computeIndex(const glm::tvec3<T> &coords,
		const glm::tvec3<T> &start, const glm::tvec3<T> &dist)
{
	return static_cast<glm::uvec3>(glm::round((coords - start) / dist));
}

static unsigned int compute1Dindex(const glm::uvec3 &index,
		const glm::uvec3 &size)
{
	return ((index.z * size.y) + index.y) * size.x + index.x;
}

static unsigned int computeOutIndex(unsigned int index, unsigned int variable,
		const glm::uvec3 &size, bool paraview)
{
	if (paraview)
		return index + variable*size.x*size.y*size.z;

	return index*3 + variable;
}

int main(int argc, char* argv[])
{
	// Parse command line arguments
	utils::Args args;
	args.addOption("max-z", 0, "The maximum z value (everything above will be set to this value)",
			utils::Args::Required, false);
	args.addOption("extend-z", 0, "Extend the grid in z direction to this value by duplication the top plane\n"
			"\tUse \"auto\" to automatically find the value",
			utils::Args::Required, false);
	args.addOption("paraview", 'p', "Viewable by Paraview", utils::Args::No, false);
	args.addAdditionalOption("input", "Input file");
	args.addAdditionalOption("output", "netCDF output file", false);

	const char* inputName;
	std::string outputName;

	// Parse/check command line arguments
	if (args.parse(argc, argv) == utils::Args::Success) {
		inputName = args.getAdditionalArgument<const char*>("input");

		std::string defaultOutputName(inputName);
		if (utils::StringUtils::endsWith(defaultOutputName, ".dat"))
			utils::StringUtils::replaceLast(defaultOutputName, ".dat", ".nc");
		else
			defaultOutputName += ".nc";

		outputName = args.getAdditionalArgument("output", defaultOutputName);
	} else {
		return 1;
	}

	// Initialize the input file
	std::ifstream input(inputName);

	// Header
	unsigned int lines;
	glm::uvec3 size;
	input >> lines;
	input >> size.x >> size.y >> size.z;

	// Read the rest of the line
	std::string header;
	getline(input, header);

	if (lines != size.x*size.y*size.z)
		logError() << "Number of lines in the file does not match with dimensions";

	logInfo() << "Size in input file:" << utils::nospace
			<< "[" << size.x << ", " << size.y << ", " << size.z << "]";

	// Allocate arrays
	glm::fvec3* vertices = new glm::fvec3[lines];
	float* rho = new float[lines];
	float* mu = new float[lines];
	float* lambda = new float[lines];

	// Read data
	logInfo() << "Reading data";
	for (unsigned int i = 0; i < lines; i++) {
		input >> vertices[i].x >> vertices[i].y >> vertices[i].z;
		input >> rho[i] >> mu[i] >> lambda[i];
	}

	float maxZ = args.getArgument("max-z", std::numeric_limits<float>::infinity());

	std::string extendZStr = args.getArgument<std::string>("extend-z", "no");
	float extendZ;
	bool autoExtendZ = false;
	if (extendZStr == "no")
		extendZ = -std::numeric_limits<float>::infinity();
	else if (extendZStr == "auto") {
		extendZ = -std::numeric_limits<float>::infinity();
		autoExtendZ = true;
	} else {
		extendZ = utils::StringUtils::parse<float>(extendZStr);
	}

	// Find max/min
	glm::fvec3 minVert(
			std::numeric_limits<float>::infinity(),
			std::numeric_limits<float>::infinity(),
			std::numeric_limits<float>::infinity());
	glm::fvec3 maxVert(
			-std::numeric_limits<float>::infinity(),
			-std::numeric_limits<float>::infinity(),
			-std::numeric_limits<float>::infinity());
	for (unsigned int i = 0; i < lines; i++) {
		// Auto extend?
		if (autoExtendZ)
			extendZ = std::max(extendZ, vertices[i].z);

		// Fix max
		if (vertices[i].z > maxZ)
			vertices[i].z = maxZ;

		minVert = glm::min(minVert, vertices[i]);
		maxVert = glm::max(maxVert, vertices[i]);
	}
	logInfo() << "Preliminary grid dimension in file:" << utils::nospace
			<< "[" << minVert.x << ", " << minVert.y << ", " << minVert.z << "] x ["
			<< maxVert.x << ", " << maxVert.y << ", " << maxVert.z << "]";

	// Compute distance between points
	glm::fvec3 avgDist = (maxVert - minVert)
			/ static_cast<glm::fvec3>(size - glm::uvec3(1, 1, 1));
	logInfo() << "Preliminary average distance:" << utils::nospace
			<< "[" << avgDist.x << ", " << avgDist.y << ", " << avgDist.z << "]";

#if 1
	// Fix real min/max
	// TODO maybe we need to do this more than once
	glm::fvec3 sumMin, sumMax;
	glm::uvec3 countMin, countMax;
	for (unsigned int i = 0; i < lines; i++) {
		for (unsigned int j = 0; j < 3; j++) {
			if (vertices[i][j] < minVert[j] + avgDist[j]) {
				sumMin[j] += vertices[i][j];
				countMin[j]++;
			}
			if (vertices[i][j] > maxVert[j] - avgDist[j]) {
				sumMax[j] += vertices[i][j];
				countMax[j]++;
			}
		}
	}

	minVert = sumMin / static_cast<glm::fvec3>(countMin);
	maxVert = sumMax / static_cast<glm::fvec3>(countMax);
	logInfo() << "Final grid dimension in file:" << utils::nospace
			<< "[" << minVert.x << ", " << minVert.y << ", " << minVert.z << "] x ["
			<< maxVert.x << ", " << maxVert.y << ", " << maxVert.z << "]";

	// Compute real distance between points
	avgDist = (maxVert - minVert)
			/ static_cast<glm::fvec3>(size - glm::uvec3(1, 1, 1));
	logInfo() << "Final average distance:" << utils::nospace
			<< "[" << avgDist.x << ", " << avgDist.y << ", " << avgDist.z << "]";
#endif

	// Extend grid in z dimension
	glm::uvec3 outputSize = size;
	if (maxVert.z < extendZ) {
		unsigned int planes = std::ceil((extendZ - maxVert.z) / avgDist.z);
		outputSize.z = outputSize.z + planes;
	}

	logInfo() << "Size in output file:" << utils::nospace
			<< "[" << outputSize.x << ", " << outputSize.y << ", " << outputSize.z << "]";

	// Open the netCDF file
	int ncFile;
	checkNcError(nc_create(outputName.c_str(), NC_NETCDF4, &ncFile));

	// Create dimensions
	int ncDims[3];
	checkNcError(nc_def_dim(ncFile, "x", outputSize.x, &ncDims[2]));
	checkNcError(nc_def_dim(ncFile, "y", outputSize.y, &ncDims[1]));
	checkNcError(nc_def_dim(ncFile, "z", outputSize.z, &ncDims[0]));

	// Create dimension variables
	int ncX, ncY, ncZ;
	checkNcError(nc_def_var(ncFile, "x", NC_FLOAT, 1, &ncDims[2], &ncX));
	checkNcError(nc_def_var(ncFile, "y", NC_FLOAT, 1, &ncDims[1], &ncY));
	checkNcError(nc_def_var(ncFile, "z", NC_FLOAT, 1, &ncDims[0], &ncZ));

	// Fill dimension variables
	float* x = new float[outputSize.x];
	for (unsigned int i = 0; i < outputSize.x; i++) {
		x[i] = minVert.x + avgDist.x * i;
	}
	checkNcError(nc_put_var_float(ncFile, ncX, x));
	delete [] x;

	float* y = new float[outputSize.y];
	for (unsigned int i = 0; i < outputSize.y; i++) {
		y[i] = minVert.y + avgDist.y * i;
	}
	checkNcError(nc_put_var_float(ncFile, ncY, y));
	delete [] y;

	float* z = new float[outputSize.z];
	for (unsigned int i = 0; i < outputSize.z; i++) {
		z[i] = minVert.z + avgDist.z * i;
	}
	checkNcError(nc_put_var_float(ncFile, ncZ, z));
	delete [] z;

	// Viewable by Paraview?
	bool paraview = args.getArgument("paraview", false);

	// Create variables
	int ncData, ncRho, ncMu, ncLambda;
	if (paraview) {
		checkNcError(nc_def_var(ncFile, "rho", NC_FLOAT, 3, ncDims, &ncRho));
		checkNcError(nc_def_var(ncFile, "mu", NC_FLOAT, 3, ncDims, &ncMu));
		checkNcError(nc_def_var(ncFile, "lambda", NC_FLOAT, 3, ncDims, &ncLambda));
	} else {
		// Create compound type
		int ncType;
		checkNcError(nc_def_compound(ncFile, 3*sizeof(float), "material", &ncType));
		checkNcError(nc_insert_compound(ncFile, ncType, "rho", 0, NC_FLOAT));
		checkNcError(nc_insert_compound(ncFile, ncType, "mu", sizeof(float), NC_FLOAT));
		checkNcError(nc_insert_compound(ncFile, ncType, "lambda", 2*sizeof(float), NC_FLOAT));

		checkNcError(nc_def_var(ncFile, "data", ncType, 3, ncDims, &ncData));
	}

	// Save variables
	const float * const data[3] = {rho, mu, lambda};

	float* outputData = new float[3*outputSize.x*outputSize.y*outputSize.z];
	memset(outputData, 0, 3*outputSize.x*outputSize.y*outputSize.z*sizeof(float));
	for (int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < lines; j++) {
			glm::uvec3 index = computeIndex(vertices[j], minVert, avgDist);

			if (glm::any(glm::greaterThanEqual(index, size))) {
				logError() << "Index" << utils::nospace
						<< "[" << index.x << ", " << index.y << ", " << index.z << "] out of range";
			}

			unsigned int pos = compute1Dindex(index, outputSize);
			assert(pos < outputSize.x*outputSize.y*outputSize.z);

			pos = computeOutIndex(pos, i, outputSize, paraview);

			if (outputData[pos] != 0) {
				logInfo() << "Value at" << utils::nospace
						<< "[" << index.x << ", " << index.y << ", " << index.z << "] already set";
			}

			outputData[pos] = data[i][j];
		}

		// Fill extend
		for (unsigned int j = size.z; j < outputSize.z; j++) {
			for (unsigned int k = 0; k < outputSize.x*outputSize.y; k++) {
				unsigned int posRead = computeOutIndex((size.z-1)*outputSize.x*outputSize.y + k,
						i, outputSize, paraview);
				unsigned int posWrite = computeOutIndex(j*outputSize.x*outputSize.y + k,
						i, outputSize, paraview);

				outputData[posWrite] = outputData[posRead];
			}
		}

	}

	if (paraview) {
		checkNcError(nc_put_var_float(ncFile, ncRho, &outputData[0]));
		checkNcError(nc_put_var_float(ncFile, ncMu, &outputData[outputSize.x*outputSize.y*outputSize.z]));
		checkNcError(nc_put_var_float(ncFile, ncLambda, &outputData[2*outputSize.x*outputSize.y*outputSize.z]));
	} else {
		checkNcError(nc_put_var(ncFile, ncData, outputData));
	}

	delete [] outputData;

	// Cleanup
	checkNcError(nc_close(ncFile));

	delete [] rho;
	delete [] mu;
	delete [] lambda;

	logInfo() << "Finished";
}
