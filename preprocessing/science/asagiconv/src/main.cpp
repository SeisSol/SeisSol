/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
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

#include <cstring>
#include <limits>
#include <string>

#include <pthread.h>

#include <pstream.h>

#include <Eigen/Dense>

#include <netcdf.h>

#include <proj_api.h>

#include "utils/args.h"
#include "utils/logger.h"
#include "utils/stringutils.h"

const static float MIN_VALID = -10000;

typedef Eigen::Array<unsigned long, 3, 1> ularr3;
typedef Eigen::Array<unsigned long, 2, 1> ularr2;

struct inputData_t
{
        Eigen::Array3d minVert;
        ularr3 gridSize;
	Eigen::Array3d gridDist;

	const char* meshProj;
	bool meshIsGeo;
	const char* queryProj;
	bool queryIsGeo;

	redi::rpstream& query;
};

/**
 * Check netCDF return value for errors
 */
static void checkNcIError(int error)
{
	if (error != NC_NOERR)
		logError() << "Error while reading netCDF file:" << nc_strerror(error);
}

/**
 * Check netCDF return value for errors
 */
static void checkNcOError(int error)
{
	if (error != NC_NOERR)
		logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

/**
 * The input thread for the query program
 */
static void* runInput(void* p)
{
	inputData_t* i = reinterpret_cast<inputData_t*>(p);

	// Project from mesh to sample code
	unsigned long totalArea = i->gridSize(0) * i->gridSize(1);
	double* projX = new double[totalArea];
	double* projY = new double[totalArea];

	const bool meshIsGeo = i->meshIsGeo;

	for (unsigned long y = 0; y < i->gridSize(1); y++) {
		for (unsigned long x = 0; x < i->gridSize(0); x++) {
			ularr2 pos(x, y);
                        Eigen::Array2d coord = i->minVert.block<2,1>(0,0) + pos.cast<double>() * i->gridDist.block<2,1>(0,0);
			if (meshIsGeo)
				coord *= Eigen::Array2d(DEG_TO_RAD, DEG_TO_RAD);

			projX[y*i->gridSize(0) + x] = coord(0);
			projY[y*i->gridSize(0) + x] = coord(1);
		}
	}

	// Project coordinates
	projPJ pjMesh = pj_init_plus(i->meshProj);
	projPJ pjQuery = pj_init_plus(i->queryProj);
	if (pj_transform(pjMesh, pjQuery, totalArea, 1, projX, projY, 0L) != 0)
		logError() << "Coordinate transformation failed";

	// Convert to degree
	if (i->queryIsGeo) {
		for (unsigned long x = 0; x < totalArea; x++) {
			projX[x] *= RAD_TO_DEG;
			projY[x] *= RAD_TO_DEG;
		}
	}

	// Send the input data
	i->query << (totalArea * i->gridSize(2)) << std::endl;
	for (unsigned long z = 0; z < i->gridSize(2); z++) {
		double depth = i->minVert(2) + z * i->gridDist(2);
		for (unsigned long x = 0; x < totalArea; x++)
			i->query << projX[x] << ' ' << projY[x] << ' ' << depth << std::endl;
	}

	// Send eof
	i->query << redi::peof;

	delete [] projX;
	delete [] projY;

	return 0L;
}

static float vs2mu(float vs, float rho)
{
	return vs*vs * rho;
}

static float vp2lambda(float vp, float rho, float mu)
{
	return vp*vp * rho - 2*mu;
}

int main(int argc, char* argv[])
{
	// Parse command line arguments
	utils::Args args;
	args.addOption("mesh", 'm', "The mesh for which the data should be collected");
	args.addOption("border", 'b', "Border around the mesh for which data should be included",
			utils::Args::Required, false);
	args.addOption("size", 's', "Size of the out grid (Format x:y:z)");
	args.addOption("query", 'q', "Path to the query program (Default: ./scripts/vx_lite_wrapper)",
			utils::Args::Required, false);
//	args.addOption("elevation", 'e', "Path to x_SRTM3, y_SRTM3, z_SRTM3 containing an elevation map\n"
//			"                              If empty, elevation instead of elevation offset is used for querying.",
//			utils::Args::Required, false);
	// See http://spatialreference.org/
	args.addOption("mesh-proj", 0, "Mesh coordinate projection (Default: UTM11)",
			utils::Args::Required, false);
	args.addOption("mesh-proj-geo", 0, "Mesh coordinates are geographic locations (Default: No)",
			utils::Args::Required, false);
	args.addOption("query-proj", 0, "Query coordinate projection (Default: WGS84)",
			utils::Args::Required, false);
	args.addOption("query-proj-geo", 0, "Query coordinates are geographic locations (Default: Yes)",
			utils::Args::Required, false);
	args.addOption("paraview", 'p', "Viewable by Paraview", utils::Args::No, false);
	args.addOption("chunk-size", 'c', "Chunk size for the netCDF file",
			utils::Args::Required, false);
	args.addAdditionalOption("output", "netCDF output file");

	// Parse/check command line arguments
	if (!args.parse(argc, argv) == utils::Args::Success)
		return 1;

	// Get the grid size
	std::vector<std::string> sizeVec = utils::StringUtils::split(args.getArgument<std::string>("size"), ':');
	if (sizeVec.size() != 3)
		logError() << "Size format must be \"x:y:z\"";
	ularr3 gridSize(utils::StringUtils::parse<unsigned long>(sizeVec[0]),
			utils::StringUtils::parse<unsigned long>(sizeVec[1]),
			utils::StringUtils::parse<unsigned long>(sizeVec[2]));

	logInfo() << "Grid size:" << utils::nospace
			<< "[" << gridSize(0) << ", " << gridSize(1) << ", " << gridSize(2) << "]";

	// Read the coordinates from the mesh
	const char* mesh = args.getArgument<const char*>("mesh");

	int ncFile;
	checkNcIError(nc_open(mesh, NC_NOWRITE, &ncFile));

	// Get the number of partitions
	int ncDimPart;
	checkNcIError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
	size_t partitions;
	checkNcIError(nc_inq_dimlen(ncFile, ncDimPart, &partitions));

	// Get max number of vertices
	int ncDimVertices;
	checkNcIError(nc_inq_dimid(ncFile, "vertices", &ncDimVertices));
	size_t maxVertices;
	checkNcIError(nc_inq_dimlen(ncFile, ncDimVertices, &maxVertices));

	// Get min/max vertices
	Eigen::Array3d minVert(
			std::numeric_limits<double>::infinity(),
			std::numeric_limits<double>::infinity(),
			std::numeric_limits<double>::infinity());
	Eigen::Array3d maxVert(
			-std::numeric_limits<double>::infinity(),
			-std::numeric_limits<double>::infinity(),
			-std::numeric_limits<double>::infinity());

	double* vertices = new double[maxVertices*3];

	int ncVarVrtxSize;
	checkNcIError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
	int ncVarVrtxCoords;
	checkNcIError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
	for (size_t i = 0; i < partitions; i++) {
		size_t start[3] = {i, 0, 0};
		int size;
		checkNcIError(nc_get_var1_int(ncFile, ncVarVrtxSize, start, &size));

		size_t count[3] = {1, static_cast<size_t>(size), 3};
		checkNcIError(nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count, vertices));

		for (int j = 0; j < size; j++) {
			Eigen::Array3d v = Eigen::Map<Eigen::Array3d>(&vertices[j*3]);

			minVert = minVert.min(v);
			maxVert = maxVert.max(v);
		}
	}

	checkNcIError(nc_close(ncFile));

	delete [] vertices;

	// Add borders
	double border = args.getArgument("border", 0.);
	minVert -= Eigen::Array3d(border, border, border);
	maxVert += Eigen::Array3d(border, border, border);

	logInfo() << "Grid dimension (incl. border):" << utils::nospace
			<< "[" << minVert(0) << ", " << minVert(1) << ", " << minVert(2) << "] x ["
			<< maxVert(0) << ", " << maxVert(1) << ", " << maxVert(2) << "]";

	// Distance between grid points
	Eigen::Vector3d gridDist = (maxVert - minVert)
			* (gridSize - ularr3(1, 1, 1)).cast<double>().inverse();

	logInfo() << "Grid interval:" << utils::nospace
			<< "[" << gridDist(0) << ", " << gridDist(1) << ", " << gridDist(2) << "]";

	// Create query process
	redi::rpstream pstream(args.getArgument("query", "./scripts/vx_lite_wrapper"),
			redi::pstreambuf::pstdin | redi::pstreambuf::pstdout);

	// Start the input and error thread
	inputData_t inputData {
		minVert,
		gridSize,
		gridDist,
		args.getArgument("mesh-proj", "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
		args.getArgument("mesh-proj-geo", false),
		args.getArgument("query-proj", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
		args.getArgument("query-proj-geo", true),
		pstream
	};

	pthread_t inputThread;
	pthread_create(&inputThread, 0L, runInput, &inputData);

	// Viewable by Paraview?
	bool paraview = args.getArgument("paraview", false);

	// Open the netCDF file
	const char* outputName = args.getAdditionalArgument<const char*>("output");
	checkNcOError(nc_create(outputName, NC_NETCDF4, &ncFile));

	// Create dimensions
	int ncDims[3];
	checkNcOError(nc_def_dim(ncFile, "x", gridSize(0), &ncDims[2]));
	checkNcOError(nc_def_dim(ncFile, "y", gridSize(1), &ncDims[1]));
	checkNcOError(nc_def_dim(ncFile, "z", gridSize(2), &ncDims[0]));

	// Create dimension variables
	int ncX, ncY, ncZ;
	checkNcOError(nc_def_var(ncFile, "x", NC_FLOAT, 1, &ncDims[2], &ncX));
	checkNcOError(nc_def_var(ncFile, "y", NC_FLOAT, 1, &ncDims[1], &ncY));
	checkNcOError(nc_def_var(ncFile, "z", NC_FLOAT, 1, &ncDims[0], &ncZ));

	// Create variables
	int ncData, ncRho, ncMu, ncLambda;
	if (paraview) {
		checkNcOError(nc_def_var(ncFile, "rho", NC_FLOAT, 3, ncDims, &ncRho));
		checkNcOError(nc_def_var(ncFile, "mu", NC_FLOAT, 3, ncDims, &ncMu));
		checkNcOError(nc_def_var(ncFile, "lambda", NC_FLOAT, 3, ncDims, &ncLambda));
	} else {
		// Create compound type
		int ncType;
		checkNcOError(nc_def_compound(ncFile, 3*sizeof(float), "material", &ncType));
		checkNcOError(nc_insert_compound(ncFile, ncType, "rho", 0, NC_FLOAT));
		checkNcOError(nc_insert_compound(ncFile, ncType, "mu", sizeof(float), NC_FLOAT));
		checkNcOError(nc_insert_compound(ncFile, ncType, "lambda", 2*sizeof(float), NC_FLOAT));

		checkNcOError(nc_def_var(ncFile, "data", ncType, 3, ncDims, &ncData));

		if (args.isSet("chunk-size")) {
			unsigned int chunkSize = args.getArgument<unsigned int>("chunk-size");

			size_t chunks[3] = {chunkSize, chunkSize, chunkSize};
			checkNcOError(nc_def_var_chunking(ncFile, ncData, NC_CHUNKED, chunks));
		}
	}

	checkNcOError(nc_enddef(ncFile));

	// Fill dimension variables
	float* x = new float[gridSize(0)];
	for (unsigned int i = 0; i < gridSize(0); i++) {
		x[i] = minVert(0) + gridDist(0) * i;
	}
	checkNcOError(nc_put_var_float(ncFile, ncX, x));
	delete [] x;

	float* y = new float[gridSize(1)];
	for (unsigned int i = 0; i < gridSize(1); i++) {
		y[i] = minVert(1) + gridDist(1) * i;
	}
	checkNcOError(nc_put_var_float(ncFile, ncY, y));
	delete [] y;

	float* z = new float[gridSize(2)];
	for (unsigned int i = 0; i < gridSize(2); i++) {
		z[i] = minVert(2) + gridDist(2) * i;
	}
	checkNcOError(nc_put_var_float(ncFile, ncZ, z));
	delete [] z;

	// Read the data from the query tool
	unsigned long totalArea = gridSize(0) * gridSize(1);
	float* data = new float[totalArea*3];

	// Initialize data with invalid values
	for (unsigned long i = 0; i < totalArea*3; i++)
		data[i] = std::numeric_limits<float>::infinity();

	for (unsigned long z = 0; z < gridSize(2); z++) {
		for (unsigned long i = 0; i < totalArea; i++) {
			// Ignore the coordinates
			float x, y, z;
			pstream.out() >> x >> y >> z;

			float vp, vs, rho;
			pstream.out() >> vp >> vs >> rho;

			float mu, lambda;
			if (vp >= MIN_VALID && vs >= MIN_VALID && rho >= MIN_VALID) {
				mu = vs2mu(vs, rho);
				lambda = vp2lambda(vp, rho, mu);

				if (paraview) {
					data[i] = rho;
					data[i + totalArea] = mu;
					data[i + totalArea*2] = lambda;
				} else {
					data[i*3] = rho;
					data[i*3 + 1] = mu;
					data[i*3 + 2] = lambda;
				}
			} else {
				if (z == 0)
					logError() << "Invalid value in deepest level.";
			}
		}

		size_t start[3] = {z, 0, 0};
		size_t count[3] = {1, gridSize(1), gridSize(0)};

		if (paraview) {
			checkNcOError(nc_put_vara_float(ncFile, ncRho, start, count, data));
			checkNcOError(nc_put_vara_float(ncFile, ncMu, start, count, &data[totalArea]));
			checkNcOError(nc_put_vara_float(ncFile, ncLambda, start, count, &data[totalArea*2]));
		} else {
			checkNcOError(nc_put_vara(ncFile, ncData, start, count, data));
		}
	}

	pthread_join(inputThread, 0L);

	// Cleanup
	checkNcOError(nc_close(ncFile));

	delete [] data;

	logInfo() << "Finished";
}
