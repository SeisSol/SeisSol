/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
#include <cmath>

#include <hdf5.h>

#include <Eigen/Dense>

#include "utils/args.h"

#include "Geometry/Refinement/RefinerUtils.h"

template<typename T>
T pow2(T v)
{
	return v*v;
}

int main(int argc, char* argv[])
{
	utils::Args args;
	args.addOption("diff", 'd', "Maximum allowed difference (Default 0.001)", utils::Args::Required, false);
	args.addOption("verbose", 'V', "Enable verbose output", utils::Args::No, false);
	args.addAdditionalOption("output.h5", "The SeisSol wave field output file");

	switch (args.parse(argc, argv)) {
	case utils::Args::Error:
		return 1;
	case utils::Args::Help:
		return 127;
	default:
		break;
	}

	bool verbose = args.getArgument("verbose", false);

	hid_t hdfFile = H5Fopen(args.getAdditionalArgument<const char*>("output.h5"),
			H5F_ACC_RDONLY, H5P_DEFAULT);

	// Read vertices
	hid_t vertData = H5Dopen(hdfFile, "/geometry", H5P_DEFAULT);
	hid_t space = H5Dget_space(vertData);
	hsize_t dims[2]; // Assuming to dimensions
	H5Sget_simple_extent_dims(space, dims, 0L);

	double* vertices = new double[dims[0]*dims[1]];
	H5Dread(vertData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertices);

	H5Sclose(space);
	H5Dclose(vertData);

	// Read cells
	hid_t cellData = H5Dopen(hdfFile, "/connect", H5P_DEFAULT);
	space = H5Dget_space(cellData);
	H5Sget_simple_extent_dims(space, dims, 0L);

	unsigned long* cells = new unsigned long[dims[0]*dims[1]];
	H5Dread(cellData, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, cells);

	hsize_t numCells = dims[0];

	H5Sclose(space);
	H5Dclose(cellData);

	// Read the 9 variables
	double* var[9];
	const char* varNames[9] = {"/sigma_xx", "/sigma_yy", "/sigma_zz", "/sigma_xy", "/sigma_yz", "/sigma_xz", "/u", "/v", "/w"};

	hid_t memSpace = H5Screate_simple(1, &numCells, &numCells);
	H5Sselect_all(memSpace);

	for (unsigned int i = 0; i < 9; i++) {
		hid_t varData = H5Dopen(hdfFile, varNames[i], H5P_DEFAULT);
		space = H5Dget_space(varData);
		H5Sget_simple_extent_dims(space, dims, 0L);

		assert(dims[1] == numCells);
		var[i] = new double[dims[1]];

		hsize_t start[2] = {0, 0};
		hsize_t count[2] = {1, numCells};
		H5Sselect_hyperslab(space, H5S_SELECT_SET, start, 0L, count, 0L);

		H5Dread(varData, H5T_NATIVE_DOUBLE, memSpace, space, H5P_DEFAULT, var[i]);

		H5Sclose(space);
		H5Dclose(varData);
	}

	H5Sclose(memSpace);

	H5Fclose(hdfFile);

	// Compute parameters
	const double maxDiff = args.getArgument("diff", 0.001);

      Eigen::Vector3d vec_n(1, 1, 1);
      const double n = 1.0/(vec_n.norm());
      const double mu = 1.0;
      const double lambda = 2.0;
      const double rho = 1.0;

      const double r1[9] = {
                      rho * (-2.0*pow2(n)*mu - 2.0*pow2(n)*mu + lambda + 2.0*mu),
			rho * (2.0*pow2(n)*mu + lambda),
			rho * (2.0*pow2(n)*mu + lambda),
			2.0 * n * mu * n * rho,
			2.0 * mu * n * rho * n,
			2.0 * mu * n * rho * n,
			n * sqrt(rho * (lambda + 2.0*mu)),
			n * sqrt(rho * (lambda + 2.0*mu)),
			sqrt(rho * (lambda + 2.0*mu)) * n
	};

	const double r8[9] = {
			2.0 * mu * n * rho * pow2(n) * n,
			-2.0 * mu * n * rho * pow2(n) * n,
			0,
			mu * rho * n * (2.0*pow2(n) + pow2(n) - 1) * n,
			-pow2(n) * mu * rho * pow2(n),
			mu * n * pow2(n) * rho * n,
			-n * n * n * sqrt(rho*mu),
			pow2(n) * n * sqrt(rho*mu),
			0
	};

	const Eigen::Vector3d k = 2.0 * M_PI / 100 * Eigen::Vector3d(1, 1, 1);

	// Check cells
	bool failed = false;

	for (unsigned int i = 0; i < numCells; i++) {
		const Eigen::Vector3d a = Eigen::Vector3d(&vertices[cells[i*4] * 3]);
		const Eigen::Vector3d b = Eigen::Vector3d(&vertices[cells[i*4 + 1] * 3]);
		const Eigen::Vector3d c = Eigen::Vector3d(&vertices[cells[i*4 + 2] * 3]);
		const Eigen::Vector3d d = Eigen::Vector3d(&vertices[cells[i*4 + 3] * 3]);

		seissol::refinement::Tetrahedron<double> tet(a, b, c, d);
		const Eigen::Vector3d center = tet.center();

		const double kx = k.dot(center);

		for (unsigned int j = 0; j < 9; j++) {
			const double value = r1[j]*sin(kx) + r8[j]*sin(kx);

			if (fabs(var[j][i] - value) > maxDiff) {
				failed = true;
				if (verbose)
					std::cout << var[j][i] << ' ' << value << ' ' << fabs(var[j][i] - value) << ' ';
			} else {
				if (verbose)
					std::cout << "ok" << ' ';
			}
		}

		if (verbose)
			std::cout << std::endl;

		if (failed && !verbose)
			// We can stop in this case
			break;
	}

	// Cleanup
	for (unsigned int i = 0; i < 9; i++) {
		delete [] var[i];
	}

	delete [] vertices;
	delete [] cells;

	return failed ? 1 : 0;
}
