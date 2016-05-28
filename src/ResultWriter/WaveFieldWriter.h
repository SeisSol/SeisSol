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

#ifndef WAVE_FIELD_WRITER_H
#define WAVE_FIELD_WRITER_H

#include "Parallel/MPI.h"

#include <string>
#include <vector>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "Geometry/MeshReader.h"
#include "Geometry/refinement/RefinerUtils.h"
#include "Geometry/refinement/MeshRefiner.h"
#include "Geometry/refinement/VariableSubSampler.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol
{

class WaveFieldWriter
{
private:
    /** True if wave field output is enabled */
    bool m_enabled;

    /** The output prefix for the filename */
    std::string m_outputPrefix;

    /** The XMDF Writer used for the wave field */
    xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_waveFieldWriter;

    /** The XDMF Writer for low order data */
    xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_lowWaveFieldWriter;

    /** The variable subsampler for the refined mesh */
    refinement::VariableSubsampler<double>* m_variableSubsampler;

    /** Original number of cells */
    unsigned int m_numCells;

    /** Number of variables */
    unsigned int m_numVariables;

    /** Pointer to the degrees of freedom */
    const double* m_dofs;

	/** Pointer to the plastic strain */
	const double* m_pstrain;

    /** Mapping from the cell order to dofs order */
    const unsigned int* m_map;

    /** Time of the last output (makes sure output is not written twice at the end */
    double m_lastTimeStep;

    /** The tolerance in the time for ignoring duplicate time steps */
    double m_timeTolerance;

    /** Buffer required to extract the output data from the unknowns */
    double* m_outputBuffer;

    /** Flag indicated which variables should be written */
    std::vector<bool> m_output_flags;

public:
    WaveFieldWriter()
     	 : m_enabled(false),
		   m_waveFieldWriter(0L), m_lowWaveFieldWriter(0L),
		   m_variableSubsampler(0L),
		   m_numCells(0),
		   m_numVariables(0),
		   m_dofs(0L), m_pstrain(0L),
		   m_map(0L),
		   m_lastTimeStep(-1),
		   m_timeTolerance(0),
		   m_outputBuffer(0L)
    {
    }

    /**
     * Activate the wave field output
     */
    void enable()
    {
        m_enabled = true;
    }

    /**
     * @return True if wave field output is enabled, false otherwise
     */
    bool isEnabled() const
    {
        return m_enabled;
    }

    /**
     * Set the output prefix for the filename
     */
    void setFilename(const char* outputPrefix)
    {
        m_outputPrefix = outputPrefix;
    }

    /**
     * Initialize the wave field ouput
     *
     * @param map The mapping from the cell order to dofs order
     * @param timeTolerance The tolerance in the time for ignoring duplicate time steps
     */
    void init(int numVars, int order, int numAlignedDOF,
            const MeshReader &meshReader,
            const double* dofs,  const double* pstrain,
			const unsigned int* map,
            int refinement, int timestep,
			double timeTolerance)
    {
        if (!m_enabled)
            return;

        const int rank = MPI::mpi.rank();

		logInfo(rank) << "Initializing HDF5 wave field output.";

		if (m_waveFieldWriter != 0L)
			logError() << "Wave field writer already initialized";

		// Get the original number of cells (currently required for pstrain output)
		m_numCells = meshReader.getElements().size();

		//
		// High order I/O
		//

		m_numVariables = numVars;

		if (numVars != 9)
			logError()
					<< "XDMF output supports currently only 9 variables. Number of variables specified:"
					<< m_numVariables;
		std::vector<const char*> variables(m_numVariables);
		variables[0] = "sigma_xx";
		variables[1] = "sigma_yy";
		variables[2] = "sigma_zz";
		variables[3] = "sigma_xy";
		variables[4] = "sigma_yz";
		variables[5] = "sigma_xz";
		variables[6] = "u";
		variables[7] = "v";
		variables[8] = "w";

		// Currently all variables have to be chosen.
		m_output_flags.resize(numVars);
		std::fill(m_output_flags.begin(), m_output_flags.end(), true);

		// Setup the tetrahedron refinement strategy
		refinement::TetrahedronRefiner<double>* tetRefiner = 0L;
		switch (refinement) {
		case 0:
			logInfo(rank) << "Refinement is turned off.";
			tetRefiner = new refinement::IdentityRefiner<double>();
			break;
		case 1:
			logInfo(rank) << "Refinement Startegy is \"Divide by 4\"";
			tetRefiner = new refinement::DivideTetrahedronBy4<double>();
			break;
		case 2:
			logInfo(rank) << "Refinement Startegy is \"Divide by 8\"";
			tetRefiner = new refinement::DivideTetrahedronBy8<double>();
			break;
		case 3:
			logInfo(rank) << "Refinement Startegy is \"Divide by 32\"";
			tetRefiner = new refinement::DivideTetrahedronBy32<double>();
			break;
		default:
			logError() << "Refinement Strategy is invalid!" << std::endl
					<< "Current value : " << refinement << std::endl
					<< "Valid options are :" << std::endl << "0 - No Refinement"
					<< std::endl << "1 - Refinement by 4" << std::endl
					<< "2 - Refinement by 8" << std::endl
					<< "3 - Refinement by 32";
		}

		// Refine the mesh
        refinement::MeshRefiner<double> meshRefiner(meshReader, *tetRefiner);

		logInfo(rank) << "Refinement class initialized";
		logDebug() << "Cells : "
				<< meshReader.getElements().size() << "refined-to ->"
				<< meshRefiner.getNumCells();
		logDebug() << "Vertices : "
				<< meshReader.getVertices().size()
				<< "refined-to ->"
				<< meshRefiner.getNumVertices();

        // Initialize the I/O handler and write the mesh
		m_waveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
				rank, m_outputPrefix.c_str(), variables, timestep);
        m_waveFieldWriter->init(
                meshRefiner.getNumCells(), meshRefiner.getCellData(),
                meshRefiner.getNumVertices(), meshRefiner.getVertexData(),
                true);

        logInfo(rank) << "WaveFieldWriter initialized";

        // Initialize the variable subsampler
        m_variableSubsampler = new refinement::VariableSubsampler<double>(
                    meshReader.getElements().size(),
					*tetRefiner, order, numVars, numAlignedDOF);

        // Delete the tetRefiner since it is no longer required
        delete tetRefiner;

        logInfo(rank) << "VariableSubsampler initialized";

        // Create output buffer
        m_outputBuffer = new double[meshRefiner.getNumCells()];
        assert(meshRefiner.getNumCells() >= meshReader.getElements().size());

        //
        //  Low order I/O
        //

        if (pstrain) {
        	logInfo(rank) << "Initialize low order output";

        	// Refinement strategy (no refinement)
        	refinement::IdentityRefiner<double> lowTetRefiner;

        	// Mesh refiner
        	refinement::MeshRefiner<double> lowMeshRefiner(meshReader, lowTetRefiner);

        	// Variables
    		std::vector<const char*> lowVariables(7);
    		lowVariables[0] = "ep_xx";
    		lowVariables[1] = "ep_yy";
    		lowVariables[2] = "ep_zz";
    		lowVariables[3] = "ep_xy";
    		lowVariables[4] = "ep_yz";
    		lowVariables[5] = "ep_xz";
    		lowVariables[6] = "eta";

			// Initialize the low order I/O handler
			m_lowWaveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
					rank, (m_outputPrefix+"-low").c_str(), lowVariables, timestep);
        	m_lowWaveFieldWriter->init(
        			lowMeshRefiner.getNumCells(), lowMeshRefiner.getCellData(),
        			lowMeshRefiner.getNumVertices(), lowMeshRefiner.getVertexData(),
        			true);

        	logInfo(rank) << "Low order output initialized";
        }

        	// Save dof/map pointer
        m_dofs = dofs;
        m_pstrain = pstrain;
        m_map = map;

        logInfo(rank) << "Initializing HDF5 wave field output. Done.";
    }

    /**
     * @return The current time step of the wave field output
     */
    unsigned int timestep() const
    {
        if (!m_enabled)
            return 0;

        return m_waveFieldWriter->timestep();
    }

    /**
     * Write a time step
     */
    void write(double time)
    {
        EPIK_TRACER("WaveFieldWriter_write");
        SCOREP_USER_REGION("WaveFieldWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION);

        if (!m_enabled)
            return;

        const int rank = MPI::mpi.rank();

        if (time <= m_lastTimeStep + m_timeTolerance) {
        	// Ignore duplicate time steps. Might happen at the end of a simulation
        	logInfo(rank) << "Ignoring duplicate time step at time " << time;
        	return;
        }

        logInfo(rank) << "Writing wave field at time" << utils::nospace << time << '.';

        // High order output
        m_waveFieldWriter->addTimeStep(time);

        for (unsigned int i = 0; i < m_numVariables; i++) {
            if (!m_output_flags[i]) continue;

            m_variableSubsampler->get(m_dofs, m_map, i, m_outputBuffer);

            m_waveFieldWriter->writeData(i, m_outputBuffer);
        }

        m_waveFieldWriter->flush();

        // Low order output
        if (m_pstrain) {
        	m_lowWaveFieldWriter->addTimeStep(time);

        	for (unsigned int i = 0; i < 7; i++) {
#ifdef _OPENMP
				#pragma omp parallel for schedule(static)
#endif // _OPENMP
        		for (unsigned int j = 0; j < m_numCells; j++)
        			m_outputBuffer[j] = m_pstrain[m_map[j] * 7 + i];

        		m_lowWaveFieldWriter->writeData(i, m_outputBuffer);
        	}

        	m_lowWaveFieldWriter->flush();
        }

        // Update last time step
        m_lastTimeStep = time;

        logInfo(rank) << "Writing wave field at time" << utils::nospace << time << ". Done.";
    }

    /**
     * Close wave field writer and free resources
     */
    void close()
    {
        if (!m_enabled)
            return;

		delete m_waveFieldWriter;
		m_waveFieldWriter = 0L;
		delete m_lowWaveFieldWriter;
		m_lowWaveFieldWriter = 0L;
		delete m_variableSubsampler;
		m_variableSubsampler = 0L;
		delete m_outputBuffer;
		m_outputBuffer = 0L;
    }
};

}

#endif // WAVE_FIELD_WRITER_H
