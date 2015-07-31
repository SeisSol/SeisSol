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

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <string>
#include <vector>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "Geometry/MeshReader.h"
#include "Geometry/refinement/RefinerUtils.h"
#include "Geometry/refinement/Refinement.h"
#include "Geometry/refinement/NoRefinement.h"
#include "Geometry/refinement/StaticRefinement.h"

namespace seissol
{

class WaveFieldWriter
{
private:
    /** True if wave field output is enabled */
    bool m_enabled;

    /** The rank of the process */
    int m_rank;

    /** The output prefix for the filename */
    std::string m_outputPrefix;

    /** The XMDF Writer used for the wave field */
    xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_waveFieldWriter;

    /** Number of variables */
    unsigned int m_numVariables;

    /** Tet refinement strategy */
    refinement::Refinement<double>* m_tetRefinement;

    /** Pointer to the degrees of freedom */
    const double* m_dofs;

    /** Mapping from the cell order to dofs order */
    const unsigned int* m_map;

    /** Buffer required to extract the output data from the unknowns */
    double *m_outputBuffer;

    /** Checks wether a refinement strategy is valid or not. **/
    static bool isRefinementStartegyValid(int refinement) {
        switch(refinement) {
            case(0): // No refinement strategy
                return true;
            case(1): // Strategy "Divide by 4"
                return true;
            case(2): // Strategy "Divide by 8"
                return true;
            case(3): // Strategy "Divide by 32"
                return true;
            default:
                return false;
        }
    }

public:
    WaveFieldWriter()
    : m_enabled(false), m_rank(0),
    m_waveFieldWriter(0L),
    m_numVariables(0),
    m_tetRefinement(0L),
    m_dofs(0L), m_map(0L),
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
     */
    void init(int numVars, int order,
            const MeshReader &meshReader,
            const double* dofs, const unsigned int* map,
            int refinement, int timestep)
    {
        if (!m_enabled)
            return;


#ifdef USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
#endif // USE_MPI

        logInfo(m_rank) << "Initializing HDF5 wave field output.";

        if (m_waveFieldWriter != 0L)
            logError() << "Wave field writer already initialized";

        m_numVariables = numVars;

        if (numVars != 9)
            logError() << "XDMF output supports currently only 9 variables. Number of variables specified:" << m_numVariables;
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

        m_waveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
                m_rank, m_outputPrefix.c_str(), variables, timestep);

        // Setup cell refinement class
        if (!isRefinementStartegyValid(refinement)) {
            logError()  << "Refinement Strategy is invalid!" << std::endl
            << "Current value : " << refinement << std::endl
            << "Valid options are :" << std::endl
            << "0 - No Refinement" << std::endl
            << "1 - Refinement by 4" << std::endl
            << "2 - Refinement by 8" << std::endl
            << "3 - Refinement by 32";
        }
        if (refinement == 0) {
            logInfo(m_rank) << "Refinement is turned off.";
            m_tetRefinement = new refinement::NoRefinement<double>(
                    meshReader, order, numVars);
        }
        else {
            switch(refinement) {
                case(1):
                    logInfo(m_rank) << "Refinement Startegy is \"Divide by 4\"";
                    m_tetRefinement = new refinement::StaticRefinement<double>(
                            meshReader,
                            refinement::DivideTetrahedronBy4<double>(),
                            order, numVars);
                    break;
                case(2):
                    logInfo(m_rank) << "Refinement Startegy is \"Divide by 8\"";
                    m_tetRefinement = new refinement::StaticRefinement<double>(
                            meshReader,
                            refinement::DivideTetrahedronBy8<double>(),
                            order, numVars);
                    break;
                case(3):
                    logInfo(m_rank) << "Refinement Startegy is \"Divide by 32\"";
                    m_tetRefinement = new refinement::StaticRefinement<double>(
                            meshReader,
                            refinement::DivideTetrahedronBy32<double>(),
                            order, numVars);
                    break;
            }
        }

        logInfo(m_rank) << "Refinement class initialized : ("
        << meshReader.getElements().size() << "Cells, "
        << meshReader.getVertices().size() << "Vertices)"
        << " -refined-to-> (" 
        << m_tetRefinement->nCells() << "Cells, " 
        << m_tetRefinement->nVertices() << " Vertices)";

        m_waveFieldWriter->init(
                m_tetRefinement->nCells(), m_tetRefinement->cells(),
                m_tetRefinement->nVertices(), m_tetRefinement->vertices(),
                true);

        logInfo(m_rank) << "WaveFieldWriter initialized";

        // Create output buffer
        m_outputBuffer = new double[m_tetRefinement->nCells()];

        // Save dof/map pointer
        m_dofs = dofs;
        m_map = map;

        logInfo(m_rank) << "Initializing HDF5 wave field output. Done.";
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

        logInfo(m_rank) << "Writing wave field at time" << utils::nospace << time << '.';

        m_waveFieldWriter->addTimeStep(time);

        for (unsigned int i = 0; i < m_numVariables; i++) {
            m_tetRefinement->get(m_dofs, m_map, i, m_outputBuffer);

            m_waveFieldWriter->writeData(i, m_outputBuffer);
        }

        m_waveFieldWriter->flush();

        logInfo(m_rank) << "Writing wave field at time" << utils::nospace << time << ". Done.";
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
        delete m_tetRefinement;
        m_tetRefinement = 0L;
        delete m_outputBuffer;
        m_outputBuffer = 0L;
    }
};

}

#endif // WAVE_FIELD_WRITER_H
