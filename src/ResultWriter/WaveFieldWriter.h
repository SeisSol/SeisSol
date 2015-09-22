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
#include <memory>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "Geometry/MeshReader.h"
#include "Geometry/refinement/RefinerUtils.h"
#include "Geometry/refinement/MeshRefiner.h"
#include "Geometry/refinement/VariableSubSampler.h"

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
    std::auto_ptr<xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON> > m_waveFieldWriter;

    /** The variable subsampler for the refined mesh **/
    std::auto_ptr<refinement::VariableSubsampler<double> > m_variableSubsampler;

    /** Number of variables */
    unsigned int m_numVariables;

    /** Pointer to the degrees of freedom */
    const double* m_dofs;

    /** Mapping from the cell order to dofs order */
    const unsigned int* m_map;

    /** Buffer required to extract the output data from the unknowns */
    std::vector<double> m_outputBuffer;

    /** Number of cells after refinement step **/
    std::size_t m_refCellCount;

    std::vector<bool> m_output_flags;

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
    m_waveFieldWriter(NULL),
    m_numVariables(0),
    m_dofs(NULL), m_map(NULL),
    m_refCellCount(0)
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
    void init(int numVars, int order, int numAlignedDOF,
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

        if (m_waveFieldWriter.get() != NULL)
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

        // Currently all variables have to be chosen.
        m_output_flags.resize(numVars);
        std::fill(m_output_flags.begin(), m_output_flags.end(), true);

        m_waveFieldWriter.reset(new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
                m_rank, m_outputPrefix.c_str(), variables, timestep));

        if (!isRefinementStartegyValid(refinement)) {
            logError()  << "Refinement Strategy is invalid!" << std::endl
            << "Current value : " << refinement << std::endl
            << "Valid options are :" << std::endl
            << "0 - No Refinement" << std::endl
            << "1 - Refinement by 4" << std::endl
            << "2 - Refinement by 8" << std::endl
            << "3 - Refinement by 32";
        }

        std::auto_ptr<refinement::TetrahedronRefiner<double> > tetRefiner(NULL);
        switch(refinement) {
            case(0):
                logInfo(m_rank) << "Refinement is turned off.";
                tetRefiner.reset(new refinement::IdentityRefiner<double>());
                break;
            case(1):
                logInfo(m_rank) << "Refinement Startegy is \"Divide by 4\"";
                tetRefiner.reset(new refinement::DivideTetrahedronBy4<double>());
                break;
            case(2):
                logInfo(m_rank) << "Refinement Startegy is \"Divide by 8\"";
                tetRefiner.reset(new refinement::DivideTetrahedronBy8<double>());
                break;
            case(3):
                logInfo(m_rank) << "Refinement Startegy is \"Divide by 32\"";
                tetRefiner.reset(new refinement::DivideTetrahedronBy32<double>());
                break;
        }

        refinement::MeshRefiner<double> meshRefiner(meshReader, *tetRefiner);

        logInfo(m_rank) << "Refinement class initialized";
        logInfo(m_rank) << "Cells : "
        << meshReader.getElements().size() << " -refined-to-> "
        << meshRefiner.getNumCells();
        logInfo(m_rank) << "Vertices : "
        << meshReader.getVertices().size()
        << " -refined-to-> " 
        << meshRefiner.getNumVertices();

        m_waveFieldWriter->init(
                meshRefiner.getNumCells(), meshRefiner.getCellData(),
                meshRefiner.getNumVertices(), meshRefiner.getVertexData(),
                true);

        // Save cellcount for later use
        m_refCellCount = meshRefiner.getNumCells();

        logInfo(m_rank) << "WaveFieldWriter initialized";

        m_variableSubsampler.reset(new refinement::VariableSubsampler<double>(
                    *tetRefiner, order, numVars, numAlignedDOF
                    ));

        logInfo(m_rank) << "VariableSubsampler initialized";

        // Create output buffer
        m_outputBuffer.resize(meshRefiner.getNumCells());

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
            if (!m_output_flags[i]) continue;

            m_variableSubsampler->getSingle(m_dofs, m_map, i,
                    m_refCellCount, m_outputBuffer.data());

            m_waveFieldWriter->writeData(i, m_outputBuffer.data());
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

        m_waveFieldWriter.reset();
        m_variableSubsampler.reset();
        std::vector<double>().swap(m_outputBuffer); // reset and clean memory
    }
};

}

#endif // WAVE_FIELD_WRITER_H
