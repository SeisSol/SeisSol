/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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

#ifndef WAVE_FIELD_WRITER_EXECUTOR_H
#define WAVE_FIELD_WRITER_EXECUTOR_H

#include "Parallel/MPI.h"

#include <Kernels/Precision.h>
#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "async/ExecInfo.h"

#include "Monitoring/Stopwatch.h"

#include "Equations/Datastructures.h"
#include "Model/Plasticity.h"

namespace seissol::writer {

/** Buffer tags for asynchronous IO */
enum BufferTags {
  OutputPrefix = 0,
  OutputFlags = 1,
  Cells = 2,
  Vertices = 3,
  Clustering = 4,
  Variables0 = 5,
  LowCells = 6,
  LowVertices = 7,
  LowOutputFlags = 8,
  LowVariables0 = 9,
  BuffertagMax = LowVariables0
};

struct WaveFieldInitParam {
  int timestep;
  int bufferIds[BuffertagMax + 1];
  xdmfwriter::BackendType backend;
  std::string backupTimeStamp;
};

struct WaveFieldParam {
  double time;
};

class WaveFieldWriterExecutor {
  private:
  /** The XMDF Writer used for the wave field */
  xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>* m_waveFieldWriter{nullptr};

  /** The XDMF Writer for low order data */
  xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>* m_lowWaveFieldWriter{nullptr};

  /** Buffer id for the first variable for high and low order output */
  unsigned int m_variableBufferIds[2];

  /** The total number of (high order) variables */
  unsigned int m_numVariables{0};

  /** Flag indicated which variables should be written */
  const bool* m_outputFlags{nullptr};

  /** Flags indicating which low order variables should be written */
  const bool* m_lowOutputFlags{nullptr};

#ifdef USE_MPI
  /** The MPI communicator for the XDMF writer */
  MPI_Comm m_comm{MPI_COMM_NULL};
#endif // USE_MPI

  /** Stopwatch for the wave field backend */
  Stopwatch m_stopwatch;

  std::shared_ptr<std::vector<std::string>> varNames;
  std::shared_ptr<std::vector<std::string>> varNamesLowRes;

  public:
  WaveFieldWriterExecutor() = default;

  /**
   * Initialize the XDMF writers
   */
  void execInit(const async::ExecInfo& info, const WaveFieldInitParam& param) {
    if (m_waveFieldWriter != nullptr) {
      logError() << "Wave field writer already initialized";
    }

    int rank = seissol::MPI::mpi.rank();

    const xdmfwriter::BackendType type = param.backend;

    const char* outputPrefix = static_cast<const char*>(info.buffer(param.bufferIds[OutputPrefix]));

    //
    // High order I/O
    //
    m_numVariables = info.bufferSize(param.bufferIds[OutputFlags]) / sizeof(bool);
    m_outputFlags = static_cast<const bool*>(info.buffer(param.bufferIds[OutputFlags]));

    varNames = std::make_shared<std::vector<std::string>>();
    varNamesLowRes = std::make_shared<std::vector<std::string>>();
    for (const auto& quantity : seissol::model::MaterialT::Quantities) {
      varNames->emplace_back(quantity);
      varNamesLowRes->emplace_back("low_" + quantity);
    }
    for (const auto& quantity : seissol::model::PlasticityData::Quantities) {
      varNames->emplace_back(quantity);
      varNamesLowRes->emplace_back("low_" + quantity);
    }

    std::vector<const char*> variables;
    for (unsigned int i = 0; i < m_numVariables; i++) {
      if (m_outputFlags[i]) {
#ifdef USE_POROELASTIC
        assert(i < 20);
#else
        assert(i < 16);
#endif
        variables.push_back(varNames->at(i).c_str());
      }
    }

#ifdef USE_MPI
    // Split the communicator into two - those containing vertices and those
    //  not containing any vertices.
    const int commColour = (info.bufferSize(param.bufferIds[Cells]) == 0) ? 0 : 1;
    MPI_Comm_split(seissol::MPI::mpi.comm(), commColour, rank, &m_comm);
    // Start the if statement
    if (info.bufferSize(param.bufferIds[Cells]) != 0) {
      // Get the new rank
      MPI_Comm_rank(m_comm, &rank);
#endif // USE_MPI

      // Initialize the I/O handler and write the mesh
      m_waveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>(
          type, outputPrefix, param.timestep);

#ifdef USE_MPI
      m_waveFieldWriter->setComm(m_comm);
#endif // USE_MPI
      m_waveFieldWriter->setBackupTimeStamp(param.backupTimeStamp);
      const std::string extraIntVarName = "clustering";

      m_waveFieldWriter->init(
          variables, std::vector<const char*>(), extraIntVarName.c_str(), true, true);
      m_waveFieldWriter->setMesh(
          info.bufferSize(param.bufferIds[Cells]) / (4 * sizeof(unsigned int)),
          static_cast<const unsigned int*>(info.buffer(param.bufferIds[Cells])),
          info.bufferSize(param.bufferIds[Vertices]) / (3 * sizeof(double)),
          static_cast<const double*>(info.buffer(param.bufferIds[Vertices])),
          param.timestep != 0);

      setClusteringData(static_cast<const unsigned int*>(info.buffer(param.bufferIds[Clustering])));
      logInfo(rank) << "High order output initialized";

      //
      // Low order I/O
      //
      if (param.bufferIds[LowCells] >= 0) {
        // Pstrain or Integrated quantities enabled
        m_lowOutputFlags = static_cast<const bool*>(info.buffer(param.bufferIds[LowOutputFlags]));

        std::vector<const char*> lowVariables;
        for (size_t i = 0; i < NumLowvariables; i++) {
          if (m_lowOutputFlags[i]) {
            lowVariables.push_back(varNamesLowRes->at(i).c_str());
          }
        }

        m_lowWaveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>(
            type, (std::string(outputPrefix) + "-low").c_str());

#ifdef USE_MPI
        m_lowWaveFieldWriter->setComm(m_comm);
#endif // USE_MPI
        m_lowWaveFieldWriter->setBackupTimeStamp(param.backupTimeStamp);

        m_lowWaveFieldWriter->init(lowVariables, std::vector<const char*>());
        m_lowWaveFieldWriter->setMesh(
            info.bufferSize(param.bufferIds[LowCells]) / (4 * sizeof(unsigned int)),
            static_cast<const unsigned int*>(info.buffer(param.bufferIds[LowCells])),
            info.bufferSize(param.bufferIds[LowVertices]) / (3 * sizeof(double)),
            static_cast<const double*>(info.buffer(param.bufferIds[LowVertices])),
            param.timestep != 0);

        logInfo(rank) << "Low order output initialized";
      }

      // Save ids for the variables
      m_variableBufferIds[0] = param.bufferIds[Variables0];
      m_variableBufferIds[1] = param.bufferIds[LowVariables0];

      logInfo(rank) << "Initializing XDMF wave field output. Done.";
#ifdef USE_MPI
    }
    // End the if statement
#endif // USE_MPI
  }

  void setClusteringData(const unsigned* clustering) {
    m_waveFieldWriter->writeExtraIntCellData(clustering);
  }

  void exec(const async::ExecInfo& info, const WaveFieldParam& param) {
#ifdef USE_MPI
    // Execute this function only if m_waveFieldWriter is initialized
    if (m_waveFieldWriter != nullptr) {
#endif // USE_MPI
      m_stopwatch.start();

      // High order output
      m_waveFieldWriter->addTimeStep(param.time);

      unsigned int nextId = 0;
      for (unsigned int i = 0; i < m_numVariables; i++) {
        if (m_outputFlags[i]) {
          m_waveFieldWriter->writeCellData(
              nextId, static_cast<const real*>(info.buffer(m_variableBufferIds[0] + nextId)));

          nextId++;
        }
      }

      m_waveFieldWriter->flush();

      // Low order output
      if (m_lowWaveFieldWriter != nullptr) {
        m_lowWaveFieldWriter->addTimeStep(param.time);

        nextId = 0;
        for (unsigned int i = 0; i < NumLowvariables; i++) {
          if (m_lowOutputFlags[i]) {
            m_lowWaveFieldWriter->writeCellData(
                nextId, static_cast<const real*>(info.buffer(m_variableBufferIds[1] + nextId)));

            nextId++;
          }
        }

        m_lowWaveFieldWriter->flush();
      }

      m_stopwatch.pause();
#ifdef USE_MPI
    }
#endif // USE_MPI
  }

  void finalize() {
    if (m_waveFieldWriter != nullptr) {
      m_stopwatch.printTime("Time wave field writer backend:"
#ifdef USE_MPI
                            ,
                            m_comm
#endif // USE_MPI
      );
    }

#ifdef USE_MPI
    if (m_comm != MPI_COMM_NULL) {
      MPI_Comm_free(&m_comm);
      m_comm = MPI_COMM_NULL;
    }
#endif // USE_MPI

    delete m_waveFieldWriter;
    m_waveFieldWriter = nullptr;
    delete m_lowWaveFieldWriter;
    m_lowWaveFieldWriter = nullptr;
  }

  static constexpr unsigned int NumPlasticityVariables = 7;
  static constexpr unsigned int NumIntegratedVariables = 9;
  static constexpr unsigned int NumLowvariables = NumIntegratedVariables;
};

} // namespace seissol::writer

#endif // WAVE_FIELD_WRITER_EXECUTOR_H
