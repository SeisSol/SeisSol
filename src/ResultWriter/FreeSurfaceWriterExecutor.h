/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 */

#ifndef FREESURFACEWRITEREXECUTOR_H
#define FREESURFACEWRITEREXECUTOR_H

#include "async/ExecInfo.h"
#include "xdmfwriter/XdmfWriter.h"

#include "Monitoring/Stopwatch.h"
#include <Kernels/Precision.h>

namespace seissol::writer {
struct FreeSurfaceInitParam {
  int timestep;
  xdmfwriter::BackendType backend;
  std::string backupTimeStamp;
};

struct FreeSurfaceParam {
  double time;
};

class FreeSurfaceWriterExecutor {
  public:
  enum BufferIds {
    OutputPrefix = 0,
    Cells = 1,
    Vertices = 2,
    LocationFlags = 3,
    Variables0 = 4,
  };

  private:
#ifdef USE_MPI
  /** The MPI communicator for the writer */
  MPI_Comm m_comm{MPI_COMM_NULL};
#endif // USE_MPI

  xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* m_xdmfWriter{nullptr};
  unsigned m_numVariables{0};

  /** Backend stopwatch */
  Stopwatch m_stopwatch;

  public:
  FreeSurfaceWriterExecutor() = default;

  /**
   * Initialize the XDMF writer
   */
  void execInit(const async::ExecInfo& info, const FreeSurfaceInitParam& param);

  void exec(const async::ExecInfo& info, const FreeSurfaceParam& param) {
    if (m_xdmfWriter == nullptr) {
      return;
    }

    m_stopwatch.start();

    m_xdmfWriter->addTimeStep(param.time);

    for (unsigned int i = 0; i < m_numVariables; i++) {
      m_xdmfWriter->writeCellData(i, static_cast<const real*>(info.buffer(Variables0 + i)));
    }

    m_xdmfWriter->flush();

    m_stopwatch.pause();
  }

  void setLocationFlagData(const unsigned int* locationFlags) {
    m_xdmfWriter->writeExtraIntCellData(locationFlags);
  }

  void finalize() {
    if (m_xdmfWriter != nullptr) {
      m_stopwatch.printTime("Time free surface writer backend:"
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

    delete m_xdmfWriter;
    m_xdmfWriter = nullptr;
  }

  private:
  /** Variable names in the output */
  static const char* const Labels[];
};

} // namespace seissol::writer

// namespace seissol

#endif // FREESURFACEWRITEREXECUTOR_H
