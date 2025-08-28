// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_

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
    GlobalIds = 4,
    Variables0 = 5,
  };

  private:
  /** The MPI communicator for the writer */
  MPI_Comm m_comm{MPI_COMM_NULL};

  xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* m_xdmfWriter{nullptr};
  unsigned m_numVariables{0};

  /** Backend stopwatch */
  Stopwatch m_stopwatch;

  bool m_enabled{false};

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
    m_xdmfWriter->writeExtraIntCellData(0, locationFlags);
  }

  void finalize() {
    if (m_enabled) {
      // note: also includes some ranks which do nothing at all
      m_stopwatch.printTime("Time free surface writer backend:");
    }

    if (m_comm != MPI_COMM_NULL) {
      MPI_Comm_free(&m_comm);
      m_comm = MPI_COMM_NULL;
    }

    delete m_xdmfWriter;
    m_xdmfWriter = nullptr;
  }

  private:
  /** Variable names in the output */
  static const char* const Labels[];
};

} // namespace seissol::writer

// namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_
