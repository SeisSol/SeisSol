// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_

#include <mpi.h>

#include "Kernels/Precision.h"
#include "Monitoring/Stopwatch.h"
#include "async/ExecInfo.h"
#include "xdmfwriter/XdmfWriter.h"

namespace seissol::writer {

struct FaultInitParam {
  static const unsigned int OutputMaskSize = 20;

  bool outputMask[OutputMaskSize];
  int timestep;
  xdmfwriter::BackendType backend;
  std::string backupTimeStamp;
};

struct FaultParam {
  double time;
};

class FaultWriterExecutor {
  public:
  enum BufferIds {
    OutputPrefix = 0,
    Cells = 1,
    Vertices = 2,
    FaultTags = 3,
    GlobalIds = 4,
    Variables0 = 5
  };

  private:
  xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* m_xdmfWriter{nullptr};

  /** The MPI communicator for the writer */
  MPI_Comm m_comm;

  /** The number of variables that should be written */
  unsigned int m_numVariables{0};

  /** Backend stopwatch */
  Stopwatch m_stopwatch;

  bool m_enabled{false};

  public:
  FaultWriterExecutor()
      : m_comm(MPI_COMM_NULL)

  {}

  /**
   * Initialize the XDMF writer
   */
  void execInit(const async::ExecInfo& info, const FaultInitParam& param);

  void exec(const async::ExecInfo& info, const FaultParam& param) {
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

  void setFaultTagsData(const unsigned int* faultTags) {
    m_xdmfWriter->writeExtraIntCellData(0, faultTags);
  }

  void finalize() {
    if (m_enabled) {
      // note: also includes some ranks which do nothing at all
      m_stopwatch.printTime("Time fault writer backend:");
    }

    if (m_comm != MPI_COMM_NULL) {
      MPI_Comm_free(&m_comm);
      m_comm = MPI_COMM_NULL;
    }

    delete m_xdmfWriter;
    m_xdmfWriter = nullptr;
  }

  static std::string getLabelName(size_t index) { return Labels[index]; }

  private:
  /** Variable names in the output */
  static const char* const Labels[];
};

} // namespace seissol::writer

// namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_
