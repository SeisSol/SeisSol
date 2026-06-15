// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_

#include "Kernels/Precision.h"
#include "Monitoring/Stopwatch.h"
#include "xdmfwriter/XdmfWriter.h"

#include <async/ExecInfo.h>
#include <mpi.h>

namespace seissol::writer {

struct FaultInitParam {
  static const unsigned int OutputMaskSize = 20;

  bool outputMask[OutputMaskSize]{};
  int timestep{};
  xdmfwriter::BackendType backend{};
  std::string backupTimeStamp;
};

struct FaultParam {
  double time{};
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
  xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* xdmfWriter_{nullptr};

  /** The MPI communicator for the writer */
  MPI_Comm comm_;

  /** The number of variables that should be written */
  unsigned int numVariables_{0};

  /** Backend stopwatch */
  Stopwatch stopwatch_;

  bool enabled_{false};

  public:
  FaultWriterExecutor()
      : comm_(MPI_COMM_NULL)

  {}

  /**
   * Initialize the XDMF writer
   */
  void execInit(const async::ExecInfo& info, const FaultInitParam& param);

  void exec(const async::ExecInfo& info, const FaultParam& param) {
    if (xdmfWriter_ == nullptr) {
      return;
    }

    stopwatch_.start();

    xdmfWriter_->addTimeStep(param.time);

    for (unsigned int i = 0; i < numVariables_; i++) {
      xdmfWriter_->writeCellData(i, static_cast<const real*>(info.buffer(Variables0 + i)));
    }

    xdmfWriter_->flush();

    stopwatch_.pause();
  }

  void setFaultTagsData(const unsigned int* faultTags) {
    xdmfWriter_->writeExtraIntCellData(0, faultTags);
  }

  void finalize() {
    if (enabled_) {
      // note: also includes some ranks which do nothing at all
      stopwatch_.printTime("Time fault writer backend:");
    }

    if (comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&comm_);
      comm_ = MPI_COMM_NULL;
    }

    delete xdmfWriter_;
    xdmfWriter_ = nullptr;
  }

  static std::string getLabelName(size_t index) { return Labels[index]; }

  private:
  /** Variable names in the output */
  static const char* const Labels[];
};

} // namespace seissol::writer

// namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_FAULTWRITEREXECUTOR_H_
