// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_

#include "Kernels/Precision.h"
#include "Monitoring/Stopwatch.h"
#include "xdmfwriter/XdmfWriter.h"

#include <async/ExecInfo.h>

namespace seissol::writer {
struct FreeSurfaceInitParam {
  int timestep{};
  xdmfwriter::BackendType backend{};
  std::string backupTimeStamp;
};

struct FreeSurfaceParam {
  double time{};
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
  MPI_Comm comm_{MPI_COMM_NULL};

  xdmfwriter::XdmfWriter<xdmfwriter::TRIANGLE, double, real>* xdmfWriter_{nullptr};
  unsigned numVariables_{0};

  /** Backend stopwatch */
  Stopwatch stopwatch_;

  bool enabled_{false};

  public:
  FreeSurfaceWriterExecutor() = default;

  /**
   * Initialize the XDMF writer
   */
  void execInit(const async::ExecInfo& info, const FreeSurfaceInitParam& param);

  void exec(const async::ExecInfo& info, const FreeSurfaceParam& param) {
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

  void setLocationFlagData(const unsigned int* locationFlags) {
    xdmfWriter_->writeExtraIntCellData(0, locationFlags);
  }

  void finalize() {
    if (enabled_) {
      // note: also includes some ranks which do nothing at all
      stopwatch_.printTime("Time free surface writer backend:");
    }

    if (comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&comm_);
      comm_ = MPI_COMM_NULL;
    }

    delete xdmfWriter_;
    xdmfWriter_ = nullptr;
  }

  private:
  /** Variable names in the output */
  static const char* const Labels[];
};

} // namespace seissol::writer

// namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITEREXECUTOR_H_
