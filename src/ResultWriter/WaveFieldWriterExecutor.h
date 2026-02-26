// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITEREXECUTOR_H_
#define SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITEREXECUTOR_H_

#include "Equations/Datastructures.h"
#include "Kernels/Precision.h"
#include "Model/Plasticity.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/MPI.h"
#include "xdmfwriter/XdmfWriter.h"

#include <async/ExecInfo.h>
#include <cassert>
#include <memory>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::writer {

/** Buffer tags for asynchronous IO */
enum BufferTags {
  OutputPrefix = 0,
  OutputFlags = 1,
  Cells = 2,
  Vertices = 3,
  Clustering = 4,
  GlobalIds = 5,
  Variables0 = 6,
  LowCells = 7,
  LowVertices = 8,
  LowOutputFlags = 9,
  LowVariables0 = 10,
  BuffertagMax = LowVariables0
};

struct WaveFieldInitParam {
  int timestep{};
  int bufferIds[BuffertagMax + 1]{};
  xdmfwriter::BackendType backend{};
  std::string backupTimeStamp;
};

struct WaveFieldParam {
  double time{};
};

class WaveFieldWriterExecutor {
  private:
  /** The XMDF Writer used for the wave field */
  xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>* waveFieldWriter_{nullptr};

  /** The XDMF Writer for low order data */
  xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>* lowWaveFieldWriter_{nullptr};

  /** Buffer id for the first variable for high and low order output */
  unsigned int variableBufferIds_[2]{};

  /** The total number of (high order) variables */
  unsigned int numVariables_{0};

  /** Flag indicated which variables should be written */
  const bool* outputFlags_{nullptr};

  /** Flags indicating which low order variables should be written */
  const bool* lowOutputFlags_{nullptr};

  /** The MPI communicator for the XDMF writer */
  MPI_Comm comm_{MPI_COMM_NULL};

  /** Stopwatch for the wave field backend */
  Stopwatch stopwatch_;

  std::shared_ptr<std::vector<std::string>> varNames_;
  std::shared_ptr<std::vector<std::string>> varNamesLowRes_;

  public:
  WaveFieldWriterExecutor() = default;

  /**
   * Initialize the XDMF writers
   */
  void execInit(const async::ExecInfo& info, const WaveFieldInitParam& param) {
    if (waveFieldWriter_ != nullptr) {
      logError() << "Wave field writer already initialized";
    }

    int rank = seissol::Mpi::mpi.rank();

    const xdmfwriter::BackendType type = param.backend;

    const char* outputPrefix = static_cast<const char*>(info.buffer(param.bufferIds[OutputPrefix]));

    //
    // High order I/O
    //
    numVariables_ = info.bufferSize(param.bufferIds[OutputFlags]) / sizeof(bool);
    outputFlags_ = static_cast<const bool*>(info.buffer(param.bufferIds[OutputFlags]));

    varNames_ = std::make_shared<std::vector<std::string>>();
    varNamesLowRes_ = std::make_shared<std::vector<std::string>>();
    for (const auto& quantity : seissol::model::MaterialT::Quantities) {
      varNames_->emplace_back(quantity);
      varNamesLowRes_->emplace_back("low_" + quantity);
    }
    for (const auto& quantity : seissol::model::PlasticityData::Quantities) {
      varNames_->emplace_back(quantity);
      varNamesLowRes_->emplace_back("low_" + quantity);
    }

    std::vector<const char*> variables;
    for (unsigned int i = 0; i < numVariables_; i++) {
      if (outputFlags_[i]) {
        assert(i < seissol::model::MaterialT::Quantities.size() +
                       seissol::model::PlasticityData::Quantities.size());
        variables.push_back(varNames_->at(i).c_str());
      }
    }

    // Split the communicator into two - those containing vertices and those
    //  not containing any vertices.
    const int commColour = (info.bufferSize(param.bufferIds[Cells]) == 0) ? 0 : 1;
    MPI_Comm_split(seissol::Mpi::mpi.comm(), commColour, rank, &comm_);
    // Start the if statement
    if (info.bufferSize(param.bufferIds[Cells]) != 0) {
      // Get the new rank
      MPI_Comm_rank(comm_, &rank);

      // Initialize the I/O handler and write the mesh
      waveFieldWriter_ = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>(
          type, outputPrefix, param.timestep);

      waveFieldWriter_->setComm(comm_);
      waveFieldWriter_->setBackupTimeStamp(param.backupTimeStamp);
      const std::string extraIntVarName = "clustering";
      const auto vertexFilter = utils::Env("").get<bool>("SEISSOL_VERTEXFILTER", true);
      waveFieldWriter_->init(variables,
                             std::vector<const char*>(),
                             {extraIntVarName, "global-id"},
                             vertexFilter,
                             true);
      waveFieldWriter_->setMesh(
          info.bufferSize(param.bufferIds[Cells]) / (4 * sizeof(unsigned int)),
          static_cast<const unsigned int*>(info.buffer(param.bufferIds[Cells])),
          info.bufferSize(param.bufferIds[Vertices]) / (3 * sizeof(double)),
          static_cast<const double*>(info.buffer(param.bufferIds[Vertices])),
          param.timestep != 0);

      setClusteringData(static_cast<const unsigned int*>(info.buffer(param.bufferIds[Clustering])));

      waveFieldWriter_->writeExtraIntCellData(
          1, static_cast<const unsigned int*>(info.buffer(param.bufferIds[GlobalIds])));
      logInfo() << "High order output initialized";

      //
      // Low order I/O
      //
      if (param.bufferIds[LowCells] >= 0) {
        // Pstrain or Integrated quantities enabled
        lowOutputFlags_ = static_cast<const bool*>(info.buffer(param.bufferIds[LowOutputFlags]));

        std::vector<const char*> lowVariables;
        for (size_t i = 0; i < NumLowvariables; i++) {
          if (lowOutputFlags_[i]) {
            lowVariables.push_back(varNamesLowRes_->at(i).c_str());
          }
        }

        lowWaveFieldWriter_ = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>(
            type, (std::string(outputPrefix) + "-low").c_str());

        lowWaveFieldWriter_->setComm(comm_);
        lowWaveFieldWriter_->setBackupTimeStamp(param.backupTimeStamp);

        lowWaveFieldWriter_->init(lowVariables, std::vector<const char*>(), {}, vertexFilter);
        lowWaveFieldWriter_->setMesh(
            info.bufferSize(param.bufferIds[LowCells]) / (4 * sizeof(unsigned int)),
            static_cast<const unsigned int*>(info.buffer(param.bufferIds[LowCells])),
            info.bufferSize(param.bufferIds[LowVertices]) / (3 * sizeof(double)),
            static_cast<const double*>(info.buffer(param.bufferIds[LowVertices])),
            param.timestep != 0);

        logInfo() << "Low order output initialized";
      }

      // Save ids for the variables
      variableBufferIds_[0] = param.bufferIds[Variables0];
      variableBufferIds_[1] = param.bufferIds[LowVariables0];

      logInfo() << "Initializing XDMF wave field output. Done.";
    }
    // End the if statement
  }

  void setClusteringData(const unsigned* clustering) {
    waveFieldWriter_->writeExtraIntCellData(0, clustering);
  }

  void exec(const async::ExecInfo& info, const WaveFieldParam& param) {
    // Execute this function only if waveFieldWriter_ is initialized
    if (waveFieldWriter_ != nullptr) {
      stopwatch_.start();

      // High order output
      waveFieldWriter_->addTimeStep(param.time);

      unsigned int nextId = 0;
      for (unsigned int i = 0; i < numVariables_; i++) {
        if (outputFlags_[i]) {
          waveFieldWriter_->writeCellData(
              nextId, static_cast<const real*>(info.buffer(variableBufferIds_[0] + nextId)));

          nextId++;
        }
      }

      waveFieldWriter_->flush();

      // Low order output
      if (lowWaveFieldWriter_ != nullptr) {
        lowWaveFieldWriter_->addTimeStep(param.time);

        nextId = 0;
        for (unsigned int i = 0; i < NumLowvariables; i++) {
          if (lowOutputFlags_[i]) {
            lowWaveFieldWriter_->writeCellData(
                nextId, static_cast<const real*>(info.buffer(variableBufferIds_[1] + nextId)));

            nextId++;
          }
        }

        lowWaveFieldWriter_->flush();
      }

      stopwatch_.pause();
    }
  }

  void finalize() {
    if (waveFieldWriter_ != nullptr) {
      stopwatch_.printTime("Time wave field writer backend:");
    }

    if (comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&comm_);
      comm_ = MPI_COMM_NULL;
    }

    delete waveFieldWriter_;
    waveFieldWriter_ = nullptr;
    delete lowWaveFieldWriter_;
    lowWaveFieldWriter_ = nullptr;
  }

  static constexpr unsigned int NumPlasticityVariables = 7;
  static constexpr unsigned int NumIntegratedVariables = 9;
  static constexpr unsigned int NumLowvariables = NumIntegratedVariables;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITEREXECUTOR_H_
