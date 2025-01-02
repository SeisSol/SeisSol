/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#include "FreeSurfaceWriter.h"

#include <Eigen/Dense>
#include <Geometry/MeshDefinition.h>
#include <Geometry/Refinement/TriangleRefiner.h>
#include <Initializer/PreProcessorMacros.h>
#include <Kernels/Precision.h>
#include <ResultWriter/FreeSurfaceWriterExecutor.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <async/Module.h>
#include <cassert>
#include <cstring>
#include <string>
#include <utils/logger.h>
#include <vector>

#include "AsyncCellIDs.h"
#include "Geometry/MeshTools.h"
#include "Modules/Modules.h"
#include "SeisSol.h"

void seissol::writer::FreeSurfaceWriter::constructSurfaceMesh(
    const seissol::geometry::MeshReader& meshReader,
    unsigned*& cells,
    double*& vertices,
    unsigned& nCells,
    unsigned& nVertices) {
  // TODO: Vertices could be pre-filtered
  nCells = m_freeSurfaceIntegrator->totalNumberOfTriangles;
  nVertices = 3 * m_freeSurfaceIntegrator->totalNumberOfTriangles;
  if (nCells == 0 || nVertices == 0) {
    cells = nullptr;
    vertices = nullptr;
    return;
  }

  cells = new unsigned[3 * nCells];
  vertices = new double[3 * nVertices];

  const std::vector<Element>& meshElements = meshReader.getElements();
  const std::vector<Vertex>& meshVertices = meshReader.getVertices();

  const unsigned numberOfSubTriangles = m_freeSurfaceIntegrator->triRefiner.subTris.size();

  unsigned idx = 0;
  unsigned* meshIds =
      m_freeSurfaceIntegrator->surfaceLtsTree.var(m_freeSurfaceIntegrator->surfaceLts.meshId);
  unsigned* sides =
      m_freeSurfaceIntegrator->surfaceLtsTree.var(m_freeSurfaceIntegrator->surfaceLts.side);
  for (unsigned fs = 0; fs < m_freeSurfaceIntegrator->totalNumberOfFreeSurfaces; ++fs) {
    const unsigned meshId = meshIds[fs];
    const unsigned side = sides[fs];
    Eigen::Vector3d x[3];
    Eigen::Vector3d a;
    Eigen::Vector3d b;
    for (unsigned vertex = 0; vertex < 3; ++vertex) {
      const unsigned tetVertex = MeshTools::FACE2NODES[side][vertex];
      const VrtxCoords& coords = meshVertices[meshElements[meshId].vertices[tetVertex]].coords;

      x[vertex](0) = coords[0];
      x[vertex](1) = coords[1];
      x[vertex](2) = coords[2];
    }
    a = x[1] - x[0];
    b = x[2] - x[0];

    for (unsigned tri = 0; tri < numberOfSubTriangles; ++tri) {
      const seissol::refinement::Triangle& subTri =
          m_freeSurfaceIntegrator->triRefiner.subTris[tri];
      for (unsigned vertex = 0; vertex < 3; ++vertex) {
        Eigen::Vector3d v = x[0] + subTri.x[vertex][0] * a + subTri.x[vertex][1] * b;
        vertices[3 * idx + 0] = v(0);
        vertices[3 * idx + 1] = v(1);
        vertices[3 * idx + 2] = v(2);
        cells[idx] = idx;
        ++idx;
      }
    }
  }
}

void seissol::writer::FreeSurfaceWriter::setUp() {
  setExecutor(m_executor);
  if (isAffinityNecessary()) {
    const auto freeCpus = seissolInstance.getPinning().getFreeCPUsMask();
    logInfo() << "Free surface writer thread affinity:"
              << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
  }
}

void seissol::writer::FreeSurfaceWriter::enable() { m_enabled = true; }

void seissol::writer::FreeSurfaceWriter::init(
    const seissol::geometry::MeshReader& meshReader,
    seissol::solver::FreeSurfaceIntegrator* freeSurfaceIntegrator,
    const char* outputPrefix,
    double interval,
    xdmfwriter::BackendType backend,
    const std::string& backupTimeStamp) {
  if (!m_enabled) {
    return;
  }

  m_freeSurfaceIntegrator = freeSurfaceIntegrator;

  logInfo() << "Initializing free surface output.";

  // Initialize the asynchronous module
  async::Module<FreeSurfaceWriterExecutor, FreeSurfaceInitParam, FreeSurfaceParam>::init();

  unsigned* cells = nullptr;
  double* vertices = nullptr;
  unsigned nCells = 0;
  unsigned nVertices = 0;
  constructSurfaceMesh(meshReader, cells, vertices, nCells, nVertices);

  const AsyncCellIDs<3> cellIds(nCells, nVertices, cells, seissolInstance);

  // Create buffer for output prefix
  unsigned int bufferId = addSyncBuffer(outputPrefix, strlen(outputPrefix) + 1, true);
  assert(bufferId == FreeSurfaceWriterExecutor::OutputPrefix);
  NDBG_UNUSED(bufferId);

  // Create mesh buffers
  bufferId = addSyncBuffer(cellIds.cells(), nCells * 3 * sizeof(unsigned));
  assert(bufferId == FreeSurfaceWriterExecutor::Cells);
  bufferId = addSyncBuffer(vertices, nVertices * 3 * sizeof(double));
  assert(bufferId == FreeSurfaceWriterExecutor::Vertices);
  bufferId =
      addSyncBuffer(m_freeSurfaceIntegrator->locationFlags.data(), nCells * sizeof(unsigned));
  assert(bufferId == FreeSurfaceWriterExecutor::LocationFlags);

  for (auto& velocity : m_freeSurfaceIntegrator->velocities) {
    addBuffer(velocity, nCells * sizeof(real));
  }
  for (auto& displacement : m_freeSurfaceIntegrator->displacements) {
    addBuffer(displacement, nCells * sizeof(real));
  }

  //
  // Send all buffers for initialization
  //
  sendBuffer(FreeSurfaceWriterExecutor::OutputPrefix);

  sendBuffer(FreeSurfaceWriterExecutor::Cells);
  sendBuffer(FreeSurfaceWriterExecutor::Vertices);
  sendBuffer(FreeSurfaceWriterExecutor::LocationFlags);

  // Initialize the executor
  FreeSurfaceInitParam param;
  param.timestep = 0;
  param.backend = backend;
  param.backupTimeStamp = backupTimeStamp;
  callInit(param);

  // Remove unused buffers
  removeBuffer(FreeSurfaceWriterExecutor::OutputPrefix);
  removeBuffer(FreeSurfaceWriterExecutor::Cells);
  removeBuffer(FreeSurfaceWriterExecutor::Vertices);
  removeBuffer(FreeSurfaceWriterExecutor::LocationFlags);

  // Register for the synchronization point hook
  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  setSyncInterval(interval);

  delete[] cells;
  delete[] vertices;
}

void seissol::writer::FreeSurfaceWriter::write(double time) {
  SCOREP_USER_REGION("FreeSurfaceWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (!m_enabled) {
    logError() << "Trying to write free surface output, but it is disabled.";
  }

  m_stopwatch.start();

  wait();

  logInfo() << "Writing free surface at time" << utils::nospace << time << ".";

  FreeSurfaceParam param;
  param.time = time;

  for (unsigned i = 0; i < 2 * FREESURFACE_NUMBER_OF_COMPONENTS; ++i) {
    sendBuffer(FreeSurfaceWriterExecutor::Variables0 + i);
  }

  call(param);

  m_stopwatch.pause();

  logInfo() << "Writing free surface at time" << utils::nospace << time << ". Done.";
}

void seissol::writer::FreeSurfaceWriter::simulationStart() { syncPoint(0.0); }

void seissol::writer::FreeSurfaceWriter::syncPoint(double currentTime) {
  SCOREP_USER_REGION("freesurfaceoutput", SCOREP_USER_REGION_TYPE_FUNCTION)

  m_freeSurfaceIntegrator->calculateOutput();
  write(currentTime);
}
