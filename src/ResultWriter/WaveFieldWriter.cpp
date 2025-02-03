// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include <Geometry/MeshDefinition.h>
#include <Geometry/Refinement/RefinerUtils.h>
#include <Geometry/Refinement/VariableSubSampler.h>
#include <Initializer/Parameters/OutputParameters.h>
#include <Kernels/Precision.h>
#include <ResultWriter/WaveFieldWriterExecutor.h>
#include <algorithm>
#include <async/Module.h>
#include <cassert>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <mpi.h>
#include <ostream>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

#include "Geometry/MeshReader.h"
#include "Geometry/Refinement/MeshRefiner.h"
#include "Modules/Modules.h"
#include "SeisSol.h"
#include "WaveFieldWriter.h"

void seissol::writer::WaveFieldWriter::setUp() {
  setExecutor(m_executor);
  if (isAffinityNecessary()) {
    const auto freeCpus = seissolInstance.getPinning().getFreeCPUsMask();
    logInfo() << "Wave field writer thread affinity:" << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
    setAffinityIfNecessary(freeCpus);
  }
}

void seissol::writer::WaveFieldWriter::enable() { m_enabled = true; }

seissol::refinement::TetrahedronRefiner<double>*
    seissol::writer::WaveFieldWriter::createRefiner(int refinement) {
  refinement::TetrahedronRefiner<double>* tetRefiner = nullptr;
  switch (refinement) {
  case 0:
    logInfo() << "Refinement is turned off.";
    tetRefiner = new refinement::IdentityRefiner<double>();
    break;
  case 1:
    logInfo() << "Refinement Strategy is \"Divide by 4\"";
    tetRefiner = new refinement::DivideTetrahedronBy4<double>();
    break;
  case 2:
    logInfo() << "Refinement Strategy is \"Divide by 8\"";
    tetRefiner = new refinement::DivideTetrahedronBy8<double>();
    break;
  case 3:
    logInfo() << "Refinement Strategy is \"Divide by 32\"";
    tetRefiner = new refinement::DivideTetrahedronBy32<double>();
    break;
  default:
    logError() << "Refinement Strategy is invalid!" << std::endl
               << "Current value : " << refinement << std::endl
               << "Valid options are :" << std::endl
               << "0 - No Refinement" << std::endl
               << "1 - Refinement by 4" << std::endl
               << "2 - Refinement by 8" << std::endl
               << "3 - Refinement by 32";
  }
  return tetRefiner;
}

const unsigned*
    seissol::writer::WaveFieldWriter::adjustOffsets(refinement::MeshRefiner<double>* meshRefiner) {
  const unsigned* constCells = nullptr;
// Cells are a bit complicated because the vertex filter will now longer work if we just use the
// buffer We will add the offset later
#ifdef USE_MPI
  // Add the offset to the cells
  MPI_Comm groupComm = seissolInstance.asyncIO().groupComm();
  unsigned int offset = meshRefiner->getNumVertices();
  MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
  offset -= meshRefiner->getNumVertices();

  // Add the offset to all cells
  auto* cells = new unsigned int[meshRefiner->getNumCells() * 4];
  for (unsigned int i = 0; i < meshRefiner->getNumCells() * 4; i++) {
    cells[i] = meshRefiner->getCellData()[i] + offset;
  }
  constCells = cells;
#else  // USE_MPI
  const_cells = meshRefiner->getCellData();
#endif // USE_MPI
  return constCells;
}

std::vector<unsigned int> seissol::writer::WaveFieldWriter::generateRefinedClusteringData(
    refinement::MeshRefiner<double>* meshRefiner,
    const std::vector<unsigned>& ltsClusteringData,
    std::map<int, int>& newToOldCellMap) const {
  // subsampling preserves the order of the cells, so we just need to repeat the LtsClusteringData
  // kSubCellsPerCell times. we also need to account for the extractRegion filter via the
  // newToOldCellMap hash map
  std::vector<unsigned int> refinedClusteringData(meshRefiner->getNumCells());

  auto kSubCellsPerCell = static_cast<size_t>(meshRefiner->getkSubCellsPerCell());
  for (size_t j = 0; j < meshRefiner->getNumCells(); j++) {
    if (isExtractRegionEnabled) {
      refinedClusteringData[j] = ltsClusteringData[newToOldCellMap[(j / kSubCellsPerCell)]];
    } else {
      refinedClusteringData[j] = ltsClusteringData[(j / kSubCellsPerCell)];
    }
  }
  return refinedClusteringData;
}

void seissol::writer::WaveFieldWriter::init(
    unsigned int numVars,
    int order,
    int numAlignedDOF,
    const seissol::geometry::MeshReader& meshReader,
    const std::vector<unsigned>& ltsClusteringData,
    const real* dofs,
    const real* pstrain,
    const real* integrals,
    const unsigned int* map,
    const seissol::initializer::parameters::WaveFieldOutputParameters& parameters,
    xdmfwriter::BackendType backend,
    const std::string& backupTimeStamp) {
  if (!m_enabled) {
    return;
  }

  // Initialize the asynchronous module
  async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>::init();

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);

  logInfo() << "Initializing XDMF wave field output.";

  /** All initialization parameters */
  WaveFieldInitParam param{};

  param.timestep = 0;

  /** List of all buffer ids */
  param.bufferIds[OutputPrefix] =
      addSyncBuffer(m_outputPrefix.c_str(), m_outputPrefix.size() + 1, true);

  param.backend = backend;
  param.backupTimeStamp = backupTimeStamp;

  //
  // High order I/O
  //
  m_numVariables = numVars + WaveFieldWriterExecutor::NumPlasticityVariables;
  m_outputFlags = new bool[m_numVariables];
  for (size_t i = 0; i < numVars; i++) {
    m_outputFlags[i] = (parameters.outputMask[i]);
  }
  for (size_t i = 0; i < WaveFieldWriterExecutor::NumPlasticityVariables; i++) {
    m_outputFlags[numVars + i] = (pstrain != nullptr) && (parameters.plasticityMask[i]);
  }

  // WARNING: The m_outputFlags memory might be directly used by the executor.
  // Do not modify this array after the following line
  param.bufferIds[OutputFlags] = addSyncBuffer(m_outputFlags, m_numVariables * sizeof(bool), true);

  // Setup the tetrahedron refinement strategy
  refinement::TetrahedronRefiner<double>* tetRefiner =
      createRefiner(static_cast<int>(parameters.refinement));

  unsigned int numElems = meshReader.getElements().size();
  unsigned int numVerts = meshReader.getVertices().size();

  // These variables are only filled if a region is to be extracted
  // Elements of the extracted region
  std::vector<const Element*> subElements;
  // The oldToNewVertexMap defines a map between old vertex index to
  // new vertex index. This is used to assign the vertex subset as well as
  // used in MeshRefiner since the elements would hold old index of vertices
  std::map<int, int> oldToNewVertexMap;
  std::map<int, int> newToOldCellMap;
  // Vertices of the extracted region
  std::vector<const Vertex*> subVertices;
  // Mesh refiner
  refinement::MeshRefiner<double>* meshRefiner = nullptr;

  // If at least one of the bounds is non-zero then extract.
  const bool isExtractBoxEnabled = parameters.bounds.enabled;

  // If at least one group is explicitly enabled, extract
  const bool isExtractGroupEnabled = !parameters.groups.empty();

  isExtractRegionEnabled = isExtractBoxEnabled || isExtractGroupEnabled;
  // isExtractRegionEnabled = true  : Extract region
  // isExtractRegionEnabled = false : Entire region

  if (isExtractRegionEnabled) {
    const std::vector<Element>& allElements = meshReader.getElements();
    const std::vector<Vertex>& allVertices = meshReader.getVertices();

    // m_map will store a new map from new cell index to dof index
    // the old map is contained in the "map" variable - which is a map from old
    //    cell index to dof index
    m_map.resize(numElems);

    // Extract elements based on the specified group or region
    for (size_t i = 0; i < numElems; i++) {
      const auto groupId = allElements[i].group;
      // Store the current number of elements to check if new was added
      const bool isInRegion =
          !isExtractBoxEnabled ||
          parameters.bounds.contains(allVertices[allElements[i].vertices[0]].coords[0],
                                     allVertices[allElements[i].vertices[0]].coords[1],
                                     allVertices[allElements[i].vertices[0]].coords[2]) ||
          parameters.bounds.contains(allVertices[allElements[i].vertices[1]].coords[0],
                                     allVertices[allElements[i].vertices[1]].coords[1],
                                     allVertices[allElements[i].vertices[1]].coords[2]) ||
          parameters.bounds.contains(allVertices[allElements[i].vertices[2]].coords[0],
                                     allVertices[allElements[i].vertices[2]].coords[1],
                                     allVertices[allElements[i].vertices[2]].coords[2]) ||
          parameters.bounds.contains(allVertices[allElements[i].vertices[3]].coords[0],
                                     allVertices[allElements[i].vertices[3]].coords[1],
                                     allVertices[allElements[i].vertices[3]].coords[2]);
      const bool isInGroup = !isExtractGroupEnabled || parameters.groups.count(groupId) > 0;
      if (isInRegion && isInGroup) {
        // Assign the new map
        const size_t iNew = subElements.size();
        m_map[iNew] = map[i];

        // Push the address of the element into the vector
        subElements.push_back(&(allElements[i]));
        newToOldCellMap.insert(std::pair<size_t, size_t>(iNew, i));

        // Push the vertices into the map which makes sure that the entries are unique
        oldToNewVertexMap.insert(
            std::pair<int, int>(allElements[i].vertices[0], oldToNewVertexMap.size()));
        oldToNewVertexMap.insert(
            std::pair<int, int>(allElements[i].vertices[1], oldToNewVertexMap.size()));
        oldToNewVertexMap.insert(
            std::pair<int, int>(allElements[i].vertices[2], oldToNewVertexMap.size()));
        oldToNewVertexMap.insert(
            std::pair<int, int>(allElements[i].vertices[3], oldToNewVertexMap.size()));
      }
    }

    subVertices.resize(oldToNewVertexMap.size());

    for (auto& it : oldToNewVertexMap) {
      subVertices.at(it.second) = &allVertices[it.first];
    }

    numElems = subElements.size();
    numVerts = subVertices.size();

    meshRefiner = new refinement::MeshRefiner<double>(
        subElements, subVertices, oldToNewVertexMap, *tetRefiner);
  } else {
    meshRefiner = new refinement::MeshRefiner<double>(meshReader, *tetRefiner);
    m_map.resize(numElems);
    std::copy_n(map, numElems, m_map.begin());
  }

  logInfo() << "Refinement class initialized";
  logDebug() << "Cells : " << numElems << "refined-to ->" << meshRefiner->getNumCells();
  logDebug() << "Vertices : " << numVerts << "refined-to ->" << meshRefiner->getNumVertices();

  m_variableSubsampler = std::make_unique<refinement::VariableSubsampler<double>>(
      numElems, *tetRefiner, order, numVars, numAlignedDOF);
  m_variableSubsamplerPStrain = std::make_unique<refinement::VariableSubsampler<double>>(
      numElems,
      *tetRefiner,
      order,
      static_cast<unsigned int>(WaveFieldWriterExecutor::NumPlasticityVariables),
      numAlignedDOF);

  logInfo() << "VariableSubsampler initialized";

  // Delete the tetRefiner since it is no longer required
  delete tetRefiner;

  const unsigned int* constCells = adjustOffsets(meshRefiner);

  // Create mesh buffers
  param.bufferIds[Cells] =
      addSyncBuffer(constCells, meshRefiner->getNumCells() * 4 * sizeof(unsigned int));
  param.bufferIds[Vertices] = addSyncBuffer(meshRefiner->getVertexData(),
                                            meshRefiner->getNumVertices() * 3 * sizeof(double));
  std::vector<unsigned int> refinedClusteringData =
      generateRefinedClusteringData(meshRefiner, ltsClusteringData, newToOldCellMap);
  param.bufferIds[Clustering] = addSyncBuffer(refinedClusteringData.data(),
                                              meshRefiner->getNumCells() * sizeof(unsigned int));

  // Create data buffers
  bool first = false;
  for (unsigned int i = 0; i < m_numVariables; i++) {
    if (m_outputFlags[i]) {
      const unsigned int id = addBuffer(nullptr, meshRefiner->getNumCells() * sizeof(real));
      if (!first) {
        param.bufferIds[Variables0] = id;
        first = true;
      }
    }
  }

  // Save number of cells
  m_numCells = meshRefiner->getNumCells();
  // Set up for low order output flags
  m_lowOutputFlags = new bool[WaveFieldWriterExecutor::NumLowvariables];
  m_numIntegratedVariables = seissolInstance.postProcessor().getNumberOfVariables();

  seissolInstance.postProcessor().getIntegrationMask(&m_lowOutputFlags[0]);
  param.bufferIds[LowOutputFlags] = addSyncBuffer(
      m_lowOutputFlags, WaveFieldWriterExecutor::NumLowvariables * sizeof(bool), true);
  //
  //  Low order I/O
  //
  refinement::MeshRefiner<double>* pLowMeshRefiner = nullptr;
  const unsigned int* constLowCells = nullptr;
  if (integrals != nullptr) {
    logInfo() << "Initialize low order output";

    // Refinement strategy (no refinement)
    const refinement::IdentityRefiner<double> lowTetRefiner;

    // Mesh refiner
    if (isExtractRegionEnabled) {
      pLowMeshRefiner = new refinement::MeshRefiner<double>(
          subElements, subVertices, oldToNewVertexMap, lowTetRefiner);
    } else {
      pLowMeshRefiner = new refinement::MeshRefiner<double>(meshReader, lowTetRefiner);
    }

    constLowCells = adjustOffsets(pLowMeshRefiner);

    // Create mesh buffers
    param.bufferIds[LowCells] =
        addSyncBuffer(constLowCells, pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
    param.bufferIds[LowVertices] = addSyncBuffer(
        pLowMeshRefiner->getVertexData(), pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));

    // Create data buffers
    param.bufferIds[LowVariables0] =
        addBuffer(nullptr, pLowMeshRefiner->getNumCells() * sizeof(real));

    const int numLowVars = m_numIntegratedVariables;

    for (int i = 1; i < numLowVars; i++) {
      addBuffer(nullptr, pLowMeshRefiner->getNumCells() * sizeof(real));
    }

    // Save number of cells
    m_numLowCells = pLowMeshRefiner->getNumCells();
  } else {
    // No low order output
    param.bufferIds[LowCells] = -1;
    param.bufferIds[LowVertices] = -1;
    param.bufferIds[LowVariables0] = -1;
  }

  //
  // Send all buffers for initialization
  //
  sendBuffer(param.bufferIds[OutputPrefix], m_outputPrefix.size() + 1);

  sendBuffer(param.bufferIds[OutputFlags], m_numVariables * sizeof(bool));

  sendBuffer(param.bufferIds[Cells], meshRefiner->getNumCells() * 4 * sizeof(unsigned int));
  sendBuffer(param.bufferIds[Vertices], meshRefiner->getNumVertices() * 3 * sizeof(double));
  sendBuffer(param.bufferIds[Clustering], meshRefiner->getNumCells() * sizeof(unsigned int));

  if (integrals != nullptr) {
    sendBuffer(param.bufferIds[LowCells],
               pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
    sendBuffer(param.bufferIds[LowVertices],
               pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));
    sendBuffer(param.bufferIds[LowOutputFlags],
               WaveFieldWriterExecutor::NumLowvariables * sizeof(bool));
  }

  // Initialize the executor
  callInit(param);

  // Remove buffers
  removeBuffer(param.bufferIds[OutputPrefix]);
  removeBuffer(param.bufferIds[Cells]);
  removeBuffer(param.bufferIds[Vertices]);
  removeBuffer(param.bufferIds[Clustering]);
  if (integrals != nullptr) {
    removeBuffer(param.bufferIds[LowCells]);
    removeBuffer(param.bufferIds[LowVertices]);
  }

  // Remove the low mesh refiner if it was setup
  delete pLowMeshRefiner;

#ifdef USE_MPI
  delete[] constCells;
  delete[] constLowCells;
#endif // USE_MPI

  // Save dof/map pointer
  m_dofs = dofs;
  m_pstrain = pstrain;
  m_integrals = integrals;

  m_variableBufferIds[0] = param.bufferIds[Variables0];
  m_variableBufferIds[1] = param.bufferIds[LowVariables0];

  delete meshRefiner;
}

void seissol::writer::WaveFieldWriter::write(double time) {
  SCOREP_USER_REGION("WaveFieldWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION);

  if (!m_enabled) {
    return;
  }

  m_stopwatch.start();

  SCOREP_USER_REGION_DEFINE(r_wait);
  SCOREP_USER_REGION_BEGIN(r_wait, "wavfieldwriter_wait", SCOREP_USER_REGION_TYPE_COMMON);
  logInfo() << "Waiting for last wave field.";
  wait();
  SCOREP_USER_REGION_END(r_wait);

  logInfo() << "Writing wave field at time" << utils::nospace << time << '.';

  unsigned int nextId = m_variableBufferIds[0];
  for (unsigned int i = 0; i < m_numVariables; i++) {
    if (!m_outputFlags[i]) {
      continue;
    }

    real* managedBuffer =
        async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>::managedBuffer<
            real*>(nextId);
    if (i < m_numVariables - WaveFieldWriterExecutor::NumPlasticityVariables) {
      m_variableSubsampler->get(m_dofs, m_map.data(), i, managedBuffer);
    } else {
      m_variableSubsamplerPStrain->get(
          m_pstrain,
          m_map.data(),
          i - (m_numVariables - WaveFieldWriterExecutor::NumPlasticityVariables),
          managedBuffer);
    }
    for (unsigned int j = 0; j < m_numCells; j++) {
      if (!std::isfinite(managedBuffer[j])) {
        logError() << "Detected Inf/NaN in volume output. Aborting.";
      }
    }
    sendBuffer(nextId, m_numCells * sizeof(real));

    nextId++;
  }

  // nextId is required in a manner similar to above for writing integrated variables
  nextId = 0;

  if (m_integrals != nullptr) {
    for (unsigned int i = 0; i < WaveFieldWriterExecutor::NumIntegratedVariables; i++) {
      if (!m_lowOutputFlags[i]) {
        continue;
      }
      real* managedBuffer =
          async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>::managedBuffer<
              real*>(m_variableBufferIds[1] + nextId);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif // _OPENMP
      for (unsigned int j = 0; j < m_numLowCells; j++) {
        managedBuffer[j] = m_integrals[m_map[j] * m_numIntegratedVariables + nextId];
      }

      sendBuffer(m_variableBufferIds[1] + nextId, m_numLowCells * sizeof(real));
      nextId++;
    }
  }

  WaveFieldParam param;
  param.time = time;
  call(param);

  m_stopwatch.pause();

  logInfo() << "Writing wave field at time" << utils::nospace << time << ". Done.";
}

void seissol::writer::WaveFieldWriter::simulationStart() { syncPoint(0.0); }

void seissol::writer::WaveFieldWriter::syncPoint(double currentTime) { write(currentTime); }
