#include "Init.hpp"
#include "InitMesh.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>

#include "utils/logger.h"

#include "SeisSol.h"
#ifdef USE_NETCDF
#include "Geometry/NetcdfReader.h"
#endif // USE_NETCDF
#if defined(USE_HDF) && defined(USE_MPI)
#include "Geometry/PUMLReader.h"
#endif // defined(USE_HDF) && defined(USE_MPI)
#include "cubeGenerator.h"
#include "Modules/Modules.h"
#include "Monitoring/instrumentation.hpp"
#include "Monitoring/Stopwatch.h"
#include "Numerical_aux/Statistics.h"
#include "Initializer/time_stepping/LtsWeights/WeightsFactory.h"
#include "Solver/time_stepping/MiniSeisSol.h"
#include "ResultWriter/MiniSeisSolWriter.h"

#include "Parallel/MPI.h"

static void postMeshread(seissol::geometry::MeshReader& meshReader,
                         bool hasFault,
                         const std::array<double, 3>& displacement,
                         const std::array<std::array<double, 3>, 3>& scalingMatrix) {
  logInfo(seissol::MPI::mpi.rank()) << "The mesh has been read. Starting post processing.";
  
  if (meshReader.getElements().empty()) {
    logWarning(seissol::MPI::mpi.rank())
        << "There are no local mesh elements on this rank. Is your mesh big enough?";
  }

  meshReader.displaceMesh(displacement);
  meshReader.scaleMesh(scalingMatrix);

  if (hasFault) {
    logInfo(seissol::MPI::mpi.rank()) << "Extracting fault information.";

    auto* drParameters = seissol::SeisSol::main.getMemoryManager().getDRParameters();
    VrtxCoords center{drParameters->referencePoint[0],
                      drParameters->referencePoint[1],
                      drParameters->referencePoint[2]};
    meshReader.extractFaultInformation(center, drParameters->refPointMethod);
  }

  logInfo(seissol::MPI::mpi.rank()) << "Exchanging ghostlayer metadata.";
  meshReader.exchangeGhostlayerMetadata();

  seissol::SeisSol::main.getLtsLayout().setMesh(meshReader);
}

static void readMeshPUML(const seissol::initializer::parameters::SeisSolParameters& seissolParams) {
#if defined(USE_HDF) && defined(USE_MPI)
  const int rank = seissol::MPI::mpi.rank();
  double nodeWeight = 1.0;

#ifdef USE_MINI_SEISSOL
  if (seissol::MPI::mpi.size() > 1) {
    logInfo(rank) << "Running mini SeisSol to determine node weight";
    auto elapsedTime = seissol::miniSeisSol(seissol::SeisSol::main.getMemoryManager(),
                                            seissolParams.model.plasticity);
    nodeWeight = 1.0 / elapsedTime;

    const auto summary = seissol::statistics::parallelSummary(nodeWeight);
    logInfo(rank) << "Node weights: mean =" << summary.mean << " std =" << summary.std
                  << " min =" << summary.min << " median =" << summary.median
                  << " max =" << summary.max;

    writer::MiniSeisSolWriter writer(seissolParams.output.prefix.c_str());
    writer.write(elapsedTime, nodeWeight);
  }
#else
  logInfo(rank) << "Skipping mini SeisSol";
#endif

  logInfo(rank) << "Reading PUML mesh";

  seissol::Stopwatch watch;
  watch.start();

  bool readPartitionFromFile = seissol::SeisSol::main.simulator().checkPointingEnabled();

  using namespace seissol::initializers::time_stepping;
  LtsWeightsConfig config{seissolParams.model.materialFileName,
                          static_cast<unsigned int>(seissolParams.timeStepping.lts.rate),
                          seissolParams.timeStepping.vertexWeight.weightElement,
                          seissolParams.timeStepping.vertexWeight.weightDynamicRupture,
                          seissolParams.timeStepping.vertexWeight.weightFreeSurfaceWithGravity};

  const auto* ltsParameters = seissol::SeisSol::main.getMemoryManager().getLtsParameters();
  auto ltsWeights =
      getLtsWeightsImplementation(seissolParams.timeStepping.lts.weighttype, config, ltsParameters);
  auto meshReader =
      new seissol::geometry::PUMLReader(seissolParams.mesh.meshFileName.c_str(),
                                        seissolParams.mesh.partitioningLib.c_str(),
                                        seissolParams.timeStepping.maxTimestepWidth,
                                        seissolParams.output.checkpointParameters.fileName.c_str(),
                                        ltsWeights.get(),
                                        nodeWeight,
                                        readPartitionFromFile);
  seissol::SeisSol::main.setMeshReader(meshReader);

  watch.pause();
  watch.printTime("PUML mesh read in:");

#else // defined(USE_HDF) && defined(USE_MPI)
#ifndef USE_MPI
  logError() << "Tried to load a PUML mesh. However, PUML is currently only supported with MPI "
                "(and this build of SeisSol does not use MPI).";
#endif
#ifndef USE_HDF
  logError() << "Tried to load a PUML mesh. However, PUML needs SeisSol to be linked against HDF5.";
#endif
#endif // defined(USE_HDF) && defined(USE_MPI)
}

size_t getNumOutgoingEdges(seissol::geometry::MeshReader& meshReader) {
  auto& mpiNeighbors = meshReader.getMPINeighbors();
  size_t numEdges{0};
  for (auto& [_, neighborInfo] : mpiNeighbors) {
    // Note: this includes the case when multiple faces
    // of an element are located at a partition boarder
    numEdges += neighborInfo.elements.size();
  }
  return numEdges;
}

static void readCubeGenerator(const seissol::initializer::parameters::SeisSolParameters& seissolParams) {
  // unpack seissolParams
  const auto cubeGeneratorParameters = seissolParams.cubeGenerator;
  const auto boundaries = cubeGeneratorParameters.boundaries;
  const auto dims = cubeGeneratorParameters.dims;
  const auto partitions = cubeGeneratorParameters.partitions;
  const auto scaling = cubeGeneratorParameters.scaling;
  const auto translation = cubeGeneratorParameters.translation;

  // get cubeGenerator parameters 
  unsigned int cubeMinX = boundaries.cubeMinX;
  unsigned int cubeMaxX = boundaries.cubeMaxX;
  unsigned int cubeMinY = boundaries.cubeMinY;
  unsigned int cubeMaxY = boundaries.cubeMaxY;
  unsigned int cubeMinZ = boundaries.cubeMinZ;
  unsigned int cubeMaxZ = boundaries.cubeMaxZ;

  unsigned int cubeX = dims.cubeX;
  unsigned int cubeY = dims.cubeY;
  unsigned int cubeZ = dims.cubeZ;

  unsigned int cubePx = partitions.cubePx;
  unsigned int cubePy = partitions.cubePy;
  unsigned int cubePz = partitions.cubePz;

  double cubeScale = scaling.cubeS;
  double cubeScaleX = scaling.cubeSx;
  double cubeScaleY = scaling.cubeSy;
  double cubeScaleZ = scaling.cubeSz;

  double cubeTx = translation.cubeTx;
  double cubeTy = translation.cubeTy;
  double cubeTz = translation.cubeTz;
 
  // create additional variables necessary for cubeGenerator()
  unsigned int numCubes[4];
  unsigned int numPartitions[4];

  // check input arguments
  numCubes[0] = cubeX;
  numCubes[1] = cubeY;
  numCubes[2] = cubeZ;

  numPartitions[0] = cubePx;
  numPartitions[1] = cubePy;
  numPartitions[2] = cubePz;

  for (int i = 0; i < 3; i++) {
    if (numCubes[i] < 2)
      logError() << "Number of cubes in" << dim2str(i) << "dimension must be at least 2";
      if (numCubes[i] % numPartitions[i] != 0)
        logError() << "Number of cubes in" << dim2str(i)
                   << "dimension can not be distribute to" << numPartitions[i] << "partitions";
      if ((numCubes[i] / numPartitions[i]) % 2 != 0)
        logError() << "Number of cubes per partition in" << dim2str(i) << "dimension must be a multiple of 2";
  }

  // Compute additional sizes
  numCubes[3] = numCubes[0] * numCubes[1] * numCubes[2];
  numPartitions[3] = numPartitions[0] * numPartitions[1] * numPartitions[2];

  unsigned int numCubesPerPart[4];
  unsigned long numElemPerPart[4];
  for (int i = 0; i < 4; i++) {
    numCubesPerPart[i] = numCubes[i] / numPartitions[i];
    numElemPerPart[i] = numCubesPerPart[i] * 5;
  }

  unsigned int numVrtxPerPart[4];
  for (int i = 0; i < 3; i++)
          numVrtxPerPart[i] = numCubesPerPart[i] + 1;
  numVrtxPerPart[3] = numVrtxPerPart[0] * numVrtxPerPart[1] * numVrtxPerPart[2];

  unsigned int numBndElements[3];
  numBndElements[0] = 2*numCubesPerPart[1]*numCubesPerPart[2];
  numBndElements[1] = 2*numCubesPerPart[0]*numCubesPerPart[2];
  numBndElements[2] = 2*numCubesPerPart[0]*numCubesPerPart[1];

  // output file name
  auto cubeOutput = seissolParams.mesh.meshFileName + ".nc";  // mesh file name doesnt include file type ending
  std::string fileName = seissolParams.mesh.meshFileName + ".nc";

  // TODO: I think passing output into the cubeGenerator function breaks the call? there was a case where the output name's string representation was suddenly empty when passing output to the funciton
  // instead, replace the outut variable with its value in the function call
  logInfo() << "Start generating a mesh using the CubeGenerator";
  cubeGenerator(numCubes, numPartitions, \
                cubeMinX, cubeMaxX, cubeMinY, cubeMaxY, cubeMinZ, cubeMaxZ, \
                numCubesPerPart, numElemPerPart, numVrtxPerPart, numBndElements, \
                cubeScale, cubeScaleX, cubeScaleY, cubeScaleZ, cubeTx, cubeTy, cubeTz, fileName.c_str());
}

void seissol::initializer::initprocedure::initMesh() {
  SCOREP_USER_REGION("init_mesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  const auto commRank = seissol::MPI::mpi.rank();
  const auto commSize = seissol::MPI::mpi.size();

  logInfo(commRank) << "Begin init mesh.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<seissol::PRE_MESH>();

  const auto meshFormat = seissolParams.mesh.meshFormat;

  logInfo(commRank) << "Mesh file:" << seissolParams.mesh.meshFileName;

  seissol::Stopwatch watch;
  watch.start();

  std::string realMeshFileName = seissolParams.mesh.meshFileName;
  switch (meshFormat) {
  case seissol::geometry::MeshFormat::Netcdf:
#if USE_NETCDF
    realMeshFileName = seissolParams.mesh.meshFileName + ".nc";
    logInfo(commRank)
        << "The Netcdf file extension \".nc\" has been appended. Updated mesh file name:"
        << realMeshFileName;
    seissol::SeisSol::main.setMeshReader(
        new seissol::geometry::NetcdfReader(commRank, commSize, realMeshFileName.c_str()));
#else
    logError()
        << "Tried to load a Netcdf mesh, however this build of SeisSol is not linked to Netcdf.";
#endif
    break;
  case seissol::geometry::MeshFormat::PUML:
    readMeshPUML(seissolParams);
    break;
  case seissol::geometry::MeshFormat::CubeGenerator:
    readCubeGenerator(seissolParams);
    realMeshFileName = seissolParams.mesh.meshFileName + ".nc";
    seissol::SeisSol::main.setMeshReader(
        new seissol::geometry::NetcdfReader(commRank, commSize, realMeshFileName.c_str()));
    break;
  default:
    logError() << "Mesh reader not implemented for format" << static_cast<int>(meshFormat);
  }

  auto& meshReader = seissol::SeisSol::main.meshReader();
  postMeshread(meshReader,
               seissolParams.dynamicRupture.hasFault,
               seissolParams.mesh.displacement,
               seissolParams.mesh.scaling);

  watch.pause();
  watch.printTime("Mesh initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<seissol::POST_MESH>();

  logInfo(commRank) << "End init mesh.";

  if ((seissolParams.mesh.showEdgeCutStatistics) && (commSize > 1)) {
    logInfo(commRank) << "Computing edge cut.";
    const auto numEdges = getNumOutgoingEdges(meshReader);
    const auto summary = statistics::parallelSummary(static_cast<double>(numEdges));
    logInfo(commRank) << "Edge cut: mean =" << summary.mean << " std =" << summary.std
                      << " min =" << summary.min << " median =" << summary.median
                      << " max =" << summary.max;
  }
}
