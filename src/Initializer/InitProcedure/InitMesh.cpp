#include "Init.hpp"
#include "InitMesh.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>

#include "Eigen/Dense"
#include "utils/logger.h"
#include "utils/env.h"

#ifdef USE_NETCDF
#include "Geometry/NetcdfReader.h"
#include "Geometry/CubeGenerator.h"
#endif // USE_NETCDF
#if defined(USE_HDF) && defined(USE_MPI)
#include "Geometry/PUMLReader.h"
#endif // defined(USE_HDF) && defined(USE_MPI)
#include "Modules/Modules.h"
#include "Monitoring/instrumentation.hpp"
#include "Monitoring/Stopwatch.h"
#include "Numerical_aux/Statistics.h"
#include "Initializer/time_stepping/LtsWeights/WeightsFactory.h"
#include "Solver/time_stepping/MiniSeisSol.h"
#include "SeisSol.h"
#include "ResultWriter/MiniSeisSolWriter.h"

#include "Parallel/MPI.h"

namespace {

static void postMeshread(seissol::geometry::MeshReader& meshReader,
                         const Eigen::Vector3d& displacement,
                         const Eigen::Matrix3d& scalingMatrix,
                         seissol::SeisSol& seissolInstance) {
  logInfo(seissol::MPI::mpi.rank()) << "The mesh has been read. Starting post processing.";

  if (meshReader.getElements().empty()) {
    logWarning(seissol::MPI::mpi.rank())
        << "There are no local mesh elements on this rank. Is your mesh big enough?";
  }

  meshReader.displaceMesh(displacement);
  meshReader.scaleMesh(scalingMatrix);

  logInfo(seissol::MPI::mpi.rank()) << "Extracting fault information.";

  auto* drParameters = seissolInstance.getMemoryManager().getDRParameters();
  VrtxCoords center{drParameters->referencePoint[0],
                    drParameters->referencePoint[1],
                    drParameters->referencePoint[2]};
  meshReader.extractFaultInformation(center, drParameters->refPointMethod);

  logInfo(seissol::MPI::mpi.rank()) << "Exchanging ghostlayer metadata.";
  meshReader.exchangeGhostlayerMetadata();

  seissolInstance.getLtsLayout().setMesh(meshReader);
}

static void readMeshPUML(const seissol::initializer::parameters::SeisSolParameters& seissolParams,
                         seissol::SeisSol& seissolInstance) {
#if defined(USE_HDF) && defined(USE_MPI)
  const int rank = seissol::MPI::mpi.rank();
  double nodeWeight = 1.0;

  if (utils::Env::get<bool>("SEISSOL_MINISEISSOL", true)) {
    if (seissol::MPI::mpi.size() > 1) {
      logInfo(rank) << "Running mini SeisSol to determine node weights.";
      auto elapsedTime = seissol::miniSeisSol(
          seissolInstance.getMemoryManager(), seissolParams.model.plasticity, seissolInstance);
      nodeWeight = 1.0 / elapsedTime;

      const auto summary = seissol::statistics::parallelSummary(nodeWeight);
      logInfo(rank) << "Node weights: mean =" << summary.mean << " std =" << summary.std
                    << " min =" << summary.min << " median =" << summary.median
                    << " max =" << summary.max;

      writer::MiniSeisSolWriter writer(seissolParams.output.prefix.c_str());
      writer.write(elapsedTime, nodeWeight);
    } else {
      logInfo(rank) << "Skipping mini SeisSol (SeisSol is used with a single rank only).";
    }
  } else {
    logInfo(rank) << "Skipping mini SeisSol (disabled).";
  }

  logInfo(rank) << "Reading PUML mesh";

  seissol::Stopwatch watch;
  watch.start();

  bool readPartitionFromFile = seissolInstance.simulator().checkPointingEnabled();

  using namespace seissol::initializer::time_stepping;
  LtsWeightsConfig config{seissolParams.model.materialFileName,
                          static_cast<unsigned int>(seissolParams.timeStepping.lts.getRate()),
                          seissolParams.timeStepping.vertexWeight.weightElement,
                          seissolParams.timeStepping.vertexWeight.weightDynamicRupture,
                          seissolParams.timeStepping.vertexWeight.weightFreeSurfaceWithGravity};

  auto ltsWeights = getLtsWeightsImplementation(
      seissolParams.timeStepping.lts.getLtsWeightsType(), config, seissolInstance);
  auto meshReader =
      new seissol::geometry::PUMLReader(seissolParams.mesh.meshFileName.c_str(),
                                        seissolParams.mesh.partitioningLib.c_str(),
                                        seissolParams.timeStepping.maxTimestepWidth,
                                        seissolParams.output.checkpointParameters.fileName.c_str(),
                                        ltsWeights.get(),
                                        nodeWeight,
                                        readPartitionFromFile);
  seissolInstance.setMeshReader(meshReader);

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

static size_t getNumOutgoingEdges(seissol::geometry::MeshReader& meshReader) {
  auto& mpiNeighbors = meshReader.getMPINeighbors();
  size_t numEdges{0};
  for (auto& [_, neighborInfo] : mpiNeighbors) {
    // Note: this includes the case when multiple faces
    // of an element are located at a partition boarder
    numEdges += neighborInfo.elements.size();
  }
  return numEdges;
}

} // namespace

static void
    readCubeGenerator(const seissol::initializer::parameters::SeisSolParameters& seissolParams,
                      seissol::SeisSol& seissolInstance) {
#if USE_NETCDF
  // unpack seissolParams
  const auto cubeParameters = seissolParams.cubeGenerator;

  const auto commRank = seissol::MPI::mpi.rank();
  const auto commSize = seissol::MPI::mpi.size();
  std::string realMeshFileName = seissolParams.mesh.meshFileName + ".nc";
  auto meshReader = new seissol::geometry::CubeGenerator(
      commRank, commSize, realMeshFileName.c_str(), cubeParameters);

  // Replace call to NetcdfReader with adapted Geometry/CubeGenerator
  seissolInstance.setMeshReader(
      new seissol::geometry::NetcdfReader(commRank, commSize, realMeshFileName.c_str()));
#else
  logError() << "Tried using CubeGenerator to read a Netcdf mesh, however this build of SeisSol is "
                "not linked to Netcdf.";
#endif
}

void seissol::initializer::initprocedure::initMesh(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_mesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const auto& seissolParams = seissolInstance.getSeisSolParameters();
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
  case seissol::initializer::parameters::MeshFormat::Netcdf:
#if USE_NETCDF
    realMeshFileName = seissolParams.mesh.meshFileName + ".nc";
    logInfo(commRank)
        << "The Netcdf file extension \".nc\" has been appended. Updated mesh file name:"
        << realMeshFileName;
    seissolInstance.setMeshReader(
        new seissol::geometry::NetcdfReader(commRank, commSize, realMeshFileName.c_str()));
#else
    logError()
        << "Tried to load a Netcdf mesh, however this build of SeisSol is not linked to Netcdf.";
#endif
    break;
  case seissol::initializer::parameters::MeshFormat::PUML:
    readMeshPUML(seissolParams, seissolInstance);
    break;
  case seissol::initializer::parameters::MeshFormat::CubeGenerator:
    readCubeGenerator(seissolParams, seissolInstance);
    break;
  default:
    logError() << "Mesh reader not implemented for format" << static_cast<int>(meshFormat);
  }

  auto& meshReader = seissolInstance.meshReader();
  postMeshread(
      meshReader, seissolParams.mesh.displacement, seissolParams.mesh.scaling, seissolInstance);

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
