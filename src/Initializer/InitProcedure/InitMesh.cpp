// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitMesh.h"

#include <Geometry/MeshDefinition.h>
#include <Initializer/Parameters/MeshParameters.h>
#include <Initializer/Parameters/SeisSolParameters.h>
#include <Initializer/TimeStepping/LtsWeights/LtsWeights.h>
#include <Solver/Estimator.h>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "utils/env.h"
#include "utils/logger.h"
#include <Eigen/Dense>
#include <math.h>
#include <mpi.h>
#include <vector>

#ifdef USE_NETCDF
#include "Geometry/CubeGenerator.h"
#include "Geometry/NetcdfReader.h"
#endif // USE_NETCDF
#if defined(USE_HDF) && defined(USE_MPI)
#include "Geometry/PUMLReader.h"
#include <hdf5.h>
#endif // defined(USE_HDF) && defined(USE_MPI)
#include "Initializer/TimeStepping/LtsWeights/WeightsFactory.h"
#include "Modules/Modules.h"
#include "Monitoring/Stopwatch.h"
#include "Numerical/Statistics.h"
#include "ResultWriter/MiniSeisSolWriter.h"
#include "SeisSol.h"

#include "Parallel/MPI.h"

namespace {

template <typename TT>
TT _checkH5Err(TT&& status, const char* file, int line) {
  if (status < 0) {
    logError() << utils::nospace << "An HDF5 error occurred (" << file << ": " << line << ")";
  }
  return std::forward<TT>(status);
}

#define _eh(status) _checkH5Err(status, __FILE__, __LINE__)

void postMeshread(seissol::geometry::MeshReader& meshReader,
                  const Eigen::Vector3d& displacement,
                  const Eigen::Matrix3d& scalingMatrix,
                  seissol::SeisSol& seissolInstance) {
  logInfo() << "The mesh has been read. Starting post processing.";

  if (meshReader.getElements().empty()) {
    logWarning(true) << "There are no local mesh elements on this rank. Is your mesh big enough?";
  }

  meshReader.displaceMesh(displacement);
  meshReader.scaleMesh(scalingMatrix);

  meshReader.computeTimestepIfNecessary(seissolInstance);

  logInfo() << "Exchanging ghostlayer metadata.";
  meshReader.exchangeGhostlayerMetadata();

  logInfo() << "Extracting fault information.";
  auto* drParameters = seissolInstance.getMemoryManager().getDRParameters();
  const VrtxCoords center{drParameters->referencePoint[0],
                          drParameters->referencePoint[1],
                          drParameters->referencePoint[2]};
  meshReader.extractFaultInformation(center, drParameters->refPointMethod);

  seissolInstance.getLtsLayout().setMesh(meshReader);

  double maxPointValue[3]{-INFINITY, -INFINITY, -INFINITY};
  double minPointValue[3]{INFINITY, INFINITY, INFINITY};

  const auto vertexCount = meshReader.getVertices().size();
#ifdef _OPENMP
#pragma omp parallel for reduction(min : minPointValue[ : 3]) reduction(max : maxPointValue[ : 3])
#endif
  for (std::size_t i = 0; i < vertexCount; ++i) {
    const auto& vertex = meshReader.getVertices()[i];
    for (int j = 0; j < 3; ++j) {
      maxPointValue[j] = std::max(maxPointValue[j], vertex.coords[j]);
      minPointValue[j] = std::min(minPointValue[j], vertex.coords[j]);
    }
  }

#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, maxPointValue, 3, MPI_DOUBLE, MPI_MAX, seissol::MPI::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE, minPointValue, 3, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm());
#endif

  logInfo() << "Smallest bounding box around the mesh: <" << minPointValue[0] << minPointValue[1]
            << minPointValue[2] << "> to <" << maxPointValue[0] << maxPointValue[1]
            << maxPointValue[2] << ">";
}

void readMeshPUML(const seissol::initializer::parameters::SeisSolParameters& seissolParams,
                  seissol::SeisSol& seissolInstance) {
#if defined(USE_HDF) && defined(USE_MPI)
  double nodeWeight = 1.0;

  if (utils::Env::get<bool>("SEISSOL_MINISEISSOL", true)) {
    if (seissol::MPI::mpi.size() > 1) {
      logInfo() << "Running mini SeisSol to determine node weights.";
      auto elapsedTime = seissol::solver::miniSeisSol();
      nodeWeight = 1.0 / elapsedTime;

      const auto summary = seissol::statistics::parallelSummary(nodeWeight);
      logInfo() << "Node weights: mean =" << summary.mean << " std =" << summary.std
                << " min =" << summary.min << " median =" << summary.median
                << " max =" << summary.max;

      writer::MiniSeisSolWriter writer(seissolParams.output.prefix.c_str());
      writer.write(elapsedTime, nodeWeight);
    } else {
      logInfo() << "Skipping mini SeisSol (SeisSol is used with a single rank only).";
    }
  } else {
    logInfo() << "Skipping mini SeisSol (disabled).";
  }

  logInfo() << "Reading PUML mesh";

  auto boundaryFormat = seissolParams.mesh.pumlBoundaryFormat;

  if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::Auto) {
    logInfo() << "Inferring boundary format.";
    MPI_Info info = MPI_INFO_NULL;
    const hid_t plistId = _eh(H5Pcreate(H5P_FILE_ACCESS));
    _eh(H5Pset_fapl_mpio(plistId, seissol::MPI::mpi.comm(), info));
    const hid_t dataFile =
        _eh(H5Fopen(seissolParams.mesh.meshFileName.c_str(), H5F_ACC_RDONLY, plistId));
    const hid_t existenceTest = _eh(H5Aexists(dataFile, "boundary-format"));
    if (existenceTest > 0) {
      logInfo() << "Boundary format given in PUML file.";
      hid_t boundaryAttribute = _eh(H5Aopen(dataFile, "boundary-format", H5P_DEFAULT));

      hid_t boundaryAttributeType = _eh(H5Aget_type(boundaryAttribute));

      auto format = [&]() {
        if (_eh(H5Tis_variable_str(boundaryAttributeType))) {
          char* formatRaw = nullptr;
          _eh(H5Aread(boundaryAttribute, boundaryAttributeType, &formatRaw));
          auto format = std::string(formatRaw);
          _eh(H5free_memory(formatRaw));
          return format;
        } else {
          auto length = H5Tget_size(boundaryAttributeType);
          std::vector<char> data(length);
          _eh(H5Aread(boundaryAttribute, boundaryAttributeType, data.data()));
          std::size_t actualLength = length;
          for (std::size_t i = 0; i < length; ++i) {
            if (data[i] == '\0') {
              actualLength = i;
              break;
            }
          }
          return std::string(data.begin(), data.begin() + actualLength);
        }
      }();

      _eh(H5Aclose(boundaryAttribute));
      _eh(H5Tclose(boundaryAttributeType));
      if (format == "i32x4") {
        boundaryFormat = seissol::initializer::parameters::BoundaryFormat::I32x4;
      } else if (format == "i64") {
        boundaryFormat = seissol::initializer::parameters::BoundaryFormat::I64;
      } else if (format == "i32") {
        boundaryFormat = seissol::initializer::parameters::BoundaryFormat::I32;
      } else {
        logError() << "Unkown boundary format given in PUML file:" << format;
      }
    } else {
      logInfo() << "Boundary format not given in PUML file; inferring from array shape.";
      const hid_t boundaryDataset = _eh(H5Dopen2(dataFile, "boundary", H5P_DEFAULT));
      const hid_t boundarySpace = _eh(H5Dget_space(boundaryDataset));
      auto boundaryTypeRank = _eh(H5Sget_simple_extent_ndims(boundarySpace));

      _eh(H5Sclose(boundarySpace));
      _eh(H5Dclose(boundaryDataset));

      // avoid checking the type size here
      if (boundaryTypeRank == 2) {
        boundaryFormat = seissol::initializer::parameters::BoundaryFormat::I32x4;
      } else if (boundaryTypeRank == 1) {
        boundaryFormat = seissol::initializer::parameters::BoundaryFormat::I32;
      } else {
        logError() << "Unknown boundary format of rank" << boundaryTypeRank;
      }
    }

    _eh(H5Fclose(dataFile));
    _eh(H5Pclose(plistId));
  }

  if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I32) {
    logInfo() << "Using boundary format: i32 (4xi8)";
  }
  if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I64) {
    logInfo() << "Using boundary format: i64 (4xi16)";
  }
  if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I32x4) {
    logInfo() << "Using boundary format: i32x4 (4xi32)";
  }

  seissol::Stopwatch watch;
  watch.start();

  using namespace seissol::initializer::time_stepping;
  const LtsWeightsConfig config{
      boundaryFormat,
      seissolParams.timeStepping.lts.getRate(),
      seissolParams.timeStepping.vertexWeight.weightElement,
      seissolParams.timeStepping.vertexWeight.weightDynamicRupture,
      seissolParams.timeStepping.vertexWeight.weightFreeSurfaceWithGravity};

  auto ltsWeights = getLtsWeightsImplementation(
      seissolParams.timeStepping.lts.getLtsWeightsType(), config, seissolInstance);
  auto* meshReader = new seissol::geometry::PUMLReader(seissolParams.mesh.meshFileName.c_str(),
                                                       seissolParams.mesh.partitioningLib.c_str(),
                                                       boundaryFormat,
                                                       ltsWeights.get(),
                                                       nodeWeight);
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

size_t getNumOutgoingEdges(seissol::geometry::MeshReader& meshReader) {
  const auto& mpiNeighbors = meshReader.getMPINeighbors();
  size_t numEdges{0};
  for (const auto& [_, neighborInfo] : mpiNeighbors) {
    // Note: this includes the case when multiple faces
    // of an element are located at a partition boarder
    numEdges += neighborInfo.elements.size();
  }
  return numEdges;
}

void readCubeGenerator(const seissol::initializer::parameters::SeisSolParameters& seissolParams,
                       seissol::SeisSol& seissolInstance) {
#if USE_NETCDF
  // unpack seissolParams
  const auto cubeParameters = seissolParams.cubeGenerator;

  const auto commRank = seissol::MPI::mpi.rank();
  const auto commSize = seissol::MPI::mpi.size();
  const std::string realMeshFileName = seissolParams.mesh.meshFileName + ".nc";

  seissolInstance.setMeshReader(
      new seissol::geometry::CubeGenerator(commRank, commSize, realMeshFileName, cubeParameters));
#else
  logError() << "Tried using CubeGenerator to read a Netcdf mesh, however this build of SeisSol is "
                "not linked to Netcdf.";
#endif
}

} // namespace

void seissol::initializer::initprocedure::initMesh(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_mesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto commRank = seissol::MPI::mpi.rank();
  const auto commSize = seissol::MPI::mpi.size();

  logInfo() << "Begin init mesh.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PreMesh>();

  const auto meshFormat = seissolParams.mesh.meshFormat;

  logInfo() << "Mesh file:" << seissolParams.mesh.meshFileName;

  seissol::Stopwatch watch;
  watch.start();

  const std::string realMeshFileName = seissolParams.mesh.meshFileName;
  [[maybe_unused]] bool addNC = true;
  if (realMeshFileName.size() >= 3) {
    const auto lastCharacters = realMeshFileName.substr(realMeshFileName.size() - 3);
    addNC = lastCharacters != ".nc";
  }

  switch (meshFormat) {
  case seissol::initializer::parameters::MeshFormat::Netcdf: {
#if USE_NETCDF
    const auto realMeshFileNameNetcdf = [&]() {
      if (addNC) {
        const auto newRealMeshFileName = realMeshFileName + ".nc";
        logInfo() << "The Netcdf file extension \".nc\" has been appended. Updated mesh file name:"
                  << newRealMeshFileName;
        return newRealMeshFileName;
      } else {
        // (suppress preference for return move)
        // NOLINTNEXTLINE
        return realMeshFileName;
      }
    }();
    seissolInstance.setMeshReader(
        new seissol::geometry::NetcdfReader(commRank, commSize, realMeshFileNameNetcdf.c_str()));
#else
    logError()
        << "Tried to load a Netcdf mesh, however this build of SeisSol is not linked to Netcdf.";
#endif
    break;
  }
  case seissol::initializer::parameters::MeshFormat::PUML: {
    readMeshPUML(seissolParams, seissolInstance);
    break;
  }
  case seissol::initializer::parameters::MeshFormat::CubeGenerator: {
    readCubeGenerator(seissolParams, seissolInstance);
    break;
  }
  default:
    logError() << "Mesh reader not implemented for format" << static_cast<int>(meshFormat);
  }

  auto& meshReader = seissolInstance.meshReader();
  postMeshread(
      meshReader, seissolParams.mesh.displacement, seissolParams.mesh.scaling, seissolInstance);

  watch.pause();
  watch.printTime("Mesh initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PostMesh>();

  logInfo() << "End init mesh.";

  if ((seissolParams.mesh.showEdgeCutStatistics) && (commSize > 1)) {
    logInfo() << "Computing edge cut.";
    const auto numEdges = getNumOutgoingEdges(meshReader);
    const auto summary = statistics::parallelSummary(static_cast<double>(numEdges));
    logInfo() << "Edge cut: mean =" << summary.mean << " std =" << summary.std
              << " min =" << summary.min << " median =" << summary.median
              << " max =" << summary.max;
  }
}
