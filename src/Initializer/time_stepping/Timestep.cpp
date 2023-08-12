#include "Timestep.hpp"

#include "Initializer/ParameterMaterialDB.hpp"
#include <vector>
#include <array>
#include <functional>
#include <Eigen/Dense>

#include "Equations/datastructures.hpp"
#include "Initializer/ParameterDB.h"
#include "Initializer/InputParameters.hpp"
#include "Initializer/ConfigFile.hpp"

#include "SeisSol.h"

#include <variant>

namespace seissol::initializer {
static double computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                                  double pWaveVel,
                                  double cfl,
                                  double maximumAllowedTimeStep) {
  // Compute insphere radius
  std::array<Eigen::Vector3d, 4> x = vertices;
  Eigen::Matrix4d A;
  A << x[0](0), x[0](1), x[0](2), 1.0, x[1](0), x[1](1), x[1](2), 1.0, x[2](0), x[2](1), x[2](2),
      1.0, x[3](0), x[3](1), x[3](2), 1.0;

  double alpha = A.determinant();
  double Nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
  double Nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
  double Nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
  double Nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
  double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);

  // Compute maximum timestep
  return std::fmin(maximumAllowedTimeStep,
                   cfl * 2.0 * insphere / (pWaveVel * (2 * ConvergenceOrder - 1)));
}

static int determineCellCluster(double timestep, double startTimestep, double rate) {
  int cluster = -1;
  double clusterTimestep = startTimestep;
  while (timestep >= clusterTimestep) {
    clusterTimestep *= rate;
    ++cluster;
  }
  return cluster;
}

static int enforceMaximumDifferenceLocal(int maxDiff,
                                         std::vector<int>& localClustering,
                                         const geometry::MeshReader& meshReader) {
  int numberOfReductions = 0;

  const auto& cells = meshReader.getElements();

  int reductionsPerStep = 0;
  do {
    reductionsPerStep = 0;
    for (unsigned cell = 0; cell < cells.size(); ++cell) {
      int timeCluster = localClustering[cell];

      for (unsigned f = 0; f < 4; ++f) {
        int difference = maxDiff;
        int boundary = cells[cell].boundaries[f];
        // Continue for regular, dynamic rupture, and periodic boundary cells
        if (boundary == 0 || boundary == 3 || boundary == 6) {
          // We treat MPI neighbours later
          if (cells[cell].neighbors[f] < cells.size()) {
            int otherTimeCluster = localClustering[cells[cell].neighbors[f]];

            if (boundary == 3) {
              difference = 0;
            }

            if (timeCluster > otherTimeCluster + difference) {
              timeCluster = otherTimeCluster + difference;
              ++reductionsPerStep;
            }
          }
        }
      }
      localClustering[cell] = timeCluster;
    }
    numberOfReductions += reductionsPerStep;
  } while (reductionsPerStep > 0);

#ifdef USE_MPI
  auto numExchanges = meshReader.getMPINeighbors().size();
  std::vector<MPI_Request> requests(2 * numExchanges);
  std::vector<std::vector<int>> ghost(numExchanges);
  std::vector<std::vector<int>> copy(numExchanges);

  auto exchange = meshReader.getMPINeighbors().begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.elements.size();
    ghost[ex].resize(exchangeSize);
    copy[ex].resize(exchangeSize);

    for (unsigned n = 0; n < exchangeSize; ++n) {
      copy[ex][n] = localClustering[cells[exchange->second.elements[n].localElement]
                                        .mpiIndices[exchange->second.elements[n].localSide]];
    }
    MPI_Isend(copy[ex].data(),
              exchangeSize,
              MPI_INT,
              exchange->first,
              0,
              seissol::MPI::mpi.comm(),
              &requests[ex]);
    MPI_Irecv(ghost[ex].data(),
              exchangeSize,
              MPI_INT,
              exchange->first,
              0,
              seissol::MPI::mpi.comm(),
              &requests[numExchanges + ex]);
    ++exchange;
  }

  MPI_Waitall(2 * numExchanges, requests.data(), MPI_STATUSES_IGNORE);

  exchange = meshReader.getMPINeighbors().begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.elements.size();
    for (unsigned n = 0; n < exchangeSize; ++n) {
      int difference = maxDiff;
      int otherTimeCluster = ghost[ex][n];

      auto cell = exchange->second.elements[n].localElement;
      auto f = exchange->second.elements[n].localSide;

      int boundary = cells[cell].boundaries[f];
      if (boundary == 3) {
        difference = 0;
      }

      if (localClustering[cell] > otherTimeCluster + difference) {
        localClustering[cell] = otherTimeCluster + difference;
        ++numberOfReductions;
      }
    }
    ++exchange;
  }

#endif // USE_MPI

  return numberOfReductions;
}

GlobalTimestep computeTimesteps(double cfl,
                                double maximumAllowedTimeStep,
                                const seissol::initializers::CellToVertexArray& cellToVertex) {
  // TODO(David): remove parsing here
  auto configs = seissol::initializer::readConfigFile(
      seissol::SeisSol::main.getSeisSolParameters().model.configFileName);
  auto [materials, plasticity] = seissol::model::queryMaterial(configs, cellToVertex, false);

  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned cell = 0; cell < cellToVertex.size; ++cell) {
    double pWaveVel = 0;
    std::visit([&](auto&& material) { pWaveVel = material.getMaxWaveSpeed(); }, materials[cell]);
    std::array<Eigen::Vector3d, 4> vertices = cellToVertex.elementCoordinates(cell);
    timestep.cellTimeStepWidths[cell] =
        computeCellTimestep(vertices, pWaveVel, cfl, maximumAllowedTimeStep);
  }

  const auto minmaxCellPosition =
      std::minmax_element(timestep.cellTimeStepWidths.begin(), timestep.cellTimeStepWidths.end());

  double localMinTimestep = *minmaxCellPosition.first;
  double localMaxTimestep = *minmaxCellPosition.second;

#ifdef USE_MPI
  MPI_Allreduce(&localMinTimestep,
                &timestep.globalMinTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                seissol::MPI::mpi.comm());
  MPI_Allreduce(&localMaxTimestep,
                &timestep.globalMaxTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                seissol::MPI::mpi.comm());
#else
  timestep.globalMinTimeStep = localMinTimestep;
  timestep.globalMaxTimeStep = localMaxTimestep;
#endif
  return timestep;
}

std::vector<int> clusterTimesteps(const GlobalTimestep& timestep, double rate, double wiggle) {
  std::vector<int> cluster(timestep.cellTimeStepWidths.size());
  const double startTimestep = wiggle * timestep.globalMinTimeStep;
  for (unsigned cell = 0; cell < cluster.size(); ++cell) {
    cluster[cell] = determineCellCluster(timestep.cellTimeStepWidths[cell], startTimestep, rate);
  }
  return cluster;
}

int enforceMaximumDifference(int maxDiff,
                             std::vector<int>& localClustering,
                             const geometry::MeshReader& meshReader) {
  int totalReductions = 0;
  int stepReductions = 0;
  do {
    stepReductions = enforceMaximumDifferenceLocal(maxDiff, localClustering, meshReader);
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, stepReductions, 1, MPI_INT, MPI_SUM, seissol::MPI::mpi.comm());
#endif
  } while (stepReductions > 0);
  return totalReductions;
}

} // namespace seissol::initializer
