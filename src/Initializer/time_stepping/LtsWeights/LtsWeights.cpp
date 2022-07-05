/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017 - 2020, SeisSol Group
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
 *
 **/
#include <Eigen/Eigenvalues>
#include <Kernels/precision.hpp>
#include <Initializer/typedefs.hpp>

#include <PUML/PUML.h>
#include <PUML/Downward.h>
#include <PUML/Upward.h>
#include "LtsWeights.h"
#include "WeightsModels.h"

#include <Initializer/ParameterDB.h>
#include <Parallel/MPI.h>

#include <generated_code/init.h>

namespace seissol::initializers::time_stepping {

class FaceSorter {
  private:
  std::vector<PUML::TETPUML::face_t> const& m_faces;

  public:
  FaceSorter(std::vector<PUML::TETPUML::face_t> const& faces) : m_faces(faces) {}

  bool operator()(unsigned int a, unsigned int b) const {
    return m_faces[a].gid() < m_faces[b].gid();
  }
};

void LtsWeights::computeNodeWeights(PUML::TETPUML const& mesh, double maximumAllowedTimeStep) {
  logInfo(seissol::MPI::mpi.rank()) << "Computing LTS weights.";

  // Note: Return value optimization is guaranteed while returning temp. objects in C++17
  m_mesh = &mesh;
  m_details = collectGlobalTimeStepDetails(maximumAllowedTimeStep);
  m_clusterIds = computeClusterIds();
  m_ncon = evaluateNumberOfConstraints();
  auto totalNumberOfReductions = enforceMaximumDifference();
  m_cellCosts = computeCostsPerTimestep();

  if (!m_vertexWeights.empty()) {
    m_vertexWeights.clear();
  }
  m_vertexWeights.resize(m_clusterIds.size() * m_ncon);

  // calling virtual functions
  setVertexWeights();
  setAllowedImbalances();

  logInfo(seissol::MPI::mpi.rank()) << "Computing LTS Node Weights. Done. " << utils::nospace << '('
                                    << totalNumberOfReductions << " reductions.)";
}

void LtsWeights::computeEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                               const std::vector<idx_t>&>& graph) {
  if (!m_edgeWeights.empty()) {
    m_edgeWeights.clear();
  }

  size_t edge_count = std::get<2>(graph).size();

  m_edgeWeights.resize(edge_count);

  setEdgeWeights(graph);

  logInfo(seissol::MPI::mpi.rank()) << "Computing LTS Edge Weights. Done. " << utils::nospace;
}

const int* LtsWeights::vertexWeights() const {
  assert(!m_vertexWeights.empty() && "vertex weights are not initialized");
  return m_vertexWeights.data();
}

const int* LtsWeights::edgeWeights() const {
  assert(!m_edgeWeights.empty() && "edge weights are not initialized");
  return m_edgeWeights.data();
}

int LtsWeights::edgeCount() const {
  assert(!m_edgeWeights.empty() && "edge weights are not initialized");
  return m_edgeWeights.size();
}

const double* LtsWeights::imbalances() const {
  assert(!m_imbalances.empty() && "weight imbalances are not initialized");
  return m_imbalances.data();
}

int LtsWeights::nWeightsPerVertex() const {
  assert(m_ncon != std::numeric_limits<int>::infinity() &&
         "num. constrains has not been initialized yet");
  return m_ncon;
}

void LtsWeights::computeMaxTimesteps(std::vector<double> const& pWaveVel,
                                     std::vector<double>& timeSteps,
                                     double maximumAllowedTimeStep) {
  std::vector<PUML::TETPUML::cell_t> const& cells = m_mesh->cells();
  std::vector<PUML::TETPUML::vertex_t> const& vertices = m_mesh->vertices();

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    // Compute insphere radius
    Eigen::Vector3d barycentre(0., 0., 0.);
    Eigen::Vector3d x[4];
    unsigned vertLids[4];
    PUML::Downward::vertices(*m_mesh, cells[cell], vertLids);
    for (unsigned vtx = 0; vtx < 4; ++vtx) {
      for (unsigned d = 0; d < 3; ++d) {
        x[vtx](d) = vertices[vertLids[vtx]].coordinate()[d];
      }
    }
    Eigen::Matrix4d A;
    A << x[0](0), x[0](1), x[0](2), 1.0, x[1](0), x[1](1), x[1](2), 1.0, x[2](0), x[2](1), x[2](2),
        1.0, x[3](0), x[3](1), x[3](2), 1.0;

    double alpha = A.determinant();
    double Nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
    double Nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
    double Nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
    double Nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
    double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);

    // Compute maximum timestep (CFL=1)
    timeSteps[cell] = std::fmin(maximumAllowedTimeStep,
                                2.0 * insphere / (pWaveVel[cell] * (2 * CONVERGENCE_ORDER - 1)));
  }
}

int LtsWeights::getCluster(double timestep, double globalMinTimestep, unsigned rate) {
  if (rate == 1) {
    return 0;
  }

  double upper = rate * globalMinTimestep;

  int cluster = 0;
  while (upper <= timestep) {
    upper *= rate;
    ++cluster;
  }
  return cluster;
}

int LtsWeights::getBoundaryCondition(int const* boundaryCond, unsigned cell, unsigned face) {
  int bcCurrentFace = ((boundaryCond[cell] >> (face * 8)) & 0xFF);
  if (bcCurrentFace > 64) {
    bcCurrentFace = 3;
  }
  return bcCurrentFace;
}

int LtsWeights::ipow(int x, int y) {
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while (--y) {
    result *= x;
  }
  return result;
}

LtsWeights::GlobalTimeStepDetails
LtsWeights::collectGlobalTimeStepDetails(double maximumAllowedTimeStep) {

  const auto& cells = m_mesh->cells();
  std::vector<double> pWaveVel;
  pWaveVel.resize(cells.size());

  seissol::initializers::ElementBarycentreGeneratorPUML queryGen(*m_mesh);
  // up to now we only distinguish between anisotropic elastic any other isotropic material
#ifdef USE_ANISOTROPIC
  std::vector<seissol::model::AnisotropicMaterial> materials(cells.size());
  seissol::initializers::MaterialParameterDB<seissol::model::AnisotropicMaterial> parameterDB;
#elif defined(USE_POROELASTIC)
  std::vector<seissol::model::PoroElasticMaterial> materials(cells.size());
  seissol::initializers::MaterialParameterDB<seissol::model::PoroElasticMaterial> parameterDB;
#else
  std::vector<seissol::model::ElasticMaterial> materials(cells.size());
  seissol::initializers::MaterialParameterDB<seissol::model::ElasticMaterial> parameterDB;
#endif
  parameterDB.setMaterialVector(&materials);
  parameterDB.evaluateModel(m_velocityModel, queryGen);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    pWaveVel[cell] = materials[cell].getMaxWaveSpeed();
  }

  GlobalTimeStepDetails details{};
  details.timeSteps.resize(cells.size());
  computeMaxTimesteps(pWaveVel, details.timeSteps, maximumAllowedTimeStep);

  double localMinTimestep = *std::min_element(details.timeSteps.begin(), details.timeSteps.end());
  double localMaxTimestep = *std::max_element(details.timeSteps.begin(), details.timeSteps.end());

#ifdef USE_MPI
  MPI_Allreduce(&localMinTimestep, &details.globalMinTimeStep, 1, MPI_DOUBLE, MPI_MIN,
                seissol::MPI::mpi.comm());
  MPI_Allreduce(&localMaxTimestep, &details.globalMaxTimeStep, 1, MPI_DOUBLE, MPI_MAX,
                seissol::MPI::mpi.comm());
#else
  details.globalMinTimeStep = localMinTimestep;
  details.globalMaxTimeStep = localMaxTimestep;
#endif
  return details;
}

std::vector<int> LtsWeights::computeClusterIds() {
  const auto& cells = m_mesh->cells();
  std::vector<int> clusterIds(cells.size(), 0);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    clusterIds[cell] = getCluster(m_details.timeSteps[cell], m_details.globalMinTimeStep, m_rate);
  }
  return clusterIds;
}

std::vector<int> LtsWeights::computeCostsPerTimestep() {
  const auto& cells = m_mesh->cells();

  std::vector<int> cellCosts(cells.size());
  int const* boundaryCond = m_mesh->cellData(1);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int dynamicRupture = 0;
    int freeSurfaceWithGravity = 0;

    unsigned int faceids[4];
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);

    for (unsigned face = 0; face < 4; ++face) {
      const auto faceType = static_cast<FaceType>(getBoundaryCondition(boundaryCond, cell, face));
      dynamicRupture += (faceType == FaceType::dynamicRupture) ? 1 : 0;
      freeSurfaceWithGravity += (faceType == FaceType::freeSurfaceGravity) ? 1 : 0;
    }

    const int costDynamicRupture = m_vertexWeightDynamicRupture * dynamicRupture;
    const int costDisplacement = m_vertexWeightFreeSurfaceWithGravity * freeSurfaceWithGravity;
    cellCosts[cell] = m_vertexWeightElement + costDynamicRupture + costDisplacement;
  }
  return cellCosts;
}

int LtsWeights::enforceMaximumDifference() {
  int totalNumberOfReductions = 0;
  int globalNumberOfReductions;
  do {
    int localNumberOfReductions = enforceMaximumDifferenceLocal();

#ifdef USE_MPI
    MPI_Allreduce(&localNumberOfReductions, &globalNumberOfReductions, 1, MPI_INT, MPI_SUM,
                  seissol::MPI::mpi.comm());
#else
    globalNumberOfReductions = localNumberOfReductions;
#endif
    totalNumberOfReductions += globalNumberOfReductions;
  } while (globalNumberOfReductions > 0);
  return totalNumberOfReductions;
}

int LtsWeights::enforceMaximumDifferenceLocal(int maxDifference) {
  int numberOfReductions = 0;

  std::vector<PUML::TETPUML::cell_t> const& cells = m_mesh->cells();
  std::vector<PUML::TETPUML::face_t> const& faces = m_mesh->faces();
  int const* boundaryCond = m_mesh->cellData(1);

#ifdef USE_MPI
  std::unordered_map<int, std::vector<int>> rankToSharedFaces;
  std::unordered_map<int, int> localFaceIdToLocalCellId;
#endif

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int timeCluster = m_clusterIds[cell];

    unsigned int faceids[4];
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);
    for (unsigned f = 0; f < 4; ++f) {
      int difference = maxDifference;
      int boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (boundary == 0 || boundary == 3 || boundary == 6) {
        // We treat MPI neighbours later
        auto const& face = faces[faceids[f]];
        if (!face.isShared()) {
          int cellIds[2];
          PUML::Upward::cells(*m_mesh, face, cellIds);

          int neighbourCell = (cellIds[0] == static_cast<int>(cell)) ? cellIds[1] : cellIds[0];
          int otherTimeCluster = m_clusterIds[neighbourCell];

          if (boundary == 3) {
            difference = 0;
          }

          if (timeCluster > otherTimeCluster + difference) {
            timeCluster = otherTimeCluster + difference;
            ++numberOfReductions;
          }
        }
#ifdef USE_MPI
        else {
          rankToSharedFaces[face.shared()[0]].push_back(faceids[f]);
          localFaceIdToLocalCellId[faceids[f]] = cell;
        }
#endif
      }
    }
    m_clusterIds[cell] = timeCluster;
  }

#ifdef USE_MPI
  FaceSorter faceSorter(faces);
  for (auto& sharedFaces : rankToSharedFaces) {
    std::sort(sharedFaces.second.begin(), sharedFaces.second.end(), faceSorter);
  }

  auto numExchanges = rankToSharedFaces.size();
  std::vector<MPI_Request> requests(2 * numExchanges);
  std::vector<std::vector<int>> ghost(numExchanges);
  std::vector<std::vector<int>> copy(numExchanges);

  auto exchange = rankToSharedFaces.begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.size();
    ghost[ex].resize(exchangeSize);
    copy[ex].resize(exchangeSize);

    for (unsigned n = 0; n < exchangeSize; ++n) {
      copy[ex][n] = m_clusterIds[localFaceIdToLocalCellId[exchange->second[n]]];
    }
    MPI_Isend(copy[ex].data(), exchangeSize, MPI_INT, exchange->first, 0, seissol::MPI::mpi.comm(),
              &requests[ex]);
    MPI_Irecv(ghost[ex].data(), exchangeSize, MPI_INT, exchange->first, 0, seissol::MPI::mpi.comm(),
              &requests[numExchanges + ex]);
    ++exchange;
  }

  MPI_Waitall(2 * numExchanges, requests.data(), MPI_STATUSES_IGNORE);

  exchange = rankToSharedFaces.begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.size();
    for (unsigned n = 0; n < exchangeSize; ++n) {
      int difference = maxDifference;
      int otherTimeCluster = ghost[ex][n];

      int cellIds[2];
      PUML::Upward::cells(*m_mesh, faces[exchange->second[n]], cellIds);
      int cell = (cellIds[0] >= 0) ? cellIds[0] : cellIds[1];

      unsigned int faceids[4];
      PUML::Downward::faces(*m_mesh, cells[cell], faceids);
      unsigned f = 0;
      for (; f < 4 && static_cast<int>(faceids[f]) != exchange->second[n]; ++f)
        ;
      assert(f != 4);

      int boundary = getBoundaryCondition(boundaryCond, cell, f);
      if (boundary == 3) {
        difference = 0;
      }

      if (m_clusterIds[cell] > otherTimeCluster + difference) {
        m_clusterIds[cell] = otherTimeCluster + difference;
        ++numberOfReductions;
      }
    }
    ++exchange;
  }

#endif

  return numberOfReductions;
}

const LtsWeights::GlobalTimeStepDetails& LtsWeights::getDetails() const { return m_details; }

const std::string& LtsWeights::getVelocityModel() const { return m_velocityModel; }

unsigned LtsWeights::getRate() const { return m_rate; }

std::vector<int>& LtsWeights::getVertexWeights() { return m_vertexWeights; }

std::vector<double>& LtsWeights::getImbalances() { return m_imbalances; }

const std::vector<int>& LtsWeights::getCellCosts() const { return m_cellCosts; }

int LtsWeights::getVertexWeightElement() const { return m_vertexWeightElement; }

int LtsWeights::getVertexWeightDynamicRupture() const { return m_vertexWeightDynamicRupture; }

int LtsWeights::getVertexWeightFreeSurfaceWithGravity() const {
  return m_vertexWeightFreeSurfaceWithGravity;
}

int LtsWeights::getNcon() const { return m_ncon; }

const PUML::TETPUML* LtsWeights::getMesh() const { return m_mesh; }

const std::vector<int>& LtsWeights::getClusterIds() const { return m_clusterIds; }

LtsWeights::LtsWeights(const LtsWeightsConfig& config) noexcept
    : m_velocityModel(config.velocityModel), m_rate(config.rate),
      m_vertexWeightElement(config.vertexWeightElement),
      m_vertexWeightDynamicRupture(config.vertexWeightDynamicRupture),
      m_vertexWeightFreeSurfaceWithGravity(config.vertexWeightFreeSurfaceWithGravity),
      m_nodeWeightModel(nullptr), m_edgeWeightModel(nullptr) {}

void LtsWeights::setVertexWeights() { m_nodeWeightModel->setVertexWeights(); }

void LtsWeights::setAllowedImbalances() { m_nodeWeightModel->setAllowedImbalances(); }

int LtsWeights::evaluateNumberOfConstraints() const {
  return m_nodeWeightModel->evaluateNumberOfConstraints();
}

void LtsWeights::setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                           const std::vector<idx_t>&>& graph) {
  m_edgeWeightModel->setEdgeWeights(graph);
}

LtsWeights::~LtsWeights() noexcept {
  delete m_nodeWeightModel;
  delete m_edgeWeightModel;
}

void LtsWeights::addWeightModels(NodeWeightModel* nwm, EdgeWeightModel* ewm) {
  m_nodeWeightModel = nwm;
  m_edgeWeightModel = ewm;
}

std::vector<int>& LtsWeights::getEdgeWeights() { return m_edgeWeights; }

int LtsWeights::find_rank(const std::vector<idx_t>& vrtxdist, idx_t elemId){
  assert(!vrtxdist.empty() && "Vrtxdist should not be empty, implementation in PUML could have a bug");
  assert(elemId < vrtxdist.back() && "Local ids should be within vrtxdist, implementation in PUML could have a bug");
  assert(elemId >= 0 && "Element id should not be empty, implementation in PUML could have a bug");
  assert(vrtxdist[0] == 0 && "Metis returns the vrtxdist where the first element is always 0, implementation in PUML could have a bug (unauthorized write");

  for (size_t i = 0; i < vrtxdist.size()-1; i++){
    if (elemId >= vrtxdist[i] && elemId < vrtxdist[i+1])
    {
      return i;
    }
  }

  throw std::runtime_error("element id is not in valid range!");
}


void LtsWeights::exchangeGhostLayer(const std::tuple<const std::vector<idx_t>&,
                                                     const std::vector<idx_t>&,
                                                     const std::vector<idx_t>&>& graph) {

  if (!ghost_layer_info.empty())
  {
    return;
  }
  
  int mpillint_size = 0;
  MPI_Type_size(MPI_LONG_LONG_INT,&mpillint_size);
  assert(sizeof(idx_t) <= static_cast<unsigned int>(mpillint_size) && "For ghost layer exchange to work the size of idx_t has to be least or equal to MPI_LONG_LONG_INT");
  
  const std::vector<idx_t>& vrtxdist = std::get<0>(graph);
  const std::vector<idx_t>& xadj = std::get<1>(graph);
  const std::vector<idx_t>& adjncy = std::get<2>(graph);

  const int rank = seissol::MPI::mpi.rank();
  const int size = seissol::MPI::mpi.size();

  assert((vrtxdist.size() == static_cast<size_t>(size + 1) && !xadj.empty() && !adjncy.empty()) && "The dual graph should be initialized for the ghost layer exchange");
  assert(vrtxdist[0] == 0 && "Metis always returns 0 for the first element of vrtxdist, either the graph is not intiliazed or there was an unllowed write to vrtxdist[0]");

  const size_t vertex_id_begin = vrtxdist[rank];

  assert(!m_clusterIds.empty() && "Cluster Ids need to updated by this point");

  // gather all the ranks that have neighbors to these ranks
  std::set<int> neighbor_ranks;

  // ghost layers of other ranks
  // ghost_layer_to_send[i] means the actors and their cluster to send to rank i!
  // using a vector because we need continous data and therefore we will have
  // a vector of idx_t the length of element count * 2
  std::vector<std::vector<idx_t>> ghost_layer_to_send;
  std::vector<std::vector<idx_t>> ghost_layer_received;
  std::vector<std::unordered_map<idx_t, int>> ghost_layer_mapped;
  std::vector<MPI_Request> send_requests;
  std::vector<MPI_Status> recv_stats;

  // the at(my_rank) will be empty but it makes coding much easier
  ghost_layer_to_send.resize(size);
  ghost_layer_received.resize(size);
  ghost_layer_mapped.resize(size);

  for (size_t i = 0; i < xadj.size() - 1; i++) {
    //[xadj[i] .. xadj[i+1]) describes to offsets neighbor ids of local vertex i
    const size_t neighbors_offset_begin = xadj[i];
    const size_t neighbors_offset_end = xadj[i + 1];

    const int cluster_id = m_clusterIds[i];

    for (size_t j = neighbors_offset_begin; j < neighbors_offset_end; j++) {
      const idx_t neighbor_id = adjncy[j];

      const int neighbor_rank = find_rank(vrtxdist, neighbor_id);

      // if neighbor's rank is different than mine
      if (neighbor_rank != rank) {
        neighbor_ranks.insert(neighbor_rank);
        ghost_layer_to_send[neighbor_rank].push_back(
            static_cast<idx_t>(i + vertex_id_begin)); // send global id
        ghost_layer_to_send[neighbor_rank].push_back(static_cast<idx_t>(cluster_id));
      }
    }
  }

  // if we add empty requests and stats waiting becomes problematic
  send_requests.resize(neighbor_ranks.size());
  recv_stats.resize(neighbor_ranks.size());

  // now send ghost_layer_to_send to neighbors
  // per default sets are stored in ascending order
  {
    int offset = 0;
    for (int neighbor_rank : neighbor_ranks) {
      const std::vector<idx_t>& layer_ref = ghost_layer_to_send[neighbor_rank];

      if (sizeof(idx_t) != 8)
      {
        if (rank == 0)
        {
          std::cerr << "If the size of idx_t is not 64 bits the ghost layer exchange will crush with a segmentation fault!" << std::endl;
        }
      } 

      MPI_Isend(layer_ref.data(), layer_ref.size(), MPI_LONG_LONG_INT, neighbor_rank, rank,
                MPI_COMM_WORLD, &send_requests[offset]);
      

      offset += 1;
    }
  }

  volatile unsigned int got = 0;
  std::vector<bool> skip;
  skip.resize(neighbor_ranks.size());
  std::fill(skip.begin(), skip.end(), false);

  // now lets probe and try to get the ghost layer
  while (got < neighbor_ranks.size()) {
    int offset = 0;
    for (int neighbor_rank : neighbor_ranks) {
      if (neighbor_rank == rank) {
        throw std::runtime_error("A neighbor's rank should not be the same as self rank!");
      }

      if (skip[offset]) {
        offset += 1;
      } else {
        int flag = 0;
        MPI_Iprobe(neighbor_rank, neighbor_rank, MPI_COMM_WORLD, &flag, &recv_stats[offset]);

        if (flag) {
          int count = 0;

          MPI_Get_count(&recv_stats[offset], MPI_LONG_LONG_INT, &count);

          assert(ghost_layer_received[neighbor_rank].empty() && "After ghost layer exchange the ghost layer map should not be empty");
          ghost_layer_received[neighbor_rank].resize(count);
          assert(ghost_layer_mapped[neighbor_rank].empty() && "After ghost layer exchange the ghost layer map should not be empty");
          MPI_Recv(ghost_layer_received[neighbor_rank].data(), count, MPI_LONG_LONG_INT, neighbor_rank,
                   neighbor_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          got += 1;

          skip[offset] = true;
        }

        offset += 1;
      }
    }
  }

  // wait for the communication to be completed
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);

  // not sure if necessary as waiting for all send implies the loop will
  // completed too and the rest of communication is totally local

  for (int i : neighbor_ranks) {
    assert(i != rank && "A neighbor should not be our own local rank");
    assert(!ghost_layer_received[i].empty() && "Ghost layer received should not be empty");
    assert(ghost_layer_received[i].size() % 2 == 0 && "Ghost layer received should have an even number elements n cell ids and n cluster ids therefore always 2n");

    for (unsigned int j = 0; j < ghost_layer_received[i].size(); j += 2) {
      ghost_layer_mapped[i].emplace(ghost_layer_received[i][j], ghost_layer_received[i][j + 1]);
    }
  }

  ghost_layer_info = std::move(ghost_layer_mapped);
}

void LtsWeights::apply_constraints(
                        const std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>& graph,
                        std::vector<idx_t> &constraint_to_update, 
                        std::function<int(idx_t, idx_t)>& factor,
                        OffsetType ot)
{
   exchangeGhostLayer(graph);

  const int rank = seissol::MPI::mpi.rank();
  
  const std::vector<idx_t>& vrtxdist = std::get<0>(graph);
  const std::vector<idx_t>& xadj = std::get<1>(graph);
  const std::vector<idx_t>& adjncy = std::get<2>(graph);
  
  const size_t vertex_id_begin = vrtxdist[rank];

  // compute edge weights with the ghost layer
  for (size_t i = 0; i < xadj.size() - 1; i++) {
    // get cluster of the cell i
    const int self_cluster_id = m_clusterIds[i];

    const size_t neighbors_offset_begin = xadj[i];
    const size_t neighbors_offset_end = xadj[i + 1];

    for (size_t j = neighbors_offset_begin; j < neighbors_offset_end; j++) {
      const idx_t neighbor_id_idx = adjncy[j];

      const int rank_of_neighbor = find_rank(vrtxdist, neighbor_id_idx);

      // if  on the same rank get from cluster_ids otherwise get it from the received ghost_layer

      const int other_cluster_id =
          (rank_of_neighbor == rank)
              ? m_clusterIds[neighbor_id_idx -
                             vertex_id_begin] // mapping from global id to local id
              : ghost_layer_info[rank_of_neighbor].find(neighbor_id_idx)->second;
      
      if (ot == OffsetType::edgeWeight)
      {
        constraint_to_update[j] = factor(self_cluster_id, other_cluster_id);
      } else if (ot == OffsetType::minMsg) {
        constexpr int constraint_beg_offset = 2;
        constraint_to_update[(m_ncon * i) + constraint_beg_offset] += factor(self_cluster_id, other_cluster_id);
      } else {
        constexpr int constraint_beg_offset = 2;
        constraint_to_update[(m_ncon * i) + constraint_beg_offset + other_cluster_id] += factor(self_cluster_id, other_cluster_id);
      }
    }
  }
}

} // namespace seissol::initializers::time_stepping