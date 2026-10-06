// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_GEOMETRY_VERTEXORDERINVARIANCE_T_H_
#define SEISSOL_TESTS_GEOMETRY_VERTEXORDERINVARIANCE_T_H_

// mesh-permuted.h5 is mesh.h5 with the vertices of every cell listed in another order (and the
// boundary tags permuted along), as written by relabel_vertices.py --mode random of
// precomputed-seissol; for ten of its cells, that also reverses the orientation. mesh-offgrid.h5
// and mesh-offgrid-permuted.h5 are the two with their vertices moved by up to 0.01, so that sums
// of coordinates are no longer exact. Neither the partition nor anything SeisSol makes of a cell
// may depend on the order in which a cell lists its vertices.

#include "Common/CompactOptional.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/PUMLReader.h"
#include "Geometry/PartitioningLib.h"
#include "Initializer/Clustering/Clustering.h"
#include "Initializer/Clustering/VertexWeights/WeightsModels.h"
#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "TestHelper.h"

#include <PUML/DataHandle.h>
#include <PUML/Hdf5Reader.h>
#include <PUML/PUML.h>
#include <PUML/PartitionGraph.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace seissol::unit_test {

namespace vertexorderinvariancetest {

// as PUMLReader reads a mesh without periodic identification
inline void readPumlMesh(seissol::geometry::PumlMesh& mesh, const std::string& file) {
  mesh.setComm(seissol::Mpi::mpi.comm());
  {
    PUML::Hdf5Reader<seissol::geometry::PumlTopology> reader(mesh);
    reader.open(file + ":/connect", file + ":/geometry");
    reader.addData<int>(
        seissol::geometry::pumldata::Group, file + ":/group", PUML::DataType::Cell, {});
    reader.addData<std::uint32_t>(
        seissol::geometry::pumldata::Boundary, file + ":/boundary", PUML::DataType::Cell, {});
  }
  mesh.generateMesh();
}

inline std::unique_ptr<seissol::geometry::PUMLReader> readMesh(const std::string& file,
                                                               const std::string& partitioningLib,
                                                               seissol::SeisSol& seissolInstance) {
  using namespace seissol::initializer;
  static const auto FaceMap = defaultFaceMap();
  const ClusteringConfig config{
      seissol::initializer::parameters::BoundaryFormat::I32, {2}, 1, 1, 1, &FaceMap};
  Clustering clustering(config, seissolInstance);
  ExponentialWeights weightModel;
  return std::make_unique<seissol::geometry::PUMLReader>(
      tpath(file),
      partitioningLib,
      FaceMap,
      seissol::initializer::parameters::BoundaryFormat::I32,
      seissol::initializer::parameters::TopologyFormat::Geometric,
      &clustering,
      &weightModel);
}

// a mesh, and the same mesh with the vertices of its cells listed in another order
inline const std::vector<std::pair<std::string, std::string>> MeshPairs = {
    {"Testing/mesh.h5", "Testing/mesh-permuted.h5"},
    {"Testing/mesh-offgrid.h5", "Testing/mesh-offgrid-permuted.h5"}};

// every partitioner SeisSol knows; the ones it was not built with are skipped. Left out are
// ParHIPEcoMesh and ParHIPEcoSocial: their evolutionary initial partitioning runs for a fixed
// wall-clock time (2048 s over the number of ranks), so their partition differs from run to run
// anyway.
inline const std::vector<std::string> PartitioningLibs = {"Parmetis",
                                                          "ParmetisGeometric",
                                                          "PtScotch",
                                                          "PtScotchQuality",
                                                          "PtScotchBalance",
                                                          "PtScotchBalanceQuality",
                                                          "PtScotchSpeed",
                                                          "PtScotchBalanceSpeed",
                                                          "ParHIPUltrafastMesh",
                                                          "ParHIPFastMesh",
                                                          "ParHIPUltrafastSocial",
                                                          "ParHIPFastSocial"};

// a neighbor on this rank (an empty optional stands for none)
inline bool isLocal(const OptionalSize& neighbor, std::size_t elementCount) {
  return neighbor.hasValue() && neighbor.value() < elementCount;
}

} // namespace vertexorderinvariancetest

using namespace vertexorderinvariancetest;

TEST_CASE("The dual graph does not depend on the vertex order within a cell" *
          doctest::test_suite("geometry")) {
  for (const auto& meshPair : MeshPairs) {
    const auto& file = meshPair.first;
    const auto& filePermuted = meshPair.second;
    CAPTURE(file);
    seissol::geometry::PumlMesh mesh;
    seissol::geometry::PumlMesh meshPermuted;
    readPumlMesh(mesh, tpath(file));
    readPumlMesh(meshPermuted, tpath(filePermuted));

    PUML::TETPartitionGraph graph(mesh);
    const PUML::TETPartitionGraph graphPermuted(meshPermuted);

    // the partitioners see nothing but these (and the weights, which are per cell)
    CHECK(graph.vertexDistribution() == graphPermuted.vertexDistribution());
    CHECK(graph.adjDisp() == graphPermuted.adjDisp());
    CHECK(graph.adj() == graphPermuted.adj());
    for (std::size_t i = 0; i + 1 < graph.adjDisp().size(); ++i) {
      CHECK(std::is_sorted(graph.adj().begin() + graph.adjDisp()[i],
                           graph.adj().begin() + graph.adjDisp()[i + 1]));
    }

    // in double, and in the float of the ParMETIS build here (real_t)
    std::vector<double> coordinates;
    std::vector<double> coordinatesPermuted;
    graph.geometricCoordinates(coordinates);
    graphPermuted.geometricCoordinates(coordinatesPermuted);
    CHECK(coordinates == coordinatesPermuted);
    std::vector<float> coordinatesFloat;
    std::vector<float> coordinatesFloatPermuted;
    graph.geometricCoordinates(coordinatesFloat);
    graphPermuted.geometricCoordinates(coordinatesFloatPermuted);
    CHECK(coordinatesFloat == coordinatesFloatPermuted);

    // the edge ids handed out while iterating over the faces address the (sorted) adjacency
    std::vector<unsigned long> gids(mesh.cells().size());
    for (std::size_t i = 0; i < gids.size(); ++i) {
      gids[i] = mesh.cells()[i].gid();
    }
    std::vector<int> visits(graph.localEdgeCount(), 0);
    graph.forEachLocalEdges<unsigned long>(
        gids,
        [&](int /*fid*/,
            int lid,
            const unsigned long& neighbor,
            const unsigned long& self,
            int eid) {
          // no REQUIRE: the faces are iterated collectively
          const bool inRange = eid >= 0 && static_cast<std::size_t>(eid) < visits.size();
          CHECK(inRange);
          if (!inRange) {
            return;
          }
          ++visits[eid];
          CHECK(graph.adjDisp()[lid] <= static_cast<unsigned long>(eid));
          CHECK(static_cast<unsigned long>(eid) < graph.adjDisp()[lid + 1]);
          CHECK(graph.adj()[eid] == neighbor);
          CHECK(self == gids[lid]);
        });
    CHECK(std::all_of(visits.begin(), visits.end(), [](int count) { return count == 1; }));
  }
}

TEST_CASE("PUMLReader does not depend on the vertex order within a cell" *
          doctest::test_suite("geometry")) {
  std::cout.setstate(std::ios_base::failbit);

  seissol::initializer::parameters::SeisSolParameters seissolParameters{};
  seissolParameters.timeStepping.lts = seissol::initializer::parameters::LtsParameters(
      {2},
      1.0,
      0.01,
      false,
      100,
      false,
      1.0,
      seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
      seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);
  seissolParameters.timeStepping.cfl = 1;
  seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
  // varies inside the cells, so that the clustering sees the vertex order if it can
  seissolParameters.model.materialFileName = tpath("Testing/material-graded.yaml");
  seissolParameters.model.useCellHomogenizedMaterial = true;
  seissolParameters.model.plasticity = false;
  const utils::Env env("SEISSOL_");
  seissol::SeisSol seissolInstance(seissolParameters, env);

  // "None" keeps the cells where they were read. On a single rank, there is nothing to partition
  // (PUML skips the partitioner for a single part), so the others only run with several ranks.
  std::vector<std::string> partitioningLibs{"None"};
  if (seissol::Mpi::mpi.size() > 1) {
    for (const auto& partitioningLib : PartitioningLibs) {
      if (toPartitionerType(partitioningLib) != PUML::PartitionerType::None) {
        partitioningLibs.push_back(partitioningLib);
      }
    }
  }
  for (const auto& partitioningLib : partitioningLibs) {
    for (const auto& meshPair : MeshPairs) {
      const auto& file = meshPair.first;
      const auto& filePermuted = meshPair.second;
      CAPTURE(partitioningLib);
      CAPTURE(file);
      std::cout.setstate(std::ios_base::failbit);
      const auto reader = readMesh(file, partitioningLib, seissolInstance);
      const auto readerPermuted = readMesh(filePermuted, partitioningLib, seissolInstance);
      std::cout.clear();

      // with several ranks, this also requires the same partition
      const auto& elements = reader->getElements();
      const auto& elementsPermuted = readerPermuted->getElements();
      const auto& vertices = reader->getVertices();
      const auto& verticesPermuted = readerPermuted->getVertices();
      // no REQUIRE from here on: a rank leaving the test case early would leave the others waiting
      // in the collectives of the next partitioner
      CHECK(elements.size() == elementsPermuted.size());
      for (std::size_t i = 0; i < std::min(elements.size(), elementsPermuted.size()); ++i) {
        CAPTURE(i);
        const auto& element = elements[i];
        const auto& elementPermuted = elementsPermuted[i];
        CHECK(element.globalId == elementPermuted.globalId);
        CHECK(element.group == elementPermuted.group);
        CHECK(element.clusterId == elementPermuted.clusterId);
        // bit for bit
        CHECK(element.timestep == elementPermuted.timestep);
        for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
          CAPTURE(k);
          const auto& coords = vertices[element.vertices[k]].coords;
          const auto& coordsPermuted = verticesPermuted[elementPermuted.vertices[k]].coords;
          CHECK(coords == coordsPermuted);
        }
        for (std::size_t j = 0; j < Cell::NumFaces; ++j) {
          CAPTURE(j);
          // as int, so that a failure prints numbers rather than characters
          CHECK(static_cast<int>(element.neighborSides[j]) ==
                static_cast<int>(elementPermuted.neighborSides[j]));
          CHECK(static_cast<int>(element.sideOrientations[j]) ==
                static_cast<int>(elementPermuted.sideOrientations[j]));
          CHECK(element.boundaries[j] == elementPermuted.boundaries[j]);
          CHECK(element.neighborRanks[j] == elementPermuted.neighborRanks[j]);
          CHECK(element.faultTags[j] == elementPermuted.faultTags[j]);
          const bool local = isLocal(element.neighbors[j], elements.size());
          const bool localPermuted = isLocal(elementPermuted.neighbors[j], elementsPermuted.size());
          CHECK(local == localPermuted);
          if (local && localPermuted) {
            CHECK(elements[element.neighbors[j].value()].globalId ==
                  elementsPermuted[elementPermuted.neighbors[j].value()].globalId);
          }
        }
      }
    }
  }
  std::cout.clear();
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_GEOMETRY_VERTEXORDERINVARIANCE_T_H_
