// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_GEOMETRY_CURVEDPUML_T_H_
#define SEISSOL_TESTS_GEOMETRY_CURVEDPUML_T_H_

// The curved cells of a mesh file as PUMGen writes it: the vertices in /geometry, the other nodes
// of every cell in /geometry_ho, in the order of gmsh, with /geometry_ho_offsets and /order.
//
// cube2-*.h5 are the cube [0,1]^3 of 2^3 hexahedra cut into six tetrahedra each, written by
// makecube_ho.py with the vertices of every cell rotated by one place (o2-permuted, which also
// turns the cells inside out) or two (o3) against the order the cube is built in, so that the
// reader has to put the nodes into its own vertex order. The vertices are moved by the map below,
// and the other nodes are where it puts the nodes of the straight cells of the unmoved cube, so
// that every node of every cell is known here without the file. o2-straight has its nodes on the
// straight cells through the moved vertices instead.

#include <doctest.h>

#include "Common/Constants.h"
#include "Config.h"
#include "Geometry/IsoparametricTransform.h"
#include "Geometry/PUMLReader.h"
#include "Geometry/PartitioningLib.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "TestHelper.h"
#include "VertexOrderInvariance.t.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace seissol::unit_test {

namespace curvedpumltest {

using Point = std::array<double, Cell::Dim>;

/// the map the vertices and nodes of the test meshes are moved by
inline auto warp(const Point& u) -> Point {
  const double bump = 0.06 * std::sin(M_PI * u[0]) * std::sin(M_PI * u[1]) * std::sin(M_PI * u[2]);
  return {u[0] + bump, u[1] + 0.7 * bump, u[2] - 0.5 * bump};
}

/// its inverse, by fixed-point iteration: the map moves a point by little against its own scale
inline auto unwarp(const Point& x) -> Point {
  Point u = x;
  for (int iteration = 0; iteration < 200; ++iteration) {
    const auto moved = warp(u);
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      u[d] = x[d] - (moved[d] - u[d]);
    }
  }
  return u;
}

} // namespace curvedpumltest

TEST_CASE("PUMLReader reads the curved cells of a mesh file of PUMGen" *
          doctest::test_suite("geometry")) {
  using namespace curvedpumltest;
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
  seissolParameters.model.materialFileName = tpath("Testing/material.yaml");
  seissolParameters.model.useCellHomogenizedMaterial = true;
  seissolParameters.model.plasticity = false;
  const utils::Env env("SEISSOL_");
  seissol::SeisSol seissolInstance(seissolParameters, env);

  // a build without curved cells takes the nodes only where the cells are straight-sided, and
  // forgets them; a curved mesh it refuses, which ends the run
  const std::vector<std::pair<std::string, std::size_t>> files =
      Config::Curvilinear
          ? std::vector<std::pair<std::string, std::size_t>>{{"Testing/cube2-o2-permuted.h5", 2},
                                                             {"Testing/cube2-o3.h5", 3},
                                                             {"Testing/cube2-o2-straight.h5", 2}}
          : std::vector<std::pair<std::string, std::size_t>>{{"Testing/cube2-o2-straight.h5", 1}};

  // "None" keeps the cells where they were read; a partitioner moves them, and their nodes with
  // them
  std::vector<std::string> partitioningLibs{"None"};
  if (seissol::Mpi::mpi.size() > 1 &&
      toPartitionerType("Parmetis") != PUML::PartitionerType::None) {
    partitioningLibs.emplace_back("Parmetis");
  }

  for (const auto& partitioningLib : partitioningLibs) {
    for (const auto& [file, order] : files) {
      CAPTURE(partitioningLib);
      CAPTURE(file);
      const auto reader =
          vertexorderinvariancetest::readMesh(file, partitioningLib, seissolInstance);
      // no REQUIRE: a rank leaving the test case early would leave the others waiting
      CHECK(reader->geometryOrder() == order);
      if (reader->geometryOrder() != order || order <= 1) {
        continue;
      }

      const auto lattice = seissol::geometry::IsoparametricTransform::latticeNodes(order);
      const auto& elements = reader->getElements();
      const auto& vertices = reader->getVertices();
      const bool straight = file.find("straight") != std::string::npos;
      double offset = 0;
      for (std::size_t cell = 0; cell < elements.size(); ++cell) {
        const auto nodes = reader->cellNodes(cell);
        CHECK(nodes.size() == lattice.size());
        std::array<Point, Cell::NumVertices> corners{};
        for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
          const auto& coords = vertices[elements[cell].vertices[k]].coords;
          corners[k] = straight ? Point{coords[0], coords[1], coords[2]}
                                : unwarp({coords[0], coords[1], coords[2]});
        }
        for (std::size_t i = 0; i < std::min(nodes.size(), lattice.size()); ++i) {
          // the node at the barycentric coordinates of the lattice node, on the local vertices
          const std::array<double, Cell::NumVertices> weights{1 - lattice[i](0) - lattice[i](1) -
                                                                  lattice[i](2),
                                                              lattice[i](0),
                                                              lattice[i](1),
                                                              lattice[i](2)};
          Point point{};
          for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
            for (std::size_t d = 0; d < Cell::Dim; ++d) {
              point[d] += weights[k] * corners[k][d];
            }
          }
          const auto expected = straight ? point : warp(point);
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            offset = std::max(offset, std::abs(nodes[i][d] - expected[d]));
          }
        }
      }
      CHECK(offset == AbsApprox(0.0).epsilon(1e-12));
    }
  }
  std::cout.clear();
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_GEOMETRY_CURVEDPUML_T_H_
