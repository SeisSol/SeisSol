// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_GEOMETRY_HIGHORDERPUML_T_H_
#define SEISSOL_TESTS_GEOMETRY_HIGHORDERPUML_T_H_

// The cells of a higher order of a mesh file as PUMGen writes it: the vertices in /geometry, the
// other nodes of every cell in /geometry_ho, in the order of gmsh, with /geometry_ho_offsets and
// /order. cube2-*.h5 are the cube [0,1]^3 of 2^3 hexahedra cut into six tetrahedra each, with its
// vertices moved by a smooth map (written by makecube_ho.py of davschneller/nodal). In
// o2-straight, the other nodes lie on the straight-sided cells through the moved vertices; in
// o2-permuted and o3, they are where the map puts the nodes of the straight cells of the unmoved
// cube, so that the cells are curved.

#include <doctest.h>

#include "Geometry/PUMLReader.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "TestHelper.h"
#include "VertexOrderInvariance.t.h"

#include <PUML/Hdf5Reader.h>
#include <PUML/PUML.h>
#include <cstddef>
#include <iostream>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

namespace seissol::unit_test {

namespace highorderpumltest {

/// what the nodes of a mesh file say about its cells, read as PUMLReader reads them before it
/// checks them
inline auto summaryOf(const std::string& file) -> seissol::geometry::HighOrderGeometrySummary {
  seissol::geometry::PumlMesh mesh;
  mesh.setComm(seissol::Mpi::mpi.comm());
  {
    PUML::Hdf5Reader<seissol::geometry::PumlTopology> reader(mesh);
    reader.open(tpath(file) + ":/connect", tpath(file) + ":/geometry");
  }
  const auto arrays = seissol::geometry::PUMLReader::readHighOrderGeometry(mesh, tpath(file));
  mesh.generateMesh();
  return seissol::geometry::summarizeHighOrderGeometry(mesh, arrays);
}

} // namespace highorderpumltest

TEST_CASE("PUMLReader tells the curved cells of a mesh file of PUMGen from the straight ones" *
          doctest::test_suite("geometry")) {
  using seissol::geometry::PUMLReader;
  std::cout.setstate(std::ios_base::failbit);

  const auto straight = highorderpumltest::summaryOf("Testing/cube2-o2-straight.h5");
  CHECK(straight.highOrderCells == 48);
  CHECK(straight.curvedCells == 0);
  CHECK(straight.order == 2);
  CHECK(straight.curvature < PUMLReader::StraightTolerance);

  const std::vector<std::pair<std::string, unsigned long>> curvedFiles = {
      {"Testing/cube2-o2-permuted.h5", 2}, {"Testing/cube2-o3.h5", 3}};
  for (const auto& curvedFile : curvedFiles) {
    const auto& file = curvedFile.first;
    CAPTURE(file);
    const auto curved = highorderpumltest::summaryOf(file);
    CHECK(curved.highOrderCells == 48);
    CHECK(curved.curvedCells > 0);
    CHECK(curved.order == curvedFile.second);
    CHECK(curved.curvature > PUMLReader::StraightTolerance);
  }

  std::cout.clear();
}

TEST_CASE("PUMLReader reads straight-sided cells of a higher order as cells of the first order" *
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
  seissolParameters.model.materialFileName = tpath("Testing/material.yaml");
  seissolParameters.model.useCellHomogenizedMaterial = true;
  seissolParameters.model.plasticity = false;
  const utils::Env env("SEISSOL_");
  seissol::SeisSol seissolInstance(seissolParameters, env);

  // the reader takes the mesh, and forgets the nodes beyond the vertices
  const auto reader =
      vertexorderinvariancetest::readMesh("Testing/cube2-o2-straight.h5", "None", seissolInstance);
  unsigned long cells = reader->getElements().size();
  MPI_Allreduce(MPI_IN_PLACE, &cells, 1, MPI_UNSIGNED_LONG, MPI_SUM, seissol::Mpi::mpi.comm());
  CHECK(cells == 48);

  std::cout.clear();
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_GEOMETRY_HIGHORDERPUML_T_H_
