/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
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

#ifdef USE_HDF

#include <generated_code/kernel.h>
#include "PUML/PUML.h"
#include "PUML/Downward.h"
#endif
#include <cmath>
#include <algorithm>
#include "ParameterDB.h"

#include "SeisSol.h"
#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "DynamicRupture/Misc.h"
#ifdef USE_ASAGI
#include "Reader/AsagiReader.h"
#endif
#include "utils/logger.h"

seissol::initializers::CellToVertexArray::CellToVertexArray(
    size_t size,
    const CellToVertexFunction& elementCoordinates,
    const CellToGroupFunction& elementGroups)
    : size(size), elementCoordinates(elementCoordinates), elementGroups(elementGroups) {}

seissol::initializers::CellToVertexArray seissol::initializers::CellToVertexArray::fromMeshReader(
    const seissol::geometry::MeshReader& meshReader) {
  const auto& elements = meshReader.getElements();
  const auto& vertices = meshReader.getVertices();

  return CellToVertexArray(
      elements.size(),
      [&](size_t index) {
        std::array<Eigen::Vector3d, 4> verts;
        for (size_t i = 0; i < 4; ++i) {
          auto vindex = elements[index].vertices[i];
          const auto& vertex = vertices[vindex];
          verts[i] << vertex.coords[0], vertex.coords[1], vertex.coords[2];
        }
        return verts;
      },
      [&](size_t index) { return elements[index].group; });
}

#ifdef USE_HDF
seissol::initializers::CellToVertexArray
    seissol::initializers::CellToVertexArray::fromPUML(const PUML::TETPUML& mesh) {
  const int* groups = mesh.cellData(0);
  const auto& elements = mesh.cells();
  const auto& vertices = mesh.vertices();
  return CellToVertexArray(
      elements.size(),
      [&](size_t cell) {
        std::array<Eigen::Vector3d, 4> x;
        unsigned vertLids[4];
        PUML::Downward::vertices(mesh, elements[cell], vertLids);
        for (unsigned vtx = 0; vtx < 4; ++vtx) {
          for (unsigned d = 0; d < 3; ++d) {
            x[vtx](d) = vertices[vertLids[vtx]].coordinate()[d];
          }
        }
        return x;
      },
      [groups](size_t cell) { return groups[cell]; });
}
#endif

seissol::initializers::CellToVertexArray seissol::initializers::CellToVertexArray::fromVectors(
    const std::vector<std::array<std::array<double, 3>, 4>>& vertices,
    const std::vector<int>& groups) {
  assert(vertices.size() == groups.size());

  return CellToVertexArray(
      vertices.size(),
      [&](size_t idx) {
        std::array<Eigen::Vector3d, 4> verts;
        for (size_t i = 0; i < 4; ++i) {
          verts[i] << vertices[idx][i][0], vertices[idx][i][1], vertices[idx][i][2];
        }
        return verts;
      },
      [&](size_t i) { return groups[i]; });
}

easi::Query seissol::initializers::FaultBarycentreGenerator::generate() const {
  std::vector<Fault> const& fault = m_meshReader.getFault();
  std::vector<Element> const& elements = m_meshReader.getElements();
  std::vector<Vertex> const& vertices = m_meshReader.getVertices();

  easi::Query query(m_numberOfPoints * fault.size(), 3);
  unsigned q = 0;
  for (Fault const& f : fault) {
    int element, side;
    if (f.element >= 0) {
      element = f.element;
      side = f.side;
    } else {
      element = f.neighborElement;
      side = f.neighborSide;
    }

    double barycentre[3] = {0.0, 0.0, 0.0};
    MeshTools::center(elements[element], side, vertices, barycentre);
    for (unsigned n = 0; n < m_numberOfPoints; ++n, ++q) {
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(q, dim) = barycentre[dim];
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

easi::Query seissol::initializers::FaultGPGenerator::generate() const {
  std::vector<Fault> const& fault = m_meshReader.getFault();
  std::vector<Element> const& elements = m_meshReader.getElements();
  auto cellToVertex = CellToVertexArray::fromMeshReader(m_meshReader);

  constexpr size_t numberOfPoints = dr::misc::numPaddedPoints;
  auto pointsView = init::quadpoints::view::create(const_cast<real*>(init::quadpoints::Values));
  easi::Query query(numberOfPoints * m_faceIDs.size(), 3);
  unsigned q = 0;
  // loop over all fault elements which are managed by this generator
  // note: we have one generator per LTS layer
  for (unsigned faultId : m_faceIDs) {
    const Fault& f = fault.at(faultId);
    int element, side, sideOrientation;
    if (f.element >= 0) {
      element = f.element;
      side = f.side;
      sideOrientation = -1;
    } else {
      element = f.neighborElement;
      side = f.neighborSide;
      sideOrientation = elements[f.neighborElement].sideOrientations[f.neighborSide];
    }

    auto coords = cellToVertex.elementCoordinates(element);
    for (unsigned n = 0; n < numberOfPoints; ++n, ++q) {
      double xiEtaZeta[3];
      double localPoints[2] = {pointsView(n, 0), pointsView(n, 1)};
      // padded points are in the middle of the tetrahedron
      if (n >= dr::misc::numberOfBoundaryGaussPoints) {
        localPoints[0] = 1.0 / 3.0;
        localPoints[1] = 1.0 / 3.0;
      }

      seissol::transformations::chiTau2XiEtaZeta(side, localPoints, xiEtaZeta, sideOrientation);
      Eigen::Vector3d xyz = seissol::transformations::tetrahedronReferenceToGlobal(
          coords[0], coords[1], coords[2], coords[3], xiEtaZeta);
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(q, dim) = xyz(dim);
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

namespace seissol {
namespace initializers {
using namespace seissol::model;

void FaultParameterDB::evaluateModel(std::string const& fileName,
                                     QueryGenerator const* const queryGen) {
  easi::Component* model = loadEasiModel(fileName);
  easi::Query query = queryGen->generate();

  easi::ArraysAdapter<real> adapter;
  for (auto& kv : m_parameters) {
    adapter.addBindingPoint(kv.first, kv.second.first, kv.second.second);
  }
  model->evaluate(query, adapter);

  delete model;
}

} // namespace initializers
} // namespace seissol

std::set<std::string>
    seissol::initializers::FaultParameterDB::faultProvides(std::string const& fileName) {
  if (fileName.length() == 0) {
    return std::set<std::string>();
  }
  easi::Component* model = loadEasiModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;
  return supplied;
}

seissol::initializers::EasiBoundary::EasiBoundary(const std::string& fileName)
    : model(loadEasiModel(fileName)) {}

seissol::initializers::EasiBoundary::EasiBoundary(EasiBoundary&& other)
    : model(std::move(other.model)) {}

seissol::initializers::EasiBoundary&
    seissol::initializers::EasiBoundary::operator=(EasiBoundary&& other) {
  std::swap(model, other.model);
  return *this;
}

seissol::initializers::EasiBoundary::~EasiBoundary() { delete model; }

void seissol::initializers::EasiBoundary::query(const real* nodes,
                                                real* mapTermsData,
                                                real* constantTermsData) const {
  if (model == nullptr) {
    logError() << "Model for easiBoundary is not initialized!";
  }
  assert(tensor::INodal::Shape[1] == 9); // only supp. for elastic currently.
  assert(mapTermsData != nullptr);
  assert(constantTermsData != nullptr);
  constexpr auto numNodes = tensor::INodal::Shape[0];
  auto query = easi::Query{numNodes, 3};
  size_t offset{0};
  for (unsigned i = 0; i < numNodes; ++i) {
    query.x(i, 0) = nodes[offset++];
    query.x(i, 1) = nodes[offset++];
    query.x(i, 2) = nodes[offset++];
    query.group(i) = 1;
  }
  const auto& supplied = model->suppliedParameters();

  // Shear stresses are irrelevant for riemann problem
  // Hence they have dummy names and won't be used for this bc.
  // We have 9 variables s.t. our tensors have the correct shape.
  const auto varNames =
      std::array<std::string, 9>{"Tn", "Ts", "Td", "unused1", "unused2", "unused3", "u", "v", "w"};

  // We read out a affine transformation s.t. val in ghost cell
  // is equal to A * val_inside + b
  // Note that easi only supports

  // Constant terms stores all terms of the vector b
  auto constantTerms = init::easiBoundaryConstant::view::create((constantTermsData));

  // Map terms stores all terms of the linear map A
  auto mapTerms = init::easiBoundaryMap::view::create(mapTermsData);

  easi::ArraysAdapter<real> adapter{};

  // Constant terms are named const_{varName}, e.g. const_u
  offset = 0;
  for (const auto& varName : varNames) {
    const auto termName = std::string{"const_"} + varName;
    if (supplied.count(termName) > 0) {
      adapter.addBindingPoint(termName, constantTermsData + offset, constantTerms.shape(0));
    }
    ++offset;
  }
  // Map terms are named map_{varA}_{varB}, e.g. map_u_v
  // Mirroring the velocity at the ghost cell would imply the param
  // map_u_u: -1
  offset = 0;
  for (size_t i = 0; i < varNames.size(); ++i) {
    const auto& varName = varNames[i];
    for (size_t j = 0; j < varNames.size(); ++j) {
      const auto& otherVarName = varNames[j];
      auto termName = std::string{"map_"};
      termName += varName;
      termName += "_";
      termName += otherVarName;
      if (supplied.count(termName) > 0) {
        adapter.addBindingPoint(
            termName, mapTermsData + offset, mapTerms.shape(0) * mapTerms.shape(1));
      } else {
        // Default: Extrapolate
        for (size_t k = 0; k < mapTerms.shape(2); ++k)
          mapTerms(i, j, k) = (varName == otherVarName) ? 1.0 : 0.0;
      }
      ++offset;
    }
  }
  model->evaluate(query, adapter);
}

easi::Component* seissol::initializers::loadEasiModel(const std::string& fileName) {
#ifdef USE_ASAGI
  seissol::asagi::AsagiReader asagiReader("SEISSOL_ASAGI");
  easi::YAMLParser parser(3, &asagiReader);
#else
  easi::YAMLParser parser(3);
#endif
  return parser.parse(fileName);
}
