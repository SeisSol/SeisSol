/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
#include <PUML/PUML.h>
#include <PUML/Downward.h>
#endif
#include "ParameterDB.h"

#include <easi/YAMLParser.h>
#include <easi/ResultAdapter.h>
#include <Numerical_aux/Transformation.h>
#ifdef USE_ASAGI
#include <Reader/AsagiReader.h>
#endif
#include <utils/logger.h>

easi::Query seissol::initializers::ElementBarycentreGenerator::generate() const {
  std::vector<Element> const& elements = m_meshReader.getElements();
  std::vector<Vertex> const& vertices = m_meshReader.getVertices();
  
  easi::Query query(elements.size(), 3);
  for (unsigned elem = 0; elem < elements.size(); ++elem) {
    // Compute barycentre for each element
    for (unsigned dim = 0; dim < 3; ++dim) {
      query.x(elem,dim) = vertices[ elements[elem].vertices[0] ].coords[dim];
    }
    for (unsigned vertex = 1; vertex < 4; ++vertex) {
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(elem,dim) += vertices[ elements[elem].vertices[vertex] ].coords[dim];
      }
    }
    for (unsigned dim = 0; dim < 3; ++dim) {
      query.x(elem,dim) *= 0.25;
    }
    // Group
    query.group(elem) = elements[elem].material;
  }
  return query;
}

#ifdef USE_HDF
easi::Query seissol::initializers::ElementBarycentreGeneratorPUML::generate() const {
  std::vector<PUML::TETPUML::cell_t> const& cells = m_mesh.cells();
	std::vector<PUML::TETPUML::vertex_t> const& vertices = m_mesh.vertices();
  int const* material = m_mesh.cellData(0);
  
  easi::Query query(cells.size(), 3);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    unsigned vertLids[4];
    PUML::Downward::vertices(m_mesh, cells[cell], vertLids);
    
    // Compute barycentre for each element
    for (unsigned dim = 0; dim < 3; ++dim) {
      query.x(cell,dim) = vertices[ vertLids[0] ].coordinate()[dim];
    }
    for (unsigned vertex = 1; vertex < 4; ++vertex) {
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(cell,dim) += vertices[ vertLids[vertex] ].coordinate()[dim];
      }
    }
    for (unsigned dim = 0; dim < 3; ++dim) {
      query.x(cell,dim) *= 0.25;
    }
    // Group
    query.group(cell) = material[cell];
  }
  return query;
}
#endif

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
    // Compute barycentre
    for (unsigned vertex = 0; vertex < 3; ++vertex) {
      for (unsigned dim = 0; dim < 3; ++dim) {
        barycentre[dim] += vertices[ elements[element].vertices[ MeshTools::FACE2NODES[side][vertex] ] ].coords[dim];
      }
    }
    for (unsigned dim = 0; dim < 3; ++dim) {
      barycentre[dim] /= 3.0;
    }
    for (unsigned n = 0; n < m_numberOfPoints; ++n, ++q) {
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(q,dim) = barycentre[dim];
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

easi::Query seissol::initializers::FaultGPGenerator::generate() const {
  std::vector<Fault> const& fault = m_meshReader.getFault();
  std::vector<Element> const& elements = m_meshReader.getElements();
  std::vector<Vertex> const& vertices = m_meshReader.getVertices();

  easi::Query query(m_numberOfPoints * fault.size(), 3);
  unsigned q = 0;
  for (Fault const& f : fault) {
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

    double const* coords[4];
    for (unsigned v = 0; v < 4; ++v) {
      coords[v] = vertices[ elements[element].vertices[ v ] ].coords;
    }
    for (unsigned n = 0; n < m_numberOfPoints; ++n, ++q) {
      double xiEtaZeta[3], xyz[3];
      seissol::transformations::chiTau2XiEtaZeta(side, m_points[n], xiEtaZeta, sideOrientation);
      seissol::transformations::tetrahedronReferenceToGlobal(coords[0], coords[1], coords[2], coords[3], xiEtaZeta, xyz);
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(q,dim) = xyz[dim];
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

easi::Component* seissol::initializers::ParameterDB::loadModel(std::string const& fileName) {
#ifdef USE_ASAGI
  seissol::asagi::AsagiReader asagiReader("SEISSOL_ASAGI");
  easi::YAMLParser parser(3, &asagiReader);
#else
  easi::YAMLParser parser(3);
#endif
  easi::Component* model = parser.parse(fileName);
  return model;
}

void seissol::initializers::ParameterDB::evaluateModel(std::string const& fileName, QueryGenerator const& queryGen) {
  easi::ArraysAdapter adapter;
  for (auto& kv : m_parameters) {
    adapter.addBindingPoint(kv.first, kv.second.first, kv.second.second);
  }
  
  easi::Query query = queryGen.generate();
  easi::Component* model = ParameterDB::loadModel(fileName);
  model->evaluate(query, adapter);  
  delete model;
}

bool seissol::initializers::ParameterDB::faultParameterizedByTraction(std::string const& fileName) {
  easi::Component* model = ParameterDB::loadModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;

  std::set<std::string> stress = {"s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"};
  std::set<std::string> traction =  {"T_n", "T_s", "T_d"};

  bool containsStress = std::includes(supplied.begin(), supplied.end(), stress.begin(), stress.end());
  bool containsTraction = std::includes(supplied.begin(), supplied.end(), traction.begin(), traction.end());

  if (containsStress == containsTraction) {
    logError() << "Both stress (s_xx, s_yy, s_zz, s_xy, s_yz, s_xz) and traction (T_n, T_s, T_d) are defined (or are missing), but only either of them must be defined.";
  }

  return containsTraction;
}
