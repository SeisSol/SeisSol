/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
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
#include "ParameterDB.h"

#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "Numerical_aux/Transformation.h"
#ifdef USE_ASAGI
#include "Reader/AsagiReader.h"
#endif
#include "utils/logger.h"


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
    MeshTools::center(elements[element], side, vertices, barycentre);
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

namespace seissol {
  namespace initializers {
    template<>
    void MaterialParameterDB<seissol::model::ElasticMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<seissol::model::ElasticMaterial> &adapter) {
      adapter.addBindingPoint("rho", &seissol::model::ElasticMaterial::rho);
      adapter.addBindingPoint("mu", &seissol::model::ElasticMaterial::mu);
      adapter.addBindingPoint("lambda", &seissol::model::ElasticMaterial::lambda);
    }

    template<>
    void MaterialParameterDB<seissol::model::ViscoElasticMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<seissol::model::ViscoElasticMaterial> &adapter) {
      adapter.addBindingPoint("rho", &seissol::model::ViscoElasticMaterial::rho);
      adapter.addBindingPoint("mu", &seissol::model::ViscoElasticMaterial::mu);
      adapter.addBindingPoint("lambda", &seissol::model::ViscoElasticMaterial::lambda);
      adapter.addBindingPoint("Qp", &seissol::model::ViscoElasticMaterial::Qp);
      adapter.addBindingPoint("Qs", &seissol::model::ViscoElasticMaterial::Qs);
    }

    template<>
    void MaterialParameterDB<seissol::model::PoroElasticMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<seissol::model::PoroElasticMaterial> &adapter) {
      adapter.addBindingPoint("bulk_solid", &seissol::model::PoroElasticMaterial::bulkSolid);
      adapter.addBindingPoint("rho", &seissol::model::PoroElasticMaterial::rho);
      adapter.addBindingPoint("lambda", &seissol::model::PoroElasticMaterial::lambda);
      adapter.addBindingPoint("mu", &seissol::model::PoroElasticMaterial::mu);
      adapter.addBindingPoint("porosity", &seissol::model::PoroElasticMaterial::porosity);
      adapter.addBindingPoint("permeability", &seissol::model::PoroElasticMaterial::permeability);
      adapter.addBindingPoint("tortuosity", &seissol::model::PoroElasticMaterial::tortuosity);
      adapter.addBindingPoint("bulk_fluid", &seissol::model::PoroElasticMaterial::bulkFluid);
      adapter.addBindingPoint("rho_fluid", &seissol::model::PoroElasticMaterial::rhoFluid);
      adapter.addBindingPoint("viscosity", &seissol::model::PoroElasticMaterial::viscosity);
    }

    template<>
    void MaterialParameterDB<seissol::model::Plasticity>::addBindingPoints(easi::ArrayOfStructsAdapter<seissol::model::Plasticity> &adapter) {
      adapter.addBindingPoint("bulkFriction", &seissol::model::Plasticity::bulkFriction);
      adapter.addBindingPoint("plastCo", &seissol::model::Plasticity::plastCo);
      adapter.addBindingPoint("s_xx", &seissol::model::Plasticity::s_xx);
      adapter.addBindingPoint("s_yy", &seissol::model::Plasticity::s_yy);
      adapter.addBindingPoint("s_zz", &seissol::model::Plasticity::s_zz);
      adapter.addBindingPoint("s_xy", &seissol::model::Plasticity::s_xy);
      adapter.addBindingPoint("s_yz", &seissol::model::Plasticity::s_yz);
      adapter.addBindingPoint("s_xz", &seissol::model::Plasticity::s_xz);
    }

    template<>
    void MaterialParameterDB<seissol::model::AnisotropicMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<seissol::model::AnisotropicMaterial> &adapter) {
      adapter.addBindingPoint("rho", &seissol::model::AnisotropicMaterial::rho);
      adapter.addBindingPoint("c11", &seissol::model::AnisotropicMaterial::c11);
      adapter.addBindingPoint("c12", &seissol::model::AnisotropicMaterial::c12);
      adapter.addBindingPoint("c13", &seissol::model::AnisotropicMaterial::c13);
      adapter.addBindingPoint("c14", &seissol::model::AnisotropicMaterial::c14);
      adapter.addBindingPoint("c15", &seissol::model::AnisotropicMaterial::c15);
      adapter.addBindingPoint("c16", &seissol::model::AnisotropicMaterial::c16);
      adapter.addBindingPoint("c22", &seissol::model::AnisotropicMaterial::c22);
      adapter.addBindingPoint("c23", &seissol::model::AnisotropicMaterial::c23);
      adapter.addBindingPoint("c24", &seissol::model::AnisotropicMaterial::c24);
      adapter.addBindingPoint("c25", &seissol::model::AnisotropicMaterial::c25);
      adapter.addBindingPoint("c26", &seissol::model::AnisotropicMaterial::c26);
      adapter.addBindingPoint("c33", &seissol::model::AnisotropicMaterial::c33);
      adapter.addBindingPoint("c34", &seissol::model::AnisotropicMaterial::c34);
      adapter.addBindingPoint("c35", &seissol::model::AnisotropicMaterial::c35);
      adapter.addBindingPoint("c36", &seissol::model::AnisotropicMaterial::c36);
      adapter.addBindingPoint("c44", &seissol::model::AnisotropicMaterial::c44);
      adapter.addBindingPoint("c45", &seissol::model::AnisotropicMaterial::c45);
      adapter.addBindingPoint("c46", &seissol::model::AnisotropicMaterial::c46);
      adapter.addBindingPoint("c55", &seissol::model::AnisotropicMaterial::c55);
      adapter.addBindingPoint("c56", &seissol::model::AnisotropicMaterial::c56);
      adapter.addBindingPoint("c66", &seissol::model::AnisotropicMaterial::c66);
    }                                                               
    
    template<class T>
    void MaterialParameterDB<T>::evaluateModel(std::string const& fileName, QueryGenerator const& queryGen) {
      easi::Component* model = loadEasiModel(fileName);
      easi::Query query = queryGen.generate();

      easi::ArrayOfStructsAdapter<T> adapter(m_materials->data());
      addBindingPoints(adapter);
      model->evaluate(query, adapter); 
      
      delete model;
    }
    
    template<>
    void MaterialParameterDB<seissol::model::AnisotropicMaterial>::evaluateModel(std::string const& fileName, QueryGenerator const& queryGen) {
      easi::Component* model = loadEasiModel(fileName);
      easi::Query query = queryGen.generate();
      auto suppliedParameters = model->suppliedParameters();
      //TODO(Sebastian): inhomogeneous materials, where in some parts only mu and lambda are given
      //                 and in other parts the full elastic tensor is given

      //if we look for an anisotropic material and only mu and lambda are supplied, 
      //assume isotropic behavior and calculate the parameters accordingly
      if (suppliedParameters.find("mu") != suppliedParameters.end() && suppliedParameters.find("lambda") != suppliedParameters.end()) {
        std::vector<seissol::model::ElasticMaterial> elasticMaterials(query.numPoints());
        easi::ArrayOfStructsAdapter<seissol::model::ElasticMaterial> adapter(elasticMaterials.data());
        MaterialParameterDB<seissol::model::ElasticMaterial>().addBindingPoints(adapter);
        unsigned numPoints = query.numPoints();
        model->evaluate(query, adapter);

        for(unsigned i = 0; i < numPoints; i++) {
          m_materials->at(i) = seissol::model::AnisotropicMaterial(elasticMaterials[i]);
        }
      }
      else {
        easi::ArrayOfStructsAdapter<seissol::model::AnisotropicMaterial> arrayOfStructsAdapter(m_materials->data());
        addBindingPoints(arrayOfStructsAdapter);
        model->evaluate(query, arrayOfStructsAdapter);  
        
      }
      delete model;
    }

    void FaultParameterDB::evaluateModel(std::string const& fileName, QueryGenerator const& queryGen) {
      easi::Component* model = loadEasiModel(fileName);
      easi::Query query = queryGen.generate();

      easi::ArraysAdapter<double> adapter;
      for (auto& kv : m_parameters) {
        adapter.addBindingPoint(kv.first, kv.second.first, kv.second.second);
      }
      model->evaluate(query, adapter); 
      
      delete model;
    }

  }
}

bool seissol::initializers::FaultParameterDB::faultParameterizedByTraction(std::string const& fileName) {
  easi::Component* model = loadEasiModel(fileName);
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

bool seissol::initializers::FaultParameterDB::nucleationParameterizedByTraction(std::string const& fileName) {
  easi::Component* model = loadEasiModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;

  std::set<std::string> stress = {"nuc_xx", "nuc_yy", "nuc_zz", "nuc_xy", "nuc_yz", "nuc_xz"};
  std::set<std::string> traction =  {"Tnuc_n", "Tnuc_s", "Tnuc_d"};

  bool containsStress = std::includes(supplied.begin(), supplied.end(), stress.begin(), stress.end());
  bool containsTraction = std::includes(supplied.begin(), supplied.end(), traction.begin(), traction.end());

  if (containsStress == containsTraction) {
    logError() << "Both nucleation stress (nuc_xx, nuc_yy, nuc_zz, nuc_xy, nuc_yz, nuc_xz) and nucleation traction (Tnuc_n, Tnuc_s, Tnuc_d) are defined (or are missing), but only either of them must be defined.";
  }

  return containsTraction;
}



seissol::initializers::EasiBoundary::EasiBoundary(const std::string& fileName)
  : model(loadEasiModel(fileName)) {
}


seissol::initializers::EasiBoundary::EasiBoundary(EasiBoundary&& other)
  : model(std::move(other.model)) {}

seissol::initializers::EasiBoundary& seissol::initializers::EasiBoundary::operator=(EasiBoundary&& other) {
  std::swap(model, other.model);
  return *this;
}

seissol::initializers::EasiBoundary::~EasiBoundary() {
  delete model;
}

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
  const auto varNames = std::array<std::string, 9>{
    "Tn", "Ts", "Td", "unused1", "unused2", "unused3", "u", "v", "w"
  };

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
      adapter.addBindingPoint(
          termName,
          constantTermsData + offset,
          constantTerms.shape(0)
      );
    }
    ++offset;
  }
  // Map terms are named map_{varA}_{varB}, e.g. map_u_v
  // Mirroring the velocity at the ghost cell would imply the param
  // map_u_u: -1
  offset = 0;
  for (size_t i = 0; i < varNames.size(); ++i){
    const auto& varName = varNames[i];
    for (size_t j = 0; j < varNames.size(); ++j)  {
      const auto& otherVarName = varNames[j];
      auto termName = std::string{"map_"};
      termName += varName;
      termName += "_";
      termName += otherVarName;
      if (supplied.count(termName) > 0) {
        adapter.addBindingPoint(
            termName,
            mapTermsData + offset,
            mapTerms.shape(0) * mapTerms.shape(1)
        );
      } else {
        // Default: Extrapolate
        for (size_t k = 0; k < mapTerms.shape(2); ++k)
          mapTerms(i,j,k) = (varName == otherVarName) ? 1.0 : 0.0;
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

template class seissol::initializers::MaterialParameterDB<seissol::model::AnisotropicMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::ElasticMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::ViscoElasticMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::PoroElasticMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::Plasticity>;


