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
    query.group(elem) = elements[elem].group;
  }
  return query;
}

seissol::initializers::ElementAverageGenerator::ElementAverageGenerator(MeshReader const& meshReader)
  : m_meshReader(meshReader)
  {
    double quadraturePoints[NUM_QUADPOINTS][3];
    double quadratureWeights[NUM_QUADPOINTS];
    seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, CONVERGENCE_ORDER);

    std::copy(std::begin(quadratureWeights), std::end(quadratureWeights), std::begin(m_quadratureWeights));
    for (int i = 0; i < NUM_QUADPOINTS; ++i) {
      std::copy(std::begin(quadraturePoints[i]), std::end(quadraturePoints[i]), std::begin(m_quadraturePoints[i]));
    }
  }

easi::Query seissol::initializers::ElementAverageGenerator::generate() const {
  std::vector<Element> const& elements = m_meshReader.getElements();
  std::vector<Vertex> const& vertices = m_meshReader.getVertices();

  // Generate query using quadrature points for each element
  easi::Query query(elements.size() * NUM_QUADPOINTS, 3);
  
  // Transform quadrature points to global coordinates for all elements
  #pragma omp parallel for
  for (unsigned elem = 0; elem < elements.size(); ++elem) {
    for (unsigned i = 0; i < NUM_QUADPOINTS; ++i) {
      std::array<double, 3> xyz{};
      seissol::transformations::tetrahedronReferenceToGlobal(vertices[ elements[elem].vertices[0] ].coords, vertices[ elements[elem].vertices[1] ].coords,
        vertices[ elements[elem].vertices[2] ].coords, vertices[ elements[elem].vertices[3] ].coords, m_quadraturePoints[i].data(), xyz.data());
      for (unsigned dim = 0; dim < 3; ++dim) {
        query.x(elem * NUM_QUADPOINTS + i,dim) = xyz[dim];
      }
      query.group(elem * NUM_QUADPOINTS + i) = elements[elem].group;
    }
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

  constexpr size_t numberOfPoints = dr::misc::numPaddedPoints;
  auto pointsView = init::quadpoints::view::create(const_cast<real *>(init::quadpoints::Values));
  easi::Query query(numberOfPoints * m_faceIDs.size(), 3);
  unsigned q = 0;
  // loop over all fault elements which are managed by this generator
  // note: we have one generator per LTS layer
  for (unsigned faultId: m_faceIDs) {
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

    double const* coords[4];
    for (unsigned v = 0; v < 4; ++v) {
      coords[v] = vertices[ elements[element].vertices[ v ] ].coords;
    }
    for (unsigned n = 0; n < numberOfPoints; ++n, ++q) {
      double xiEtaZeta[3], xyz[3];
      double localPoints[2] = {pointsView(n,0), pointsView(n,1)};
      // padded points are in the middle of the tetrahedron
      if (n >= dr::misc::numberOfBoundaryGaussPoints) {
        localPoints[0] = 1.0/3.0;
        localPoints[1] = 1.0/3.0;
      }

      seissol::transformations::chiTau2XiEtaZeta(side, localPoints, xiEtaZeta, sideOrientation);
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
    using namespace seissol::model;

    template<>
    void MaterialParameterDB<ElasticMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<ElasticMaterial> &adapter) {
      adapter.addBindingPoint("rho", &ElasticMaterial::rho);
      adapter.addBindingPoint("mu", &ElasticMaterial::mu);
      adapter.addBindingPoint("lambda", &ElasticMaterial::lambda);
    }

    template<>
    void MaterialParameterDB<ViscoElasticMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<ViscoElasticMaterial> &adapter) {
      adapter.addBindingPoint("rho", &ViscoElasticMaterial::rho);
      adapter.addBindingPoint("mu", &ViscoElasticMaterial::mu);
      adapter.addBindingPoint("lambda", &ViscoElasticMaterial::lambda);
      adapter.addBindingPoint("Qp", &ViscoElasticMaterial::Qp);
      adapter.addBindingPoint("Qs", &ViscoElasticMaterial::Qs);
    }

    template<>
    void MaterialParameterDB<PoroElasticMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<PoroElasticMaterial> &adapter) {
      adapter.addBindingPoint("bulk_solid", &PoroElasticMaterial::bulkSolid);
      adapter.addBindingPoint("rho", &PoroElasticMaterial::rho);
      adapter.addBindingPoint("lambda", &PoroElasticMaterial::lambda);
      adapter.addBindingPoint("mu", &PoroElasticMaterial::mu);
      adapter.addBindingPoint("porosity", &PoroElasticMaterial::porosity);
      adapter.addBindingPoint("permeability", &PoroElasticMaterial::permeability);
      adapter.addBindingPoint("tortuosity", &PoroElasticMaterial::tortuosity);
      adapter.addBindingPoint("bulk_fluid", &PoroElasticMaterial::bulkFluid);
      adapter.addBindingPoint("rho_fluid", &PoroElasticMaterial::rhoFluid);
      adapter.addBindingPoint("viscosity", &PoroElasticMaterial::viscosity);
    }

    template<>
    void MaterialParameterDB<Plasticity>::addBindingPoints(easi::ArrayOfStructsAdapter<Plasticity> &adapter) {
      adapter.addBindingPoint("bulkFriction", &Plasticity::bulkFriction);
      adapter.addBindingPoint("plastCo", &Plasticity::plastCo);
      adapter.addBindingPoint("s_xx", &Plasticity::s_xx);
      adapter.addBindingPoint("s_yy", &Plasticity::s_yy);
      adapter.addBindingPoint("s_zz", &Plasticity::s_zz);
      adapter.addBindingPoint("s_xy", &Plasticity::s_xy);
      adapter.addBindingPoint("s_yz", &Plasticity::s_yz);
      adapter.addBindingPoint("s_xz", &Plasticity::s_xz);
    }

    template<>
    void MaterialParameterDB<AnisotropicMaterial>::addBindingPoints(easi::ArrayOfStructsAdapter<AnisotropicMaterial> &adapter) {
      adapter.addBindingPoint("rho", &AnisotropicMaterial::rho);
      adapter.addBindingPoint("c11", &AnisotropicMaterial::c11);
      adapter.addBindingPoint("c12", &AnisotropicMaterial::c12);
      adapter.addBindingPoint("c13", &AnisotropicMaterial::c13);
      adapter.addBindingPoint("c14", &AnisotropicMaterial::c14);
      adapter.addBindingPoint("c15", &AnisotropicMaterial::c15);
      adapter.addBindingPoint("c16", &AnisotropicMaterial::c16);
      adapter.addBindingPoint("c22", &AnisotropicMaterial::c22);
      adapter.addBindingPoint("c23", &AnisotropicMaterial::c23);
      adapter.addBindingPoint("c24", &AnisotropicMaterial::c24);
      adapter.addBindingPoint("c25", &AnisotropicMaterial::c25);
      adapter.addBindingPoint("c26", &AnisotropicMaterial::c26);
      adapter.addBindingPoint("c33", &AnisotropicMaterial::c33);
      adapter.addBindingPoint("c34", &AnisotropicMaterial::c34);
      adapter.addBindingPoint("c35", &AnisotropicMaterial::c35);
      adapter.addBindingPoint("c36", &AnisotropicMaterial::c36);
      adapter.addBindingPoint("c44", &AnisotropicMaterial::c44);
      adapter.addBindingPoint("c45", &AnisotropicMaterial::c45);
      adapter.addBindingPoint("c46", &AnisotropicMaterial::c46);
      adapter.addBindingPoint("c55", &AnisotropicMaterial::c55);
      adapter.addBindingPoint("c56", &AnisotropicMaterial::c56);
      adapter.addBindingPoint("c66", &AnisotropicMaterial::c66);
    }                                                               
    
    template<class T>
    void MaterialParameterDB<T>::evaluateModel(std::string const& fileName, QueryGenerator const * const queryGen) {
      easi::Component* model = loadEasiModel(fileName);
      easi::Query query = queryGen->generate();
      const unsigned numPoints = query.numPoints();

      std::vector<T> materialsFromQuery(numPoints);
      easi::ArrayOfStructsAdapter<T> adapter(materialsFromQuery.data());
      MaterialParameterDB<T>().addBindingPoints(adapter);
      model->evaluate(query, adapter);

      // Only use homogenization when ElementAverageGenerator has been supplied
      if (const ElementAverageGenerator* gen = dynamic_cast<const ElementAverageGenerator*>(queryGen)) {
        const unsigned numElems = numPoints / NUM_QUADPOINTS;
        std::array<double, NUM_QUADPOINTS> quadratureWeights{ gen->getQuadratureWeights() };

        // Compute homogenized material parameters for every element in a specialization for the particular material
        #pragma omp parallel for
        for (unsigned elementIdx = 0; elementIdx < numElems; ++elementIdx) {
          m_materials->at(elementIdx) = this->computeAveragedMaterial(elementIdx, quadratureWeights, materialsFromQuery); 
        }
      } else {
      // Usual behavior without homogenization
        for (unsigned i = 0; i < numPoints; ++i) {
          m_materials->at(i) = T(materialsFromQuery[i]);
        }
      }
      delete model;
    }

    // Computes the averaged material, assuming that materialsFromQuery, stores 
    // NUM_QUADPOINTS material samples per mesh element.
    // We assume that materialsFromQuery[i * NUM_QUADPOINTS, ..., (i+1)*NUM_QUADPOINTS-1] 
    // stores samples from element i.
    // Fallback variant: if no specialization is implemented return materialsFromQuery[i*NUM_QUADPOINTS]
    // Note: this base variant should actually never be called, but only the specializations!
    // The specializations should implement proper averaging.
    template<class T>
    T MaterialParameterDB<T>::computeAveragedMaterial(unsigned elementIdx,
        std::array<double, NUM_QUADPOINTS> const& quadratureWeights,
        std::vector<T> const& materialsFromQuery) {
      logWarning() << "You want me to compute an average material for a generic type. In general, this function should never be called, but always a proper specialization!";
      unsigned globalPointIdx = NUM_QUADPOINTS * elementIdx;
      return T(materialsFromQuery[globalPointIdx]);
    }

    template<>
    ElasticMaterial MaterialParameterDB<ElasticMaterial>::computeAveragedMaterial (unsigned elementIdx, 
     std::array<double, NUM_QUADPOINTS> const& quadratureWeights,
     std::vector<ElasticMaterial> const& materialsFromQuery) {
      double muMeanInv = 0.0;
      double rhoMean = 0.0;
      // Average of v / E with v: Poisson's ratio, E: Young's modulus
      double vERatioMean = 0.0;

      for (unsigned quadPointIdx = 0; quadPointIdx < NUM_QUADPOINTS; ++quadPointIdx) {
        // Divide by volume of reference tetrahedron (1/6)
        const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
        const unsigned globalPointIdx = NUM_QUADPOINTS * elementIdx + quadPointIdx;
        const auto& elementMaterial = materialsFromQuery[globalPointIdx];
        muMeanInv += 1.0 / elementMaterial.mu * quadWeight;
        rhoMean += elementMaterial.rho * quadWeight;
        vERatioMean += elementMaterial.lambda / (2.0 * elementMaterial.mu * (3.0 * elementMaterial.lambda + 2.0 * elementMaterial.mu)) * quadWeight;
      }

      // Harmonic average is used for mu, so take the reciprocal
      double muMean = 1.0 / muMeanInv;
      // Derive lambda from averaged mu and (Poisson ratio / elastic modulus)
      double lambdaMean = (4.0 * std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);

      ElasticMaterial result{};
      result.rho = rhoMean;
      result.mu = muMean;
      result.lambda = lambdaMean;
      return result;
    }

    template<>
    ViscoElasticMaterial MaterialParameterDB<ViscoElasticMaterial>::computeAveragedMaterial (unsigned elementIdx,
        std::array<double, NUM_QUADPOINTS> const& quadratureWeights,
        std::vector<ViscoElasticMaterial> const& materialsFromQuery) {
      double muMeanInv = 0.0;
      double rhoMean = 0.0;
      double vERatioMean = 0.0;
      double QpMean = 0.0;
      double QsMean = 0.0;

      for (unsigned quadPointIdx = 0; quadPointIdx < NUM_QUADPOINTS; ++quadPointIdx) {
        const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
        const unsigned globalPointIdx = NUM_QUADPOINTS * elementIdx + quadPointIdx;
        const auto& elementMaterial = materialsFromQuery[globalPointIdx];
        muMeanInv += 1.0 / elementMaterial.mu * quadWeight;
        rhoMean += elementMaterial.rho * quadWeight;
        vERatioMean += elementMaterial.lambda / (2.0 * elementMaterial.mu * (3.0 * elementMaterial.lambda + 2.0 * elementMaterial.mu)) * quadWeight;
        QpMean += elementMaterial.Qp * quadWeight;
        QsMean += elementMaterial.Qs * quadWeight;
      }

      // Harmonic average is used for mu, so take the reciprocal
      double muMean = 1.0 / muMeanInv;
      // Derive lambda from averaged mu and (Poisson ratio / elastic modulus)
      double lambdaMean = (4.0 * std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);

      ViscoElasticMaterial result{};
      result.rho = rhoMean;
      result.mu = muMean;
      result.lambda = lambdaMean;
      result.Qp = QpMean;
      result.Qs = QsMean;
      
      return result;
    }
    
    template<>
    void MaterialParameterDB<AnisotropicMaterial>::evaluateModel(std::string const& fileName, QueryGenerator const * const queryGen) {
      easi::Component* model = loadEasiModel(fileName);
      easi::Query query = queryGen->generate();
      auto suppliedParameters = model->suppliedParameters();
      //TODO(Sebastian): inhomogeneous materials, where in some parts only mu and lambda are given
      //                 and in other parts the full elastic tensor is given

      //if we look for an anisotropic material and only mu and lambda are supplied, 
      //assume isotropic behavior and calculate the parameters accordingly
      if (suppliedParameters.find("mu") != suppliedParameters.end() && suppliedParameters.find("lambda") != suppliedParameters.end()) {
        std::vector<ElasticMaterial> elasticMaterials(query.numPoints());
        easi::ArrayOfStructsAdapter<ElasticMaterial> adapter(elasticMaterials.data());
        MaterialParameterDB<ElasticMaterial>().addBindingPoints(adapter);
        unsigned numPoints = query.numPoints();
        model->evaluate(query, adapter);

        for(unsigned i = 0; i < numPoints; i++) {
          m_materials->at(i) = AnisotropicMaterial(elasticMaterials[i]);
        }
      }
      else {
        easi::ArrayOfStructsAdapter<AnisotropicMaterial> arrayOfStructsAdapter(m_materials->data());
        addBindingPoints(arrayOfStructsAdapter);
        model->evaluate(query, arrayOfStructsAdapter);  
        
      }
      delete model;
    }

    void FaultParameterDB::evaluateModel(std::string const& fileName, QueryGenerator const * const queryGen) {
      easi::Component* model = loadEasiModel(fileName);
      easi::Query query = queryGen->generate();

      easi::ArraysAdapter<real> adapter;
      for (auto& kv : m_parameters) {
        adapter.addBindingPoint(kv.first, kv.second.first, kv.second.second);
      }
      model->evaluate(query, adapter); 
      
      delete model;
    }

  }
}

std::set<std::string> seissol::initializers::FaultParameterDB::faultProvides(std::string const& fileName) {
  if (fileName.length() == 0) {
    return std::set<std::string>();
  }
  easi::Component* model = loadEasiModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;
  return supplied;
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

namespace seissol::initializers {
QueryGenerator* getBestQueryGenerator(bool anelasticity,
    bool plasticity,
    bool anisotropy,
    bool poroelasticity,
    bool useCellHomogenizedMaterial,
    MeshReader const& meshReader) {
  QueryGenerator* queryGen = nullptr;
  if (!useCellHomogenizedMaterial) {
    queryGen = new ElementBarycentreGenerator(meshReader);
  } else {
    if (anisotropy) {
      logWarning() << "Material Averaging is not implemented for anisotropic materials. Falling back to material properties sampled from the element barycenters instead.";
      queryGen = new ElementBarycentreGenerator(meshReader);
    } else if (plasticity) {
      logWarning() << "Material Averaging is not implemented for plastic materials. Falling back to material properties sampled from the element barycenters instead.";
      queryGen = new ElementBarycentreGenerator(meshReader);
    } else if (poroelasticity) {
      logWarning() << "Material Averaging is not implemented for poroelastic materials. Falling back to material properties sampled from the element barycenters instead.";
      queryGen = new ElementBarycentreGenerator(meshReader);
    } else {
      queryGen = new ElementAverageGenerator(meshReader);
    }
  }
  return queryGen;
}
} // namespace seissol::initializers
    

template class seissol::initializers::MaterialParameterDB<seissol::model::AnisotropicMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::ElasticMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::ViscoElasticMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::PoroElasticMaterial>;
template class seissol::initializers::MaterialParameterDB<seissol::model::Plasticity>;


