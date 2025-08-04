// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include <Common/Constants.h>
#include <Equations/Datastructures.h>
#include <Equations/acoustic/Model/Datastructures.h>
#include <Equations/anisotropic/Model/Datastructures.h>
#include <Equations/elastic/Model/Datastructures.h>
#include <Equations/poroelastic/Model/Datastructures.h>
#include <Equations/viscoelastic2/Model/Datastructures.h>
#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshTools.h>
#include <Kernels/Precision.h>
#include <Model/CommonDatastructures.h>
#include <Model/Datastructures.h>
#include <Solver/MultipleSimulations.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <easi/Query.h>
#include <init.h>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tensor.h>
#include <utility>
#include <vector>
#ifdef USE_HDF
// PUML.h needs to be included before Downward.h

#include "PUML/PUML.h"

#include "PUML/Downward.h"

#include "generated_code/kernel.h"
#endif
#include "ParameterDB.h"
#include <algorithm>
#include <cmath>

#include "DynamicRupture/Misc.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include "easi/ResultAdapter.h"
#include "easi/YAMLParser.h"
#ifdef USE_ASAGI
#include "Reader/AsagiReader.h"
#endif
#include "utils/logger.h"

namespace seissol::initializer {

CellToVertexArray::CellToVertexArray(size_t size,
                                     const CellToVertexFunction& elementCoordinates,
                                     const CellToGroupFunction& elementGroups)
    : size(size), elementCoordinates(elementCoordinates), elementGroups(elementGroups) {}

CellToVertexArray
    CellToVertexArray::fromMeshReader(const seissol::geometry::MeshReader& meshReader) {
  const auto& elements = meshReader.getElements();
  const auto& vertices = meshReader.getVertices();

  return CellToVertexArray(
      elements.size(),
      [&](size_t index) {
        std::array<Eigen::Vector3d, 4> verts;
        for (size_t i = 0; i < Cell::NumVertices; ++i) {
          auto vindex = elements[index].vertices[i];
          const auto& vertex = vertices[vindex];
          verts[i] << vertex.coords[0], vertex.coords[1], vertex.coords[2];
        }
        return verts;
      },
      [&](size_t index) { return elements[index].group; });
}

#ifdef USE_HDF
CellToVertexArray CellToVertexArray::fromPUML(const PUML::TETPUML& mesh) {
  const int* groups = reinterpret_cast<const int*>(mesh.cellData(0));
  const auto& elements = mesh.cells();
  const auto& vertices = mesh.vertices();
  return CellToVertexArray(
      elements.size(),
      [&](size_t cell) {
        std::array<Eigen::Vector3d, 4> x;
        unsigned vertLids[Cell::NumVertices];
        PUML::Downward::vertices(mesh, elements[cell], vertLids);
        for (std::size_t vtx = 0; vtx < Cell::NumVertices; ++vtx) {
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            x[vtx](d) = vertices[vertLids[vtx]].coordinate()[d];
          }
        }
        return x;
      },
      [groups](size_t cell) { return groups[cell]; });
}
#endif

CellToVertexArray CellToVertexArray::fromVectors(
    const std::vector<std::array<std::array<double, Cell::Dim>, Cell::NumVertices>>& vertices,
    const std::vector<int>& groups) {
  assert(vertices.size() == groups.size());

  return CellToVertexArray(
      vertices.size(),
      [&](size_t idx) {
        std::array<Eigen::Vector3d, Cell::NumVertices> verts;
        for (size_t i = 0; i < Cell::NumVertices; ++i) {
          verts[i] << vertices[idx][i][0], vertices[idx][i][1], vertices[idx][i][2];
        }
        return verts;
      },
      [&](size_t i) { return groups[i]; });
}

CellToVertexArray CellToVertexArray::join(std::vector<CellToVertexArray> arrays) {
  std::size_t totalSize = 0;
  std::vector<std::size_t> sizes(arrays.size());
  std::vector<std::size_t> offsets(arrays.size());
  for (std::size_t i = 0; i < arrays.size(); ++i) {
    offsets[i] = totalSize;
    totalSize += arrays[i].size;
    sizes[i] = totalSize;
  }
  return CellToVertexArray(
      totalSize,
      [=](size_t idx) {
        for (std::size_t i = 0; i < sizes.size(); ++i) {
          if (idx < sizes[i]) {
            return arrays[i].elementCoordinates(idx - offsets[i]);
          }
        }
        throw std::out_of_range(std::to_string(idx) + " vs " + std::to_string(totalSize));
      },
      [=](size_t idx) {
        for (std::size_t i = 0; i < sizes.size(); ++i) {
          if (idx < sizes[i]) {
            return arrays[i].elementGroups(idx - offsets[i]);
          }
        }
        throw std::out_of_range(std::to_string(idx) + " vs " + std::to_string(totalSize));
      });
}

easi::Query ElementBarycenterGenerator::generate() const {
  easi::Query query(m_cellToVertex.size, Cell::Dim);

#pragma omp parallel for schedule(static)
  for (unsigned elem = 0; elem < m_cellToVertex.size; ++elem) {
    auto vertices = m_cellToVertex.elementCoordinates(elem);
    Eigen::Vector3d barycenter = (vertices[0] + vertices[1] + vertices[2] + vertices[3]) * 0.25;
    query.x(elem, 0) = barycenter(0);
    query.x(elem, 1) = barycenter(1);
    query.x(elem, 2) = barycenter(2);
    query.group(elem) = m_cellToVertex.elementGroups(elem);
  }
  return query;
}

ElementAverageGenerator::ElementAverageGenerator(const CellToVertexArray& cellToVertex)
    : m_cellToVertex(cellToVertex) {
  double quadraturePoints[NumQuadpoints][3];
  double quadratureWeights[NumQuadpoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, ConvergenceOrder);

  std::copy(
      std::begin(quadratureWeights), std::end(quadratureWeights), std::begin(m_quadratureWeights));
  for (std::size_t i = 0; i < NumQuadpoints; ++i) {
    std::copy(std::begin(quadraturePoints[i]),
              std::end(quadraturePoints[i]),
              std::begin(m_quadraturePoints[i]));
  }
}

easi::Query ElementAverageGenerator::generate() const {
  // Generate query using quadrature points for each element
  easi::Query query(m_cellToVertex.size * NumQuadpoints, 3);

// Transform quadrature points to global coordinates for all elements
#pragma omp parallel for schedule(static) collapse(2)
  for (unsigned elem = 0; elem < m_cellToVertex.size; ++elem) {
    for (unsigned i = 0; i < NumQuadpoints; ++i) {
      auto vertices = m_cellToVertex.elementCoordinates(elem);
      Eigen::Vector3d transformed = seissol::transformations::tetrahedronReferenceToGlobal(
          vertices[0], vertices[1], vertices[2], vertices[3], m_quadraturePoints[i].data());
      query.x(elem * NumQuadpoints + i, 0) = transformed(0);
      query.x(elem * NumQuadpoints + i, 1) = transformed(1);
      query.x(elem * NumQuadpoints + i, 2) = transformed(2);
      query.group(elem * NumQuadpoints + i) = m_cellToVertex.elementGroups(elem);
    }
  }

  return query;
}

easi::Query FaultBarycenterGenerator::generate() const {
  const std::vector<Fault>& fault = m_meshReader.getFault();
  const std::vector<Element>& elements = m_meshReader.getElements();
  const std::vector<Vertex>& vertices = m_meshReader.getVertices();

  easi::Query query(m_numberOfPoints * fault.size(), Cell::Dim);
  unsigned q = 0;
  for (const Fault& f : fault) {
    int element = 0;
    int side = 0;
    if (f.element >= 0) {
      element = f.element;
      side = f.side;
    } else {
      element = f.neighborElement;
      side = f.neighborSide;
    }

    double barycenter[3] = {0.0, 0.0, 0.0};
    MeshTools::center(elements[element], side, vertices, barycenter);
    for (unsigned n = 0; n < m_numberOfPoints; ++n, ++q) {
      for (unsigned dim = 0; dim < Cell::Dim; ++dim) {
        query.x(q, dim) = barycenter[dim];
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

easi::Query FaultGPGenerator::generate() const {
  const std::vector<Fault>& fault = m_meshReader.getFault();
  const std::vector<Element>& elements = m_meshReader.getElements();
  auto cellToVertex = CellToVertexArray::fromMeshReader(m_meshReader);

  constexpr size_t NumPoints = dr::misc::NumPaddedPointsSingleSim;
  auto pointsView = init::quadpoints::view::create(const_cast<real*>(init::quadpoints::Values));
  easi::Query query(NumPoints * m_faceIDs.size(), Cell::Dim);
  unsigned q = 0;
  // loop over all fault elements which are managed by this generator
  // note: we have one generator per LTS layer
  for (const unsigned faultId : m_faceIDs) {
    const Fault& f = fault.at(faultId);
    int element = 0;
    int side = 0;
    int sideOrientation = 0;
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
    for (unsigned n = 0; n < NumPoints; ++n, ++q) {
      double xiEtaZeta[3];
      double localPoints[2] = {seissol::multisim::multisimTranspose(pointsView, n, 0),
                               seissol::multisim::multisimTranspose(pointsView, n, 1)};
      // padded points are in the middle of the tetrahedron
      if (n >= dr::misc::NumBoundaryGaussPoints) {
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

using namespace seissol::model;

template <class T>
void MaterialParameterDB<T>::evaluateModel(const std::string& fileName,
                                           const QueryGenerator& queryGen) {
  easi::Component* model = loadEasiModel(fileName);
  auto suppliedParameters = model->suppliedParameters();
  for (const auto& [name, pointer] : T::ParameterMap) {
  }

  easi::Query query = queryGen.generate();
  const std::size_t numPoints = query.numPoints();

  std::vector<T> materialsFromQuery(numPoints);
  easi::ArrayOfStructsAdapter<T> adapter(materialsFromQuery.data());
  for (const auto& [name, pointer] : T::ParameterMap) {
    adapter.addBindingPoint(name, pointer);
  }
  model->evaluate(query, adapter);

  // Only use homogenization when ElementAverageGenerator has been supplied
  if (const auto* gen = dynamic_cast<const ElementAverageGenerator*>(&queryGen)) {
    const std::size_t numElems = numPoints / NumQuadpoints;
    const std::vector<double> quadratureWeights(gen->getQuadratureWeights().begin(),
                                                gen->getQuadratureWeights().end());

// Compute homogenized material parameters for every element in a specialization for the
// particular material
#pragma omp parallel for schedule(static)
    for (std::size_t elementIdx = 0; elementIdx < numElems; ++elementIdx) {
      m_materials->at(elementIdx) = MaterialAverager<T>::computeAveragedMaterial(
          elementIdx, quadratureWeights, materialsFromQuery);
    }
  } else {
    // Usual behavior without homogenization
    for (std::size_t i = 0; i < numPoints; ++i) {
      m_materials->at(i) = T(materialsFromQuery[i]);
    }
  }
  delete model;
}

template <>
void MaterialParameterDB<AnisotropicMaterial>::evaluateModel(const std::string& fileName,
                                                             const QueryGenerator& queryGen) {
  easi::Component* model = loadEasiModel(fileName);
  easi::Query query = queryGen.generate();
  auto suppliedParameters = model->suppliedParameters();
  // TODO(Sebastian): inhomogeneous materials, where in some parts only mu and lambda are given
  //                  and in other parts the full elastic tensor is given

  // if we look for an anisotropic material and only mu and lambda are supplied,
  // assume isotropic behavior and calculate the parameters accordingly
  if (suppliedParameters.find("mu") != suppliedParameters.end() &&
      suppliedParameters.find("lambda") != suppliedParameters.end()) {
    MaterialParameterDB<ElasticMaterial> edb;
    std::vector<ElasticMaterial> preMaterials(m_materials->size());
    edb.setMaterialVector(&preMaterials);
    edb.evaluateModel(fileName, queryGen);

    for (std::size_t i = 0; i < m_materials->size(); i++) {
      m_materials->at(i) = AnisotropicMaterial(preMaterials[i]);
    }
  } else {
    easi::ArrayOfStructsAdapter<AnisotropicMaterial> arrayOfStructsAdapter(m_materials->data());
    for (const auto& [name, pointer] : AnisotropicMaterial::ParameterMap) {
      arrayOfStructsAdapter.addBindingPoint(name, pointer);
    }
    model->evaluate(query, arrayOfStructsAdapter);
  }
  delete model;
}

// Computes the averaged material, assuming that materialsFromQuery, stores
// NUM_QUADPOINTS material samples per mesh element.
// We assume that materialsFromQuery[i * NUM_QUADPOINTS, ..., (i+1)*NUM_QUADPOINTS-1]
// stores samples from element i.

template <>
struct MaterialAverager<ElasticMaterial> {
  static constexpr bool Implemented = true;
  static ElasticMaterial
      computeAveragedMaterial(std::size_t elementIdx,
                              const std::vector<double>& quadratureWeights,
                              const std::vector<ElasticMaterial>& materialsFromQuery) {
    double muMeanInv = 0.0;
    double rhoMean = 0.0;
    // Average of v / E with v: Poisson's ratio, E: Young's modulus
    double vERatioMean = 0.0;

    // Acoustic material has zero mu. This is a special case because the harmonic mean of a set
    // of numbers that includes zero is defined as zero.
    // Hence: If part of the element is acoustic, the entire element is considered to be acoustic!
    bool isAcoustic = false;

    // Average of the bulk modulus, used for acoustic material
    double kMeanInv = 0.0;

    for (unsigned quadPointIdx = 0; quadPointIdx < NumQuadpoints; ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const unsigned globalPointIdx = NumQuadpoints * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery[globalPointIdx];
      isAcoustic |= elementMaterial.mu == 0.0;
      if (!isAcoustic) {
        muMeanInv += 1.0 / elementMaterial.mu * quadWeight;
      }
      rhoMean += elementMaterial.rho * quadWeight;
      vERatioMean +=
          elementMaterial.lambda /
          (2.0 * elementMaterial.mu * (3.0 * elementMaterial.lambda + 2.0 * elementMaterial.mu)) *
          quadWeight;
      kMeanInv += 1.0 / (elementMaterial.lambda + (2.0 / 3.0) * elementMaterial.mu) * quadWeight;
    }

    ElasticMaterial result{};
    result.rho = rhoMean;

    // Harmonic average is used for mu/K, so take the reciprocal
    if (isAcoustic) {
      result.lambda = 1.0 / kMeanInv;
      result.mu = 0.0;
    } else {
      const auto muMean = 1.0 / muMeanInv;
      // Derive lambda from averaged mu and (Poisson ratio / elastic modulus)
      result.lambda =
          (4.0 * std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);
      result.mu = muMean;
    }

    return result;
  }
};

template <>
struct MaterialAverager<ViscoElasticMaterial> {
  static constexpr bool Implemented = true;
  static ViscoElasticMaterial
      computeAveragedMaterial(std::size_t elementIdx,
                              const std::vector<double>& quadratureWeights,
                              const std::vector<ViscoElasticMaterial>& materialsFromQuery) {
    double muMeanInv = 0.0;
    double rhoMean = 0.0;
    double vERatioMean = 0.0;
    double qpMean = 0.0;
    double qsMean = 0.0;

    for (unsigned quadPointIdx = 0; quadPointIdx < NumQuadpoints; ++quadPointIdx) {
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const unsigned globalPointIdx = NumQuadpoints * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery[globalPointIdx];
      muMeanInv += 1.0 / elementMaterial.mu * quadWeight;
      rhoMean += elementMaterial.rho * quadWeight;
      vERatioMean +=
          elementMaterial.lambda /
          (2.0 * elementMaterial.mu * (3.0 * elementMaterial.lambda + 2.0 * elementMaterial.mu)) *
          quadWeight;
      qpMean += elementMaterial.Qp * quadWeight;
      qsMean += elementMaterial.Qs * quadWeight;
    }

    // Harmonic average is used for mu, so take the reciprocal
    const double muMean = 1.0 / muMeanInv;
    // Derive lambda from averaged mu and (Poisson ratio / elastic modulus)
    const double lambdaMean =
        (4.0 * std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);

    ViscoElasticMaterial result{};
    result.rho = rhoMean;
    result.mu = muMean;
    result.lambda = lambdaMean;
    result.Qp = qpMean;
    result.Qs = qsMean;

    return result;
  }
};

void FaultParameterDB::evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) {
  easi::Component* model = loadEasiModel(fileName);
  easi::Query query = queryGen.generate();

  easi::ArraysAdapter<real> adapter;
  for (auto& kv : m_parameters) {
    adapter.addBindingPoint(
        kv.first, kv.second.first + simid, kv.second.second * multisim::NumSimulations);
  }
  model->evaluate(query, adapter);

  delete model;
}

std::set<std::string> FaultParameterDB::faultProvides(const std::string& fileName) {
  if (fileName.empty()) {
    return {};
  }
  easi::Component* model = loadEasiModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;
  return supplied;
}

EasiBoundary::EasiBoundary(const std::string& fileName) : model(loadEasiModel(fileName)) {}

EasiBoundary::EasiBoundary(EasiBoundary&& other) noexcept : model(other.model) {}

EasiBoundary& EasiBoundary::operator=(EasiBoundary&& other) noexcept {
  std::swap(model, other.model);
  return *this;
}

EasiBoundary::~EasiBoundary() { delete model; }

void EasiBoundary::query(const real* nodes, real* mapTermsData, real* constantTermsData) const {
  if (model == nullptr) {
    logError() << "Model for easi-provided boundary is not initialized.";
  }
  if (tensor::INodal::Shape[1] != 9) {
    logError() << "easi-provided boundary data is only supported for elastic material at the "
                  "moment currently.";
  }
  assert(mapTermsData != nullptr);
  assert(constantTermsData != nullptr);
  constexpr auto NumNodes = tensor::INodal::Shape[0];
  auto query = easi::Query{NumNodes, 3};
  size_t offset{0};
  for (unsigned i = 0; i < NumNodes; ++i) {
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
        for (size_t k = 0; k < mapTerms.shape(2); ++k) {
          mapTerms(i, j, k) = (varName == otherVarName) ? 1.0 : 0.0;
        }
      }
      ++offset;
    }
  }
  model->evaluate(query, adapter);
}

easi::Component* loadEasiModel(const std::string& fileName) {
#ifdef USE_ASAGI
  seissol::asagi::AsagiReader asagiReader;
  easi::YAMLParser parser(3, &asagiReader);
#else
  easi::YAMLParser parser(3);
#endif
  return parser.parse(fileName);
}

std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool plasticity,
                                                      bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex) {
  std::shared_ptr<QueryGenerator> queryGen;
  if (!useCellHomogenizedMaterial) {
    queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
  } else {
    if (!MaterialAverager<MaterialT>::Implemented) {
      logWarning() << "Material Averaging is not implemented for " << MaterialT::Text
                   << " materials. Falling back to "
                      "material properties sampled from the element barycenters instead.";
      queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
    } else if (plasticity) {
      logWarning()
          << "Material Averaging is not implemented for plastic materials. Falling back to "
             "material properties sampled from the element barycenters instead.";
      queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
    } else {
      queryGen = std::make_shared<ElementAverageGenerator>(cellToVertex);
    }
  }
  return queryGen;
}

template class MaterialParameterDB<seissol::model::AnisotropicMaterial>;
template class MaterialParameterDB<seissol::model::ElasticMaterial>;
template class MaterialParameterDB<seissol::model::AcousticMaterial>;
template class MaterialParameterDB<seissol::model::ViscoElasticMaterial>;
template class MaterialParameterDB<seissol::model::PoroElasticMaterial>;
template class MaterialParameterDB<seissol::model::Plasticity>;

} // namespace seissol::initializer
