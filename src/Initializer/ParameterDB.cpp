// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Equations/acoustic/Model/Datastructures.h"
#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/viscoelastic2/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Geometry/PUMLReader.h"
#include "Kernels/Precision.h"
#include "Model/CommonDatastructures.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <cstddef>
#include <easi/Component.h>
#include <easi/Query.h>
#include <exception>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
#ifdef USE_HDF
// PUML.h needs to be included before Downward.h

#include "GeneratedCode/kernel.h"

#include <PUML/Downward.h>
#endif
#include "DynamicRupture/Misc.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "ParameterDB.h"
#include "SeisSol.h"
#include "easi/ResultAdapter.h"
#include "easi/YAMLParser.h"

#include <algorithm>
#include <cmath>
#ifdef USE_ASAGI
#include "Reader/AsagiReader.h"
#endif
#include <utils/logger.h>

using namespace seissol::model;

namespace seissol::initializer {

namespace {

void easiEvalSafe(easi::Component* model,
                  easi::Query& query,
                  easi::ResultAdapter& adapter,
                  const std::string& hint) {
  try {
    model->evaluate(query, adapter);
  } catch (const std::exception& error) {
    logError() << "Error while evaluating an easi model for" << hint.c_str() << ":"
               << std::string(error.what());
  }
}

} // namespace

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
CellToVertexArray CellToVertexArray::fromPUML(const seissol::geometry::PumlMesh& mesh) {
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

easi::Query ElementBarycenterGenerator::generate() const {
  easi::Query query(cellToVertex_.size, Cell::Dim);

#pragma omp parallel for schedule(static)
  for (unsigned elem = 0; elem < cellToVertex_.size; ++elem) {
    auto vertices = cellToVertex_.elementCoordinates(elem);
    Eigen::Vector3d barycenter = (vertices[0] + vertices[1] + vertices[2] + vertices[3]) * 0.25;
    query.x(elem, 0) = barycenter(0);
    query.x(elem, 1) = barycenter(1);
    query.x(elem, 2) = barycenter(2);
    query.group(elem) = cellToVertex_.elementGroups(elem);
  }
  return query;
}

ElementAverageGenerator::ElementAverageGenerator(const CellToVertexArray& cellToVertex)
    : cellToVertex_(cellToVertex) {
  double quadraturePoints[NumQuadpoints][3];
  double quadratureWeights[NumQuadpoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, ConvergenceOrder);

  std::copy(
      std::begin(quadratureWeights), std::end(quadratureWeights), std::begin(quadratureWeights_));
  for (std::size_t i = 0; i < NumQuadpoints; ++i) {
    std::copy(std::begin(quadraturePoints[i]),
              std::end(quadraturePoints[i]),
              std::begin(quadraturePoints_[i]));
  }
}

easi::Query ElementAverageGenerator::generate() const {
  // Generate query using quadrature points for each element
  easi::Query query(cellToVertex_.size * NumQuadpoints, 3);

// Transform quadrature points to global coordinates for all elements
#pragma omp parallel for schedule(static) collapse(2)
  for (unsigned elem = 0; elem < cellToVertex_.size; ++elem) {
    for (unsigned i = 0; i < NumQuadpoints; ++i) {
      auto vertices = cellToVertex_.elementCoordinates(elem);
      Eigen::Vector3d transformed = seissol::transformations::tetrahedronReferenceToGlobal(
          vertices[0], vertices[1], vertices[2], vertices[3], quadraturePoints_[i].data());
      query.x(elem * NumQuadpoints + i, 0) = transformed(0);
      query.x(elem * NumQuadpoints + i, 1) = transformed(1);
      query.x(elem * NumQuadpoints + i, 2) = transformed(2);
      query.group(elem * NumQuadpoints + i) = cellToVertex_.elementGroups(elem);
    }
  }

  return query;
}

easi::Query FaultBarycenterGenerator::generate() const {
  const std::vector<Fault>& fault = meshReader_.getFault();
  const std::vector<Element>& elements = meshReader_.getElements();
  const std::vector<Vertex>& vertices = meshReader_.getVertices();

  easi::Query query(numberOfPoints_ * fault.size(), Cell::Dim);
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
    for (unsigned n = 0; n < numberOfPoints_; ++n, ++q) {
      for (unsigned dim = 0; dim < Cell::Dim; ++dim) {
        query.x(q, dim) = barycenter[dim];
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

easi::Query FaultGPGenerator::generate() const {
  const std::vector<Fault>& fault = meshReader_.getFault();
  const std::vector<Element>& elements = meshReader_.getElements();
  auto cellToVertex = CellToVertexArray::fromMeshReader(meshReader_);

  constexpr size_t NumPoints = dr::misc::NumPaddedPointsSingleSim;
  auto pointsView = init::quadpoints::view::create(const_cast<real*>(init::quadpoints::Values));
  easi::Query query(NumPoints * faceIDs_.size(), Cell::Dim);
  unsigned q = 0;
  // loop over all fault elements which are managed by this generator
  // note: we have one generator per LTS layer
  for (const unsigned faultId : faceIDs_) {
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

template <>
void MaterialParameterDB<ElasticMaterial>::addBindingPoints(
    easi::ArrayOfStructsAdapter<ElasticMaterial>& adapter) {
  adapter.addBindingPoint("rho", &ElasticMaterial::rho);
  adapter.addBindingPoint("mu", &ElasticMaterial::mu);
  adapter.addBindingPoint("lambda", &ElasticMaterial::lambda);
}

template <>
void MaterialParameterDB<AcousticMaterial>::addBindingPoints(
    easi::ArrayOfStructsAdapter<AcousticMaterial>& adapter) {
  adapter.addBindingPoint("rho", &AcousticMaterial::rho);
  adapter.addBindingPoint("lambda", &AcousticMaterial::lambda);
}

template <>
void MaterialParameterDB<ViscoElasticMaterial>::addBindingPoints(
    easi::ArrayOfStructsAdapter<ViscoElasticMaterial>& adapter) {
  adapter.addBindingPoint("rho", &ViscoElasticMaterial::rho);
  adapter.addBindingPoint("mu", &ViscoElasticMaterial::mu);
  adapter.addBindingPoint("lambda", &ViscoElasticMaterial::lambda);
  adapter.addBindingPoint("Qp", &ViscoElasticMaterial::qp);
  adapter.addBindingPoint("Qs", &ViscoElasticMaterial::qs);
}

template <>
void MaterialParameterDB<PoroElasticMaterial>::addBindingPoints(
    easi::ArrayOfStructsAdapter<PoroElasticMaterial>& adapter) {
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

template <>
void MaterialParameterDB<Plasticity>::addBindingPoints(
    easi::ArrayOfStructsAdapter<Plasticity>& adapter) {
  adapter.addBindingPoint("bulkFriction", &Plasticity::bulkFriction);
  adapter.addBindingPoint("plastCo", &Plasticity::plastCo);
  adapter.addBindingPoint("s_xx", &Plasticity::sXX);
  adapter.addBindingPoint("s_yy", &Plasticity::sYY);
  adapter.addBindingPoint("s_zz", &Plasticity::sZZ);
  adapter.addBindingPoint("s_xy", &Plasticity::sXY);
  adapter.addBindingPoint("s_yz", &Plasticity::sYZ);
  adapter.addBindingPoint("s_xz", &Plasticity::sXZ);
}

template <>
void MaterialParameterDB<AnisotropicMaterial>::addBindingPoints(
    easi::ArrayOfStructsAdapter<AnisotropicMaterial>& adapter) {
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

template <class T>
void MaterialParameterDB<T>::evaluateModel(const std::string& fileName,
                                           const QueryGenerator& queryGen) {
  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  easi::Query query = queryGen.generate();
  const unsigned numPoints = query.numPoints();

  std::vector<T> materialsFromQuery(numPoints);
  easi::ArrayOfStructsAdapter<T> adapter(materialsFromQuery.data());
  MaterialParameterDB<T>().addBindingPoints(adapter);

  easiEvalSafe(model, query, adapter, "volume material");

  // Only use homogenization when ElementAverageGenerator has been supplied
  if (const auto* gen = dynamic_cast<const ElementAverageGenerator*>(&queryGen)) {
    const unsigned numElems = numPoints / NumQuadpoints;
    const std::array<double, NumQuadpoints> quadratureWeights{gen->getQuadratureWeights()};

// Compute homogenized material parameters for every element in a specialization for the
// particular material
#pragma omp parallel for
    for (unsigned elementIdx = 0; elementIdx < numElems; ++elementIdx) {
      materials_->at(elementIdx) =
          this->computeAveragedMaterial(elementIdx, quadratureWeights, materialsFromQuery);
    }
  } else {
    // Usual behavior without homogenization
    for (unsigned i = 0; i < numPoints; ++i) {
      materials_->at(i) = T(materialsFromQuery[i]);
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
template <class T>
T MaterialParameterDB<T>::computeAveragedMaterial(
    unsigned elementIdx,
    const std::array<double, NumQuadpoints>& /*quadratureWeights*/,
    const std::vector<T>& materialsFromQuery) {
  logWarning() << "You want me to compute an average material for a generic type. In general, this "
                  "function should never be called, but always a proper specialization!";
  const unsigned globalPointIdx = NumQuadpoints * elementIdx;
  return T(materialsFromQuery[globalPointIdx]);
}

template <>
ElasticMaterial MaterialParameterDB<ElasticMaterial>::computeAveragedMaterial(
    unsigned elementIdx,
    const std::array<double, NumQuadpoints>& quadratureWeights,
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
    result.lambda = (4.0 * std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);
    result.mu = muMean;
  }

  return result;
}

template <>
ViscoElasticMaterial MaterialParameterDB<ViscoElasticMaterial>::computeAveragedMaterial(
    unsigned elementIdx,
    const std::array<double, NumQuadpoints>& quadratureWeights,
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
    qpMean += elementMaterial.qp * quadWeight;
    qsMean += elementMaterial.qs * quadWeight;
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
  result.qp = qpMean;
  result.qs = qsMean;

  return result;
}

template <>
void MaterialParameterDB<AnisotropicMaterial>::evaluateModel(const std::string& fileName,
                                                             const QueryGenerator& queryGen) {
  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  easi::Query query = queryGen.generate();
  auto suppliedParameters = model->suppliedParameters();
  // TODO(Sebastian): inhomogeneous materials, where in some parts only mu and lambda are given
  //                  and in other parts the full elastic tensor is given

  // if we look for an anisotropic material and only mu and lambda are supplied,
  // assume isotropic behavior and calculate the parameters accordingly
  if (suppliedParameters.find("mu") != suppliedParameters.end() &&
      suppliedParameters.find("lambda") != suppliedParameters.end()) {
    std::vector<ElasticMaterial> elasticMaterials(query.numPoints());
    easi::ArrayOfStructsAdapter<ElasticMaterial> adapter(elasticMaterials.data());
    MaterialParameterDB<ElasticMaterial>().addBindingPoints(adapter);
    const unsigned numPoints = query.numPoints();
    easiEvalSafe(model, query, adapter, "volume material (anisotropic -> elastic)");

    for (unsigned i = 0; i < numPoints; i++) {
      materials_->at(i) = AnisotropicMaterial(elasticMaterials[i]);
    }
  } else {
    easi::ArrayOfStructsAdapter<AnisotropicMaterial> arrayOfStructsAdapter(materials_->data());
    addBindingPoints(arrayOfStructsAdapter);
    easiEvalSafe(model, query, arrayOfStructsAdapter, "volume material (anisotropic)");
  }
  delete model;
}

void FaultParameterDB::evaluateModel(const std::string& fileName, const QueryGenerator& queryGen) {
  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  easi::Query query = queryGen.generate();

  easi::ArraysAdapter<real> adapter;
  for (auto& kv : parameters_) {
    adapter.addBindingPoint(
        kv.first, kv.second.first + simid_, kv.second.second * multisim::NumSimulations);
  }

  easiEvalSafe(model, query, adapter, "fault material");

  delete model;
}

std::set<std::string> FaultParameterDB::faultProvides(const std::string& fileName) {
  if (fileName.empty()) {
    return {};
  }

  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;
  return supplied;
}

EasiBoundary::EasiBoundary(const std::string& fileName) : model_(loadEasiModel(fileName)) {}

EasiBoundary::EasiBoundary(EasiBoundary&& other) noexcept : model_(other.model_) {}

EasiBoundary& EasiBoundary::operator=(EasiBoundary&& other) noexcept {
  std::swap(model_, other.model_);
  return *this;
}

EasiBoundary::~EasiBoundary() { delete model_; }

void EasiBoundary::query(const real* nodes, real* mapTermsData, real* constantTermsData) const {
  if (model_ == nullptr) {
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
  const auto& supplied = model_->suppliedParameters();

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
  easiEvalSafe(model_, query, adapter, "Dirichlet BC data");
}

easi::Component* loadEasiModel(const std::string& fileName) {
#ifdef USE_ASAGI
  seissol::asagi::AsagiReader asagiReader;
  easi::YAMLParser parser(3, &asagiReader);
#else
  easi::YAMLParser parser(3);
#endif
  try {
    return parser.parse(fileName);
  } catch (const std::exception& error) {
    logError() << "Error while parsing easi file" << fileName << ":" << std::string(error.what());
    // silence no-return warnings
    return nullptr;
  }
}

std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool plasticity,
                                                      bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex) {
  std::shared_ptr<QueryGenerator> queryGen;
  if (!useCellHomogenizedMaterial) {
    queryGen = std::make_shared<ElementBarycenterGenerator>(cellToVertex);
  } else {
    if (MaterialT::Type != MaterialType::Viscoelastic && MaterialT::Type != MaterialType::Elastic) {
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
