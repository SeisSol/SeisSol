// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "ParameterDB.h"

#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
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
#include "Model/Plasticity.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"
#include "easi/ResultAdapter.h"
#include "easi/YAMLParser.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <easi/Component.h>
#include <easi/Query.h>
#include <exception>
#include <iterator>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef USE_HDF
// PUML.h needs to be included before Downward.h

#include "GeneratedCode/kernel.h"

#include <PUML/Downward.h>
#endif

#ifdef USE_ASAGI
#include "Reader/AsagiReader.h"
#endif

using namespace seissol::model;

namespace seissol::initializer {

namespace {

template <typename SurrogateMaterialT>
bool canEvaluateFor(const std::set<std::string>& parameters) {
  return std::all_of(
      SurrogateMaterialT::ParameterMap.begin(),
      SurrogateMaterialT::ParameterMap.end(),
      [&](const auto& ptr) { return parameters.find(ptr.first) != parameters.end(); });
}

template <typename MaterialT, typename SurrogateMaterialT, typename... SurrogatesT>
bool surrogateEvaluate(const std::string& fileName,
                       const QueryGenerator& queryGen,
                       std::vector<MaterialT>* materials,
                       const std::set<std::string>& parameters) {
  if constexpr (std::is_constructible_v<MaterialT, SurrogateMaterialT>) {
    if (canEvaluateFor<SurrogateMaterialT>(parameters)) {
      MaterialParameterDB<SurrogateMaterialT> edb;
      std::vector<SurrogateMaterialT> preMaterials(materials->size());
      edb.setMaterialVector(&preMaterials);
      edb.evaluateModel(fileName, queryGen);

      for (std::size_t i = 0; i < materials->size(); i++) {
        materials->at(i) = MaterialT(preMaterials[i]);
      }

      return true;
    }
  }
  if constexpr (sizeof...(SurrogatesT) > 0) {
    return surrogateEvaluate<MaterialT, SurrogatesT...>(fileName, queryGen, materials, parameters);
  } else {
    return false;
  }
}

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
        unsigned vertLids[Cell::NumVertices]{};
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
  easi::Query query(cellToVertex_.size, Cell::Dim);

#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < cellToVertex_.size; ++elem) {
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
  for (std::size_t elem = 0; elem < cellToVertex_.size; ++elem) {
    for (std::size_t i = 0; i < NumQuadpoints; ++i) {
      const auto vertices = cellToVertex_.elementCoordinates(elem);

      const Eigen::Vector3d transformed = seissol::transformations::tetrahedronReferenceToGlobal(
          vertices[0], vertices[1], vertices[2], vertices[3], quadraturePoints_[i].data());

      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        query.x(elem * NumQuadpoints + i, j) = transformed(j);
      }

      query.group(elem * NumQuadpoints + i) = cellToVertex_.elementGroups(elem);
    }
  }

  return query;
}

std::size_t PlasticityPointGenerator::outputPerCell() const {
  constexpr auto PlasticityPoints = model::PlasticityData::PointCount;
  return pointwise_ ? PlasticityPoints : 1;
}

easi::Query PlasticityPointGenerator::generate() const {

  const auto pointsPerCell = outputPerCell();

  // Generate query using quadrature points for each element
  easi::Query query(cellToVertex_.size * pointsPerCell, Cell::Dim);

  const auto nodes = init::vNodes::view::create(init::vNodes::Values);

// Transform quadrature points to global coordinates for all elements
#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < cellToVertex_.size; ++elem) {

    const auto vertices = cellToVertex_.elementCoordinates(elem);

    for (std::size_t i = 0; i < pointsPerCell; ++i) {

      std::array<double, Cell::Dim> point{};

      if (pointwise_) {
        for (std::size_t j = 0; j < Cell::Dim; ++j) {
          if (nodes.isInRange(i, j)) {
            point[j] = nodes(i, j);
          }
        }
      } else {
        point = {1 / 4., 1 / 4., 1 / 4.};
      }

      const auto pointIdx = elem * pointsPerCell + i;

      const Eigen::Vector3d transformed = seissol::transformations::tetrahedronReferenceToGlobal(
          vertices[0], vertices[1], vertices[2], vertices[3], point.data());

      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        query.x(pointIdx, j) = transformed(j);
      }
      query.group(pointIdx) = cellToVertex_.elementGroups(elem);
    }
  }

  return query;
}

easi::Query FaultBarycenterGenerator::generate() const {
  const std::vector<Fault>& fault = meshReader_.getFault();
  const std::vector<Element>& elements = meshReader_.getElements();
  const std::vector<Vertex>& vertices = meshReader_.getVertices();

  easi::Query query(numberOfPoints_ * fault.size(), Cell::Dim);
  std::size_t q = 0;
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
    for (std::size_t n = 0; n < numberOfPoints_; ++n, ++q) {
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
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
  const auto pointsView = init::quadpoints::view::create(init::quadpoints::Values);
  easi::Query query(NumPoints * faceIDs_.size(), Cell::Dim);
  std::size_t q = 0;
  // loop over all fault elements which are managed by this generator
  // note: we have one generator per LTS layer
  for (const auto& faultId : faceIDs_) {
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
    for (std::size_t n = 0; n < NumPoints; ++n, ++q) {
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
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
        query.x(q, dim) = xyz(dim);
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

namespace {

template <typename T>
struct MaterialAverager {
  [[maybe_unused]] static constexpr bool Implemented = false;
  static T computeAveragedMaterial(std::size_t elementIdx,
                                   const std::vector<double>& quadratureWeights,
                                   const std::vector<T>& materialsFromQuery) {
    const auto numQuadPoints = quadratureWeights.size();
    return materialsFromQuery[elementIdx * numQuadPoints];
  }
};

// Computes the averaged material, assuming that materialsFromQuery, stores
// NUM_QUADPOINTS material samples per mesh element.
// We assume that materialsFromQuery[i * NUM_QUADPOINTS, ..., (i+1)*NUM_QUADPOINTS-1]
// stores samples from element i.

template <>
struct MaterialAverager<AcousticMaterial> {
  [[maybe_unused]] static constexpr bool Implemented = true;
  static AcousticMaterial
      computeAveragedMaterial(std::size_t elementIdx,
                              const std::vector<double>& quadratureWeights,
                              const std::vector<AcousticMaterial>& materialsFromQuery) {

    // (code originally extracted from the ElasticMaterial specialization below)

    double rhoMean = 0.0;

    // Average of the bulk modulus, used for acoustic material
    double kMeanInv = 0.0;

    for (std::size_t quadPointIdx = 0; quadPointIdx < NumQuadpoints; ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = NumQuadpoints * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery[globalPointIdx];
      rhoMean += elementMaterial.rho * quadWeight;
      kMeanInv += 1.0 / elementMaterial.lambda * quadWeight;
    }

    AcousticMaterial result{};
    result.rho = rhoMean;

    // Harmonic average is used for mu/K, so take the reciprocal
    result.lambda = 1.0 / kMeanInv;

    return result;
  }
};

template <>
struct MaterialAverager<ElasticMaterial> {
  [[maybe_unused]] static constexpr bool Implemented = true;
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

    for (std::size_t quadPointIdx = 0; quadPointIdx < NumQuadpoints; ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = NumQuadpoints * elementIdx + quadPointIdx;
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

template <std::size_t Mechanisms>
struct MaterialAverager<ViscoElasticMaterialParametrized<Mechanisms>> {
  [[maybe_unused]] static constexpr bool Implemented = true;
  static ViscoElasticMaterialParametrized<Mechanisms> computeAveragedMaterial(
      std::size_t elementIdx,
      const std::vector<double>& quadratureWeights,
      const std::vector<ViscoElasticMaterialParametrized<Mechanisms>>& materialsFromQuery) {
    double muMeanInv = 0.0;
    double rhoMean = 0.0;
    double vERatioMean = 0.0;
    double qpMean = 0.0;
    double qsMean = 0.0;

    for (std::size_t quadPointIdx = 0; quadPointIdx < NumQuadpoints; ++quadPointIdx) {
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = NumQuadpoints * elementIdx + quadPointIdx;
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

    ViscoElasticMaterialParametrized<Mechanisms> result{};
    result.rho = rhoMean;
    result.mu = muMean;
    result.lambda = lambdaMean;
    result.qp = qpMean;
    result.qs = qsMean;

    return result;
  }
};
} // namespace

template <class T>
void MaterialParameterDB<T>::evaluateModel(const std::string& fileName,
                                           const QueryGenerator& queryGen) {
  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  auto suppliedParameters = model->suppliedParameters();

  // the following code does:
  // * try to evaluate the model just normally
  // * if not (e.g. poroelastic or anisotropic with just elastic parameters), parse elastic or
  // acoustic first, then convert to the target material

  const auto evaluateModel = [&]() {
    easi::Query query = queryGen.generate();
    const std::size_t numPoints = query.numPoints();

    std::vector<T> materialsFromQuery(numPoints);
    easi::ArrayOfStructsAdapter<T> adapter(materialsFromQuery.data());
    for (const auto& [name, pointer] : T::ParameterMap) {
      adapter.addBindingPoint(name, pointer);
    }

    // allocate points if not done
    if (materials_->size() < numPoints) {
      materials_->resize(numPoints);
    }

    easiEvalSafe(model, query, adapter, "volume material:" + T::Text);

    return materialsFromQuery;
  };

  if (canEvaluateFor<T>(suppliedParameters)) {
    const auto materialsFromQuery = evaluateModel();
    const std::size_t numPoints = materialsFromQuery.size();

    // Only use homogenization when ElementAverageGenerator has been supplied
    if (const auto* gen = dynamic_cast<const ElementAverageGenerator*>(&queryGen)) {
      const std::size_t numElems = numPoints / NumQuadpoints;
      const std::vector<double> quadratureWeights(gen->getQuadratureWeights().begin(),
                                                  gen->getQuadratureWeights().end());

// Compute homogenized material parameters for every element in a specialization for the
// particular material
#pragma omp parallel for schedule(static)
      for (std::size_t elementIdx = 0; elementIdx < numElems; ++elementIdx) {
        materials_->at(elementIdx) = MaterialAverager<T>::computeAveragedMaterial(
            elementIdx, quadratureWeights, materialsFromQuery);
      }
    } else {
      // Usual behavior without homogenization
      for (std::size_t i = 0; i < numPoints; ++i) {
        materials_->at(i) = T(materialsFromQuery[i]);
      }
    }
  } else {
    // hard-code tests for elastic or acoustic material here (only if a conversion constructor
    // exists)
    if (!surrogateEvaluate<T, ElasticMaterial, AcousticMaterial>(
            fileName, queryGen, materials_, suppliedParameters)) {

      // no surrogate worked
      // fail gracefully by just trying to evaluate the original model and fail there
      (void)evaluateModel();
    }
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
  for (std::size_t i = 0; i < NumNodes; ++i) {
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
  auto constantTerms = init::easiBoundaryConstant::view::create(constantTermsData);

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

std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool useCellHomogenizedMaterial,
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
