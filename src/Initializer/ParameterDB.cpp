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
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "Equations/acoustic/Model/Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/viscoacoustic/Model/Datastructures.h"
#include "Equations/viscoelastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/PUMLReader.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/Quadrature.h"
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
#include <cstdint>
#include <easi/Component.h>
#include <easi/Query.h>
#include <exception>
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef USE_HDF
// PUML.h needs to be included before Downward.h

#include <PUML/Downward.h>
#endif

#ifdef USE_ASAGI
#include "Common/Real.h"
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
      std::vector<SurrogateMaterialT> preMaterials;
      edb.setMaterialVector(&preMaterials);
      edb.evaluateModel(fileName, queryGen);

      // the surrogate database sizes its own vector; callers that pre-size the target expect one
      // material per query, so a mismatch here means query and evaluation disagree
      assert(materials->empty() || materials->size() == preMaterials.size());
      materials->resize(preMaterials.size());
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

easi::Component* loadEasiModel(const std::string& fileName) {
  logInfo() << "Loading easi file" << fileName;
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
CellToVertexArray
    CellToVertexArray::fromPUML(const seissol::geometry::PumlMesh& mesh,
                                const std::vector<seissol::geometry::VertexOrder>& vertexOrders) {
  const int* groups = reinterpret_cast<const int*>(mesh.cellData(0));
  const auto& elements = mesh.cells();
  const auto& vertices = mesh.vertices();
  assert(vertexOrders.size() == elements.size());
  return CellToVertexArray(
      elements.size(),
      [&](size_t cell) {
        std::array<Eigen::Vector3d, 4> x;
        unsigned vertLids[Cell::NumVertices]{};
        PUML::Downward::vertices(mesh, elements[cell], vertLids);
        const auto& order = vertexOrders[cell];
        for (std::size_t vtx = 0; vtx < Cell::NumVertices; ++vtx) {
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            x[vtx](d) = vertices[vertLids[order[vtx]]].coordinate()[d];
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

CellToVertexArray CellToVertexArray::subset(const CellToVertexArray& array,
                                            std::vector<std::size_t> indices) {
  const auto shared = std::make_shared<const std::vector<std::size_t>>(std::move(indices));
  return CellToVertexArray(
      shared->size(),
      [array, shared](size_t idx) { return array.elementCoordinates((*shared)[idx]); },
      [array, shared](size_t idx) { return array.elementGroups((*shared)[idx]); });
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

ElementAverageGenerator::ElementAverageGenerator(const CellToVertexArray& cellToVertex,
                                                 std::size_t convergenceOrder)
    : cellToVertex_(cellToVertex) {
  const auto [quadraturePoints, quadratureWeights] =
      seissol::quadrature::simplexRule<3>(convergenceOrder);

  quadratureWeights_.assign(std::begin(quadratureWeights), std::end(quadratureWeights));
  quadraturePoints_.resize(quadratureWeights_.size());
  for (std::size_t i = 0; i < quadraturePoints_.size(); ++i) {
    std::copy(std::begin(quadraturePoints[i]),
              std::end(quadraturePoints[i]),
              std::begin(quadraturePoints_[i]));
  }
}

easi::Query ElementAverageGenerator::generate() const {
  const auto numQuadpoints = quadraturePoints_.size();

  // Generate query using quadrature points for each element
  easi::Query query(cellToVertex_.size * numQuadpoints, Cell::Dim);

// Transform quadrature points to global coordinates for all elements
#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < cellToVertex_.size; ++elem) {
    auto vertices = cellToVertex_.elementCoordinates(elem);
    const auto transform = seissol::geometry::AffineTransform(vertices);
    for (std::size_t i = 0; i < numQuadpoints; ++i) {
      const auto transformed = transform.refToSpace(quadraturePoints_[i]);
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        query.x(elem * numQuadpoints + i, d) = transformed[d];
      }
      query.group(elem * numQuadpoints + i) = cellToVertex_.elementGroups(elem);
    }
  }

  return query;
}

std::size_t PlasticityPointGenerator::outputPerCell() const {
  return pointwise_ ? nodes_.size() : 1;
}

easi::Query PlasticityPointGenerator::generate() const {

  const auto pointsPerCell = outputPerCell();

  // Generate query using quadrature points for each element
  easi::Query query(cellToVertex_.size * pointsPerCell, Cell::Dim);

// Transform quadrature points to global coordinates for all elements
#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < cellToVertex_.size; ++elem) {

    const auto vertices = cellToVertex_.elementCoordinates(elem);
    const auto transform = seissol::geometry::AffineTransform(vertices);

    for (std::size_t i = 0; i < pointsPerCell; ++i) {

      std::array<double, Cell::Dim> point{};

      if (pointwise_) {
        point = nodes_[i];
      } else {
        point = {1 / 4., 1 / 4., 1 / 4.};
      }

      const auto pointIdx = elem * pointsPerCell + i;

      const auto transformed = transform.refToSpace(point);

      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        query.x(pointIdx, d) = transformed[d];
      }
      query.group(pointIdx) = cellToVertex_.elementGroups(elem);
    }
  }

  return query;
}

easi::Query FaultBarycenterGenerator::generate() const {
  const std::vector<Fault>& fault = meshReader_.getFault();
  const std::vector<Element>& elements = meshReader_.getElements();

  easi::Query query(numberOfPoints_ * fault.size(), Cell::Dim);
  std::size_t q = 0;
  for (const Fault& f : fault) {
    std::size_t element = 0;
    std::int8_t side = 0;
    if (f.element.hasValue()) {

      element = f.element.value();
      side = f.side;
    } else {

      assert(f.neighborElement.hasValue());

      element = f.neighborElement.value();
      side = f.neighborSide;
    }

    const auto barycenter =
        seissol::geometry::AffineFaceTransform::fromMeshCell(element, side, meshReader_).center();
    for (std::size_t n = 0; n < numberOfPoints_; ++n, ++q) {
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
        query.x(q, dim) = barycenter(dim);
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

template <typename Cfg>
easi::Query FaultGPGenerator<Cfg>::generate() const {
  const std::vector<Fault>& fault = meshReader_.getFault();
  const std::vector<Element>& elements = meshReader_.getElements();
  auto cellToVertex = CellToVertexArray::fromMeshReader(meshReader_);

  constexpr size_t NumPoints = dr::misc::NumPaddedPointsSingleSim<Cfg>;
  const auto pointsView = init::quadpoints<Cfg>::view::create(init::quadpoints<Cfg>::Values);
  easi::Query query(NumPoints * faceIDs_.size(), Cell::Dim);
  std::size_t q = 0;
  // loop over all fault elements which are managed by this generator
  // note: we have one generator per LTS layer
  for (const auto& faultId : faceIDs_) {
    const Fault& f = fault.at(faultId);
    std::size_t element = 0;
    std::int8_t side = 0;
    auto sideOrientation = seissol::geometry::FaceOrientation::Local;
    if (f.element.hasValue()) {
      element = f.element.value();
      side = f.side;
    } else {
      assert(f.neighborElement.hasValue());

      element = f.neighborElement.value();
      side = f.neighborSide;
      // the canonical vertex numbering pins the face orientation index to zero
      sideOrientation = seissol::geometry::FaceOrientation::Rotate0;
    }

    auto coords = cellToVertex.elementCoordinates(element);
    const auto face = seissol::geometry::AffineFaceTransform(
        seissol::geometry::AffineTransform(coords),
        seissol::geometry::ReferenceFaceMap(side, sideOrientation));
    for (std::size_t n = 0; n < NumPoints; ++n, ++q) {
      auto localPoints = seissol::geometry::FaceTransform::FaceVectorT(
          seissol::multisim::multisimTranspose<Cfg>(pointsView, n, 0),
          seissol::multisim::multisimTranspose<Cfg>(pointsView, n, 1));
      // padded points are in the middle of the tetrahedron
      if (n >= dr::misc::NumBoundaryGaussPoints<Cfg>) {
        localPoints =
            seissol::geometry::FaceTransform::FaceVectorT(Face::ReferenceBarycenter.data());
      }

      const auto xyz = face.refToSpace(localPoints);
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
        query.x(q, dim) = xyz(dim);
      }
      query.group(q) = elements[element].faultTags[side];
    }
  }
  return query;
}

#define SEISSOL_INSTANTIATE(Cfg) template class FaultGPGenerator<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

namespace {

template <typename T>
struct MaterialAverager {
  [[maybe_unused]] static constexpr bool Implemented = false;
  static T computeAveragedMaterial(std::size_t elementIdx,
                                   const std::vector<double>& quadratureWeights,
                                   const std::function<const T&(std::size_t)>& materialsFromQuery) {
    const auto numQuadPoints = quadratureWeights.size();
    return materialsFromQuery(elementIdx * numQuadPoints);
  }
};

// Computes the averaged material, assuming that materialsFromQuery, stores
// NUM_QUADPOINTS material samples per mesh element.
// We assume that materialsFromQuery[i * NUM_QUADPOINTS, ..., (i+1)*NUM_QUADPOINTS-1]
// stores samples from element i.

template <>
struct MaterialAverager<AcousticMaterial> {
  [[maybe_unused]] static constexpr bool Implemented = true;
  static AcousticMaterial computeAveragedMaterial(
      std::size_t elementIdx,
      const std::vector<double>& quadratureWeights,
      const std::function<const AcousticMaterial&(std::size_t)>& materialsFromQuery) {

    // (code originally extracted from the ElasticMaterial specialization below)

    double rhoMean = 0.0;

    // Average of the bulk modulus, used for acoustic material
    double kMeanInv = 0.0;

    for (std::size_t quadPointIdx = 0; quadPointIdx < quadratureWeights.size(); ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = quadratureWeights.size() * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery(globalPointIdx);
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
  static ElasticMaterial computeAveragedMaterial(
      std::size_t elementIdx,
      const std::vector<double>& quadratureWeights,
      const std::function<const ElasticMaterial&(std::size_t)>& materialsFromQuery) {
    double muMeanInv = 0.0;
    double rhoMean = 0.0;
    // Average of v / E with v: Poisson's ratio, E: Young's modulus
    double vERatioMean = 0.0;

    // Average of the bulk modulus, used for acoustic material
    double kMeanInv = 0.0;

    // Acoustic material has zero mu. This is a special case because the harmonic mean of a set
    // of numbers that includes zero is defined as zero.
    // Hence: If part of the element is acoustic, the entire element is considered to be acoustic.
    bool isAcoustic = false;

    // important: scan for acousticity _first_.
    for (std::size_t quadPointIdx = 0; quadPointIdx < quadratureWeights.size(); ++quadPointIdx) {
      const std::size_t globalPointIdx = quadratureWeights.size() * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery(globalPointIdx);
      isAcoustic |= elementMaterial.mu == 0.0;
    }

    for (std::size_t quadPointIdx = 0; quadPointIdx < quadratureWeights.size(); ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = quadratureWeights.size() * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery(globalPointIdx);
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
struct MaterialAverager<ViscoElasticMaterial<Mechanisms>> {
  [[maybe_unused]] static constexpr bool Implemented = true;
  static ViscoElasticMaterial<Mechanisms> computeAveragedMaterial(
      std::size_t elementIdx,
      const std::vector<double>& quadratureWeights,
      const std::function<const ViscoElasticMaterial<Mechanisms>&(std::size_t)>&
          materialsFromQuery) {
    double qpMean = 0.0;
    double qsMean = 0.0;

    for (std::size_t quadPointIdx = 0; quadPointIdx < quadratureWeights.size(); ++quadPointIdx) {
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = quadratureWeights.size() * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery(globalPointIdx);
      qpMean += elementMaterial.qp * quadWeight;
      qsMean += elementMaterial.qs * quadWeight;
    }

    const auto base = MaterialAverager<ElasticMaterial>::computeAveragedMaterial(
        elementIdx,
        quadratureWeights,
        [materialsFromQuery](std::size_t index) -> const ElasticMaterial& {
          return materialsFromQuery(index);
        });

    auto result = ViscoElasticMaterial<Mechanisms>(base);
    result.qp = qpMean;
    result.qs = qsMean;

    return result;
  }
};

template <std::size_t Mechanisms>
struct MaterialAverager<ViscoAcousticMaterial<Mechanisms>> {
  [[maybe_unused]] static constexpr bool Implemented = true;
  static ViscoAcousticMaterial<Mechanisms> computeAveragedMaterial(
      std::size_t elementIdx,
      const std::vector<double>& quadratureWeights,
      const std::function<const ViscoAcousticMaterial<Mechanisms>&(std::size_t)>&
          materialsFromQuery) {
    double qpMean = 0.0;

    for (std::size_t quadPointIdx = 0; quadPointIdx < quadratureWeights.size(); ++quadPointIdx) {
      const double quadWeight = 6.0 * quadratureWeights[quadPointIdx];
      const std::size_t globalPointIdx = quadratureWeights.size() * elementIdx + quadPointIdx;
      const auto& elementMaterial = materialsFromQuery(globalPointIdx);
      qpMean += elementMaterial.qp * quadWeight;
    }

    const auto base = MaterialAverager<AcousticMaterial>::computeAveragedMaterial(
        elementIdx,
        quadratureWeights,
        [materialsFromQuery](std::size_t index) -> const AcousticMaterial& {
          return materialsFromQuery(index);
        });

    auto result = ViscoAcousticMaterial<Mechanisms>(base);
    result.qp = qpMean;

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

    easiEvalSafe(model, query, adapter, "volume material:" + T::Text);

    return materialsFromQuery;
  };

  if (canEvaluateFor<T>(suppliedParameters)) {
    const auto materialsFromQuery = evaluateModel();
    const std::size_t numPoints = materialsFromQuery.size();

    // Only use homogenization when ElementAverageGenerator has been supplied
    if (const auto* gen = dynamic_cast<const ElementAverageGenerator*>(&queryGen)) {
      const auto& quadratureWeights = gen->getQuadratureWeights();
      const std::size_t numElems = numPoints / quadratureWeights.size();

      // allocate output array
      materials_->resize(numElems);

      // Compute homogenized material parameters for every element in a specialization for the
      // particular material

#pragma omp parallel for schedule(static)
      for (std::size_t elementIdx = 0; elementIdx < numElems; ++elementIdx) {
        materials_->at(elementIdx) = MaterialAverager<T>::computeAveragedMaterial(
            elementIdx, quadratureWeights, [&materialsFromQuery](std::size_t index) -> const T& {
              return materialsFromQuery[index];
            });
      }
    } else {
      // allocate output array
      materials_->resize(numPoints);

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

template <typename T>
void FaultParameterDB<T>::evaluateModel(const std::string& fileName,
                                        const QueryGenerator& queryGen) {
  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  easi::Query query = queryGen.generate();

  easi::ArraysAdapter<T> adapter;
  for (auto& kv : parameters_) {
    adapter.addBindingPoint(kv.first, kv.second.first + simid_, kv.second.second * numSimulations_);
  }

  easiEvalSafe(model, query, adapter, "fault material");

  delete model;
}

template class FaultParameterDB<float>;
template class FaultParameterDB<double>;

std::set<std::string> faultProvides(const std::string& fileName) {
  if (fileName.empty()) {
    return {};
  }

  // NOLINTNEXTLINE(misc-const-correctness)
  easi::Component* const model = loadEasiModel(fileName);
  std::set<std::string> supplied = model->suppliedParameters();
  delete model;
  return supplied;
}

DirichletCondition::DirichletCondition(const std::string& fileName)
    : model_(loadEasiModel(fileName)) {}

// the model is owned by one condition at a time, and the one moved from no longer deletes it
DirichletCondition::DirichletCondition(DirichletCondition&& other) noexcept
    : model_(std::exchange(other.model_, nullptr)) {}

DirichletCondition& DirichletCondition::operator=(DirichletCondition&& other) noexcept {
  std::swap(model_, other.model_);
  return *this;
}

DirichletCondition::~DirichletCondition() { delete model_; }

template <typename Cfg>
BoundaryFrame DirichletCondition::query(const double* barycenter,
                                        Real<Cfg>* mapTermsData,
                                        Real<Cfg>* constantTermsData) const {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  if (model_ == nullptr) {
    logError() << "Model for easi-provided boundary is not initialized.";
  }
  assert(mapTermsData != nullptr);
  assert(constantTermsData != nullptr);

  // The boundary condition is constant over the face, so it is sampled at the
  // face barycenter.
  auto query = easi::Query{1, 3};
  query.x(0, 0) = barycenter[0];
  query.x(0, 1) = barycenter[1];
  query.x(0, 2) = barycenter[2];
  query.group(0) = 1;

  const auto& supplied = model_->suppliedParameters();

  // The ghost cell state is an affine function of the interior state, given in
  // global coordinates: q_ghost = A q_inside + b. The entries of A are named
  // map_{to}_{from}, those of b const_{to}, where the quantity names are the
  // ones of the material at hand. Mirroring the x velocity at the ghost cell is
  // therefore map_v1_v1: -1.
  const auto& varNames = model::MaterialOf<Cfg>::Quantities;

  auto mapTerms = init::dirichletMapGlobal<Cfg>::view::create(mapTermsData);
  auto constantTerms = init::dirichletOffsetGlobal<Cfg>::view::create(constantTermsData);

  easi::ArraysAdapter<real> adapter{};
  std::unordered_set<std::string> known;

  // easi supplies numbers, so the frame is stated as one: 0 for global, 1 for face-aligned.
  real frame = 0.0;
  known.insert("frame");
  if (supplied.count("frame") > 0) {
    adapter.addBindingPoint("frame", &frame);
  }

  for (size_t i = 0; i < varNames.size(); ++i) {
    const auto termName = std::string{"const_"} + varNames[i];
    known.insert(termName);
    auto& term = multisim::multisimWrap<Cfg>(constantTerms, 0, i);
    if (supplied.count(termName) > 0) {
      adapter.addBindingPoint(termName, &term);
    } else {
      term = 0.0;
    }
  }
  for (size_t i = 0; i < varNames.size(); ++i) {
    for (size_t j = 0; j < varNames.size(); ++j) {
      auto termName = std::string{"map_"};
      termName += varNames[i];
      termName += "_";
      termName += varNames[j];
      known.insert(termName);
      if (supplied.count(termName) > 0) {
        adapter.addBindingPoint(termName, &mapTerms(i, j));
      } else {
        // Default: Extrapolate
        mapTerms(i, j) = (i == j) ? 1.0 : 0.0;
      }
    }
  }

  for (const auto& termName : supplied) {
    if (known.count(termName) == 0) {
      std::ostringstream valid;
      for (size_t i = 0; i < varNames.size(); ++i) {
        valid << (i == 0 ? "" : ", ") << varNames[i];
      }
      logError() << "The boundary condition file supplies" << termName
                 << "which is not a term of the boundary condition. Terms are named"
                 << "map_{to}_{from} and const_{to}, where both quantity names are one of:"
                 << valid.str() << ".";
    }
  }

  easiEvalSafe(model_, query, adapter, "Dirichlet BC data");

  if (frame != 0.0 && frame != 1.0) {
    logError() << "The boundary condition file supplies a frame of" << frame
               << "-- it has to be 0 for a condition stated in global coordinates, or 1 for one "
                  "stated in the face-aligned basis.";
  }

  // The condition does not depend on the simulation index, so every fused
  // simulation gets the same one.
  for (std::size_t sim = 1; sim < Cfg::NumSimulations; ++sim) {
    for (size_t i = 0; i < varNames.size(); ++i) {
      multisim::multisimWrap<Cfg>(constantTerms, sim, i) =
          multisim::multisimWrap<Cfg>(constantTerms, 0, i);
    }
  }

  return frame == 0.0 ? BoundaryFrame::Global : BoundaryFrame::FaceAligned;
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template BoundaryFrame DirichletCondition::query<Cfg>(const double*, Real<Cfg>*, Real<Cfg>*)     \
      const;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

template <typename MaterialT>
std::shared_ptr<QueryGenerator> getBestQueryGenerator(bool useCellHomogenizedMaterial,
                                                      const CellToVertexArray& cellToVertex,
                                                      std::size_t convergenceOrder) {
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
      queryGen = std::make_shared<ElementAverageGenerator>(cellToVertex, convergenceOrder);
    }
  }
  return queryGen;
}

// the argument is a type, which the check takes for an expression in template arguments
// NOLINTBEGIN(bugprone-macro-parentheses)
#define SEISSOL_MATERIAL_INSTANTIATE(Cfg)                                                          \
  template std::shared_ptr<QueryGenerator> getBestQueryGenerator<seissol::model::MaterialOf<Cfg>>( \
      bool, const CellToVertexArray&, std::size_t);                                                \
  template class MaterialParameterDB<seissol::model::MaterialOf<Cfg>>;
// NOLINTEND(bugprone-macro-parentheses)
SEISSOL_FOR_EACH_MATERIAL(SEISSOL_MATERIAL_INSTANTIATE)
#undef SEISSOL_MATERIAL_INSTANTIATE
template class MaterialParameterDB<seissol::model::Plasticity>;

} // namespace seissol::initializer
