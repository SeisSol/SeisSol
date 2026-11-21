// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "ParameterDB.h"

#include "Config.h"
#include "Core/Common/Constants.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Data/DataReader.h"
#include "Core/Data/DataTable.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/FaceTransform.h"
#include "Core/Geometry/MeshDefinition.h"
#include "Core/Geometry/PUMLReader.h"
#include "Core/Numerical/Quadrature.h"
#include "Features/DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"
#include "Simulation/Reader/ReaderBuilder.h"
#include "Simulation/SeisSol.h"
#include "WavePropagation/Physics/CommonDatastructures.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Physics/acoustic/Datastructures.h"
#include "WavePropagation/Physics/elastic/Datastructures.h"
#include "WavePropagation/Physics/viscoacoustic/Datastructures.h"
#include "WavePropagation/Physics/viscoelastic/Datastructures.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
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
#include <PUML/Types.h>
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
bool surrogateEvaluate(data::DataReader& model,
                       const QueryGenerator& queryGen,
                       std::vector<MaterialT>* materials,
                       const std::set<std::string>& parameters) {
  if constexpr (std::is_constructible_v<MaterialT, SurrogateMaterialT>) {
    if (canEvaluateFor<SurrogateMaterialT>(parameters)) {
      MaterialParameterDB<SurrogateMaterialT> edb;
      std::vector<SurrogateMaterialT> preMaterials;
      edb.setMaterialVector(&preMaterials);
      edb.evaluate(model, queryGen);

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
    return surrogateEvaluate<MaterialT, SurrogatesT...>(model, queryGen, materials, parameters);
  } else {
    return false;
  }
}

void evaluateSafe(data::DataReader& model, const data::DataTable& table, const std::string& hint) {
  try {
    model.call(table);
  } catch (const std::exception& error) {
    logError() << "Error while evaluating a model for" << hint.c_str() << ":"
               << std::string(error.what());
  }
}

std::set<std::string> suppliedParameters(data::DataReader& model) {
  const auto& outputs = model.outputVars();
  return {outputs.begin(), outputs.end()};
}

/// Binds the point set of a query, computed a range of points at a time: `coordinates(first,
/// count, xyz)` writes the coordinates of the points [first, first + count) to xyz[i * 3 + d], and
/// `groups(first, count, out)` their groups. Nothing is materialised; the table refers to whatever
/// the two callbacks capture.
template <typename CoordinateFn, typename GroupFn>
void bindPointSet(data::DataTable& table, CoordinateFn coordinates, GroupFn groups) {
  const auto shared = std::make_shared<CoordinateFn>(std::move(coordinates));
  for (std::size_t d = 0; d < Cell::Dim; ++d) {
    table.bindComputedBatch<double>(std::string(1, "xyz"[d]),
                                    [shared, d](std::size_t first, std::size_t count, double* out) {
                                      thread_local std::vector<double> xyz;
                                      xyz.resize(count * Cell::Dim);
                                      (*shared)(first, count, xyz.data());
                                      for (std::size_t i = 0; i < count; ++i) {
                                        out[i] = xyz[i * Cell::Dim + d];
                                      }
                                    });
  }
  table.bindComputedBatch<std::int32_t>("group", std::move(groups));
}

/// Points given per cell in reference coordinates, `pointsPerCell` of them for every cell of
/// `cells`: the transformation of a cell is set up once per range and cell rather than per point.
void cellPointCoordinates(const CellToVertexArray& cells,
                          const std::vector<std::array<double, Cell::Dim>>& referencePoints,
                          std::size_t first,
                          std::size_t count,
                          double* xyz) {
  const auto pointsPerCell = referencePoints.size();
  std::optional<seissol::geometry::AffineTransform> transform;
  std::size_t transformCell = 0;
  for (std::size_t i = 0; i < count; ++i) {
    const auto point = first + i;
    const auto cell = point / pointsPerCell;
    if (!transform.has_value() || transformCell != cell) {
      transform.emplace(cells.elementCoordinates(cell));
      transformCell = cell;
    }
    const auto transformed = transform->refToSpace(referencePoints[point % pointsPerCell]);
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      xyz[i * Cell::Dim + d] = transformed[d];
    }
  }
}

} // namespace

std::unique_ptr<data::DataReader> ParameterDB::loadModel(const std::string& fileName) {
  return data::buildReader(fileName, {"x", "y", "z"});
}

CellToVertexArray::CellToVertexArray(size_t size,
                                     CellToVertexFunction elementCoordinates,
                                     CellToGroupFunction elementGroups)
    : size(size), elementCoordinates(std::move(elementCoordinates)),
      elementGroups(std::move(elementGroups)) {}

CellToVertexArray
    CellToVertexArray::fromMeshReader(const seissol::geometry::MeshReader& meshReader) {
  const auto& elements = meshReader.getElements();
  const auto& vertices = meshReader.getVertices();

  return CellToVertexArray(
      elements.size(),
      [&](size_t index) {
        std::array<Eigen::Vector3d, 4> verts;
        for (size_t i = 0; i < Cell::NumVertices; ++i) {
          const auto vindex = elements[index].vertices[i];
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
                                const seissol::geometry::PumlCellArrays& arrays,
                                const std::vector<seissol::geometry::VertexOrder>& vertexOrders) {
  const auto groups = mesh.data(arrays.group);
  const auto& elements = mesh.cells();
  const auto& vertices = mesh.vertices();
  assert(vertexOrders.size() == elements.size());
  return CellToVertexArray(
      elements.size(),
      [&](size_t cell) {
        std::array<Eigen::Vector3d, 4> x;
        PUML::LocalId vertLids[Cell::NumVertices]{};
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

data::DataTable ElementBarycenterGenerator::generate() const {
  data::DataTable table(cellToVertex_.size);
  bindPointSet(
      table,
      [this](std::size_t first, std::size_t count, double* xyz) {
        for (std::size_t i = 0; i < count; ++i) {
          const auto vertices = cellToVertex_.elementCoordinates(first + i);
          const Eigen::Vector3d barycenter =
              (vertices[0] + vertices[1] + vertices[2] + vertices[3]) * 0.25;
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            xyz[i * Cell::Dim + d] = barycenter(d);
          }
        }
      },
      [this](std::size_t first, std::size_t count, std::int32_t* out) {
        for (std::size_t i = 0; i < count; ++i) {
          out[i] = cellToVertex_.elementGroups(first + i);
        }
      });
  return table;
}

ElementAverageGenerator::ElementAverageGenerator(CellToVertexArray cellToVertex,
                                                 std::size_t convergenceOrder)
    : cellToVertex_(std::move(cellToVertex)) {
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

data::DataTable ElementAverageGenerator::generate() const {
  const auto numQuadpoints = quadraturePoints_.size();

  // the quadrature points of every element
  data::DataTable table(cellToVertex_.size * numQuadpoints);
  bindPointSet(
      table,
      [this](std::size_t first, std::size_t count, double* xyz) {
        cellPointCoordinates(cellToVertex_, quadraturePoints_, first, count, xyz);
      },
      [this, numQuadpoints](std::size_t first, std::size_t count, std::int32_t* out) {
        for (std::size_t i = 0; i < count; ++i) {
          out[i] = cellToVertex_.elementGroups((first + i) / numQuadpoints);
        }
      });
  return table;
}

std::size_t PlasticityPointGenerator::outputPerCell() const {
  return pointwise_ ? nodes_.size() : 1;
}

data::DataTable PlasticityPointGenerator::generate() const {
  const auto pointsPerCell = outputPerCell();

  // the plasticity nodes of every element, or its barycenter only
  auto referencePoints =
      pointwise_ ? nodes_ : std::vector<std::array<double, Cell::Dim>>{{1 / 4., 1 / 4., 1 / 4.}};

  data::DataTable table(cellToVertex_.size * pointsPerCell);
  bindPointSet(
      table,
      [this, referencePoints = std::move(referencePoints)](
          std::size_t first, std::size_t count, double* xyz) {
        cellPointCoordinates(cellToVertex_, referencePoints, first, count, xyz);
      },
      [this, pointsPerCell](std::size_t first, std::size_t count, std::int32_t* out) {
        for (std::size_t i = 0; i < count; ++i) {
          out[i] = cellToVertex_.elementGroups((first + i) / pointsPerCell);
        }
      });
  return table;
}

template <typename Cfg>
data::DataTable FaultGPGenerator<Cfg>::generate() const {
  constexpr size_t NumPoints = dr::misc::NumPaddedPointsSingleSim<Cfg>;

  // element, side and the face transform of each fault face managed by this generator (we have
  // one generator per LTS layer), set up once rather than per point and column
  struct FaultFace {
    std::size_t element{0};
    std::int8_t side{0};
    seissol::geometry::AffineFaceTransform transform;
  };
  const std::vector<Fault>& fault = meshReader_.getFault();
  const auto cellToVertex = CellToVertexArray::fromMeshReader(meshReader_);
  auto faces = std::make_shared<std::vector<FaultFace>>();
  faces->reserve(faceIDs_.size());
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
    faces->push_back(
        FaultFace{.element = element,
                  .side = side,
                  .transform = seissol::geometry::AffineFaceTransform(
                      seissol::geometry::AffineTransform(cellToVertex.elementCoordinates(element)),
                      seissol::geometry::ReferenceFaceMap(side, sideOrientation))});
  }

  data::DataTable table(NumPoints * faceIDs_.size());
  bindPointSet(
      table,
      [faces](std::size_t first, std::size_t count, double* xyz) {
        const auto pointsView = init::quadpoints<Cfg>::view::create(init::quadpoints<Cfg>::Values);
        for (std::size_t i = 0; i < count; ++i) {
          const auto q = first + i;
          const auto n = q % NumPoints;
          auto localPoints = seissol::geometry::FaceTransform::FaceVectorT(
              seissol::multisim::multisimTranspose<Cfg>(pointsView, n, 0),
              seissol::multisim::multisimTranspose<Cfg>(pointsView, n, 1));
          // padded points are in the middle of the tetrahedron
          if (n >= dr::misc::NumBoundaryGaussPoints<Cfg>) {
            localPoints =
                seissol::geometry::FaceTransform::FaceVectorT(Face::ReferenceBarycenter.data());
          }
          const auto point = (*faces)[q / NumPoints].transform.refToSpace(localPoints);
          for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
            xyz[i * Cell::Dim + dim] = point(dim);
          }
        }
      },
      [faces, &elements = meshReader_.getElements()](
          std::size_t first, std::size_t count, std::int32_t* out) {
        for (std::size_t i = 0; i < count; ++i) {
          const auto& face = (*faces)[(first + i) / NumPoints];
          out[i] = elements[face.element].faultTags[face.side];
        }
      });
  return table;
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

template <typename T>
std::vector<std::string> missingMaterialParameters(const std::set<std::string>& supplied) {
  // as MaterialParameterDB::evaluateModel takes the model: for the parameters of T, or else for
  // those of an elastic or acoustic material that T can be made of
  if (canEvaluateFor<T>(supplied)) {
    return {};
  }
  if constexpr (std::is_constructible_v<T, ElasticMaterial>) {
    if (canEvaluateFor<ElasticMaterial>(supplied)) {
      return {};
    }
  }
  if constexpr (std::is_constructible_v<T, AcousticMaterial>) {
    if (canEvaluateFor<AcousticMaterial>(supplied)) {
      return {};
    }
  }
  std::vector<std::string> missing;
  for (const auto& [name, pointer] : T::ParameterMap) {
    if (!supplied.contains(name)) {
      missing.push_back(name);
    }
  }
  std::ranges::sort(missing);
  return missing;
}

template <class T>
void MaterialParameterDB<T>::evaluateModel(const std::string& fileName,
                                           const QueryGenerator& queryGen) {
  const auto model = ParameterDB::loadModel(fileName);
  evaluate(*model, queryGen);
}

template <class T>
void MaterialParameterDB<T>::evaluate(data::DataReader& model, const QueryGenerator& queryGen) {
  const auto supplied = suppliedParameters(model);

  // the following code does:
  // * try to evaluate the model just normally
  // * if not (e.g. poroelastic or anisotropic with just elastic parameters), parse elastic or
  // acoustic first, then convert to the target material

  const auto evaluateModel = [&] {
    auto table = queryGen.generate();
    const std::size_t numPoints = table.numPoints();

    std::vector<T> materialsFromQuery(numPoints);
    for (const auto& [name, pointer] : T::ParameterMap) {
      table.bindMemberView(name, data::Direction::Out, materialsFromQuery.data(), pointer);
    }
    // a parameter the material file may leave out keeps its default where it does
    if constexpr (requires { T::OptionalParameterMap; }) {
      for (const auto& [name, pointer] : T::OptionalParameterMap) {
        if (supplied.contains(name)) {
          table.bindMemberView(name, data::Direction::Out, materialsFromQuery.data(), pointer);
        }
      }
    }

    evaluateSafe(model, table, "volume material:" + T::Text);

    return materialsFromQuery;
  };

  if (canEvaluateFor<T>(supplied)) {
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
            model, queryGen, materials_, supplied)) {

      // no surrogate worked
      // fail gracefully by just trying to evaluate the original model and fail there
      (void)evaluateModel();
    }
  }
}

template <typename T>
void FaultParameterDB<T>::evaluateModel(const std::string& fileName,
                                        const QueryGenerator& queryGen) {
  const auto model = ParameterDB::loadModel(fileName);
  auto table = queryGen.generate();

  // the point q of the query lies on the face faces[q / pointsPerFace] of the arrays, if the faces
  // are given; i.e. the arrays may hold the points of faces which the query does not contain, in
  // between the ones it does
  std::vector<std::uint32_t> faceIndex;
  if (faces_.has_value()) {
    faceIndex.assign(faces_->begin(), faces_->end());
  }
  for (const auto& [name, binding] : parameters_) {
    const auto [memory, stride] = binding;
    const auto pointStride = static_cast<std::size_t>(stride) * numSimulations_;
    if (faces_.has_value()) {
      table.bindCellPointView(name,
                              data::Direction::Out,
                              memory,
                              pointsPerFace_,
                              pointsPerFace_ * pointStride,
                              pointStride,
                              simid_,
                              faceIndex.data());
    } else {
      table.bindView(name, data::Direction::Out, memory, pointStride, simid_);
    }
  }

  evaluateSafe(*model, table, "fault material");
}

template class FaultParameterDB<float>;
template class FaultParameterDB<double>;

std::set<std::string> faultProvides(const std::string& fileName) {
  if (fileName.empty()) {
    return {};
  }

  const auto model = ParameterDB::loadModel(fileName);
  return suppliedParameters(*model);
}

DirichletCondition::DirichletCondition(const std::string& fileName)
    : model_(ParameterDB::loadModel(fileName)) {}

template <typename Cfg>
BoundaryFrame DirichletCondition::query(const double* barycenter,
                                        Real<Cfg>* mapTermsData,
                                        Real<Cfg>* constantTermsData,
                                        const std::vector<std::size_t>& simulations) const {
  using real = Real<Cfg>;
  if (model_ == nullptr) {
    logError() << "Model for easi-provided boundary is not initialized.";
  }
  assert(mapTermsData != nullptr);
  assert(constantTermsData != nullptr);
  assert(!simulations.empty());

  // The boundary condition is constant over the face, so it is sampled at the
  // face barycenter.
  data::DataTable table(1);
  table.bindViewConst("x", data::Direction::In, barycenter, 3, 0);
  table.bindViewConst("y", data::Direction::In, barycenter, 3, 1);
  table.bindViewConst("z", data::Direction::In, barycenter, 3, 2);
  table.bindConstant("group", std::int32_t{1});

  const auto supplied = suppliedParameters(*model_);

  // The ghost cell state is an affine function of the interior state, given in
  // global coordinates: q_ghost = A q_inside + b. The entries of A are named
  // map_{to}_{from}, those of b const_{to}, where the quantity names are the
  // ones of the material at hand. Mirroring the x velocity at the ghost cell is
  // therefore map_v1_v1: -1.
  const auto& varNames = model::MaterialOf<Cfg>::Quantities;

  auto mapTerms = init::dirichletMapGlobal<Cfg>::view::create(mapTermsData);
  auto constantTerms = init::dirichletOffsetGlobal<Cfg>::view::create(constantTermsData);

  std::unordered_set<std::string> known;

  // a model supplies numbers, so the frame is stated as one: 0 for global, 1 for face-aligned.
  real frame = 0.0;
  known.insert("frame");
  if (supplied.contains("frame")) {
    table.bindView("frame", data::Direction::Out, &frame);
  }

  for (size_t i = 0; i < varNames.size(); ++i) {
    const auto termName = std::string{"const_"} + varNames[i];
    known.insert(termName);
    auto& term = multisim::multisimWrap<Cfg>(constantTerms, simulations.front(), i);
    if (supplied.count(termName) > 0) {
      table.bindView(termName, data::Direction::Out, &term);
    } else {
      term = 0.0;
    }
  }
  for (size_t i = 0; i < varNames.size(); ++i) {
    for (size_t j = 0; j < varNames.size(); ++j) {
      auto termName = std::string{"map_"};
      termName += varNames[i];
      termName += '_';
      termName += varNames[j];
      known.insert(termName);
      if (supplied.contains(termName)) {
        table.bindView(termName, data::Direction::Out, &mapTerms(i, j));
      } else {
        // Default: Extrapolate
        mapTerms(i, j) = (i == j) ? 1.0 : 0.0;
      }
    }
  }

  for (const auto& termName : supplied) {
    if (!known.contains(termName)) {
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

  evaluateSafe(*model_, table, "Dirichlet BC data");

  if (frame != 0.0 && frame != 1.0) {
    logError() << "The boundary condition file supplies a frame of" << frame
               << "-- it has to be 0 for a condition stated in global coordinates, or 1 for one "
                  "stated in the face-aligned basis.";
  }

  // every fused simulation that takes the condition gets its constant terms
  for (const auto simulation : simulations) {
    for (size_t i = 0; i < varNames.size(); ++i) {
      multisim::multisimWrap<Cfg>(constantTerms, simulation, i) =
          multisim::multisimWrap<Cfg>(constantTerms, simulations.front(), i);
    }
  }

  return frame == 0.0 ? BoundaryFrame::Global : BoundaryFrame::FaceAligned;
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template BoundaryFrame DirichletCondition::query<Cfg>(                                           \
      const double*, Real<Cfg>*, Real<Cfg>*, const std::vector<std::size_t>&) const;
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
  template class MaterialParameterDB<seissol::model::MaterialOf<Cfg>>;                             \
  template std::vector<std::string> missingMaterialParameters<seissol::model::MaterialOf<Cfg>>(    \
      const std::set<std::string>&);
// NOLINTEND(bugprone-macro-parentheses)
SEISSOL_FOR_EACH_MATERIAL(SEISSOL_MATERIAL_INSTANTIATE)
#undef SEISSOL_MATERIAL_INSTANTIATE
template class MaterialParameterDB<seissol::model::Plasticity>;
template std::vector<std::string>
    missingMaterialParameters<seissol::model::Plasticity>(const std::set<std::string>&);

} // namespace seissol::initializer
