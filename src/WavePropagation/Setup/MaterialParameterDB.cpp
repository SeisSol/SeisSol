// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "MaterialParameterDB.h"

#include "Config.h"
#include "Core/Common/Constants.h"
#include "Core/Data/DataReader.h"
#include "Core/Data/DataTable.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/QueryGenerator.h"
#include "Core/Numerical/Quadrature.h"
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
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <utils/logger.h>
#include <vector>

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

} // namespace

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
void MaterialParameterDB<T>::evaluate(data::DataReader& model, const QueryGenerator& queryGen) {
  const auto supplied = data::suppliedParameters(model);

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

    data::evaluate(model, table, "volume material:" + T::Text);

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
