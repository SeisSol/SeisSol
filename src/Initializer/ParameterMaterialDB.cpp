#include "ParameterMaterialDB.hpp"
#include "Common/cellconfigconv.hpp"
#include "Equations/viscoelastic/Model/datastructures.hpp"
#include "Initializer/ParameterDB.h"
#include "Model/common_datastructures.hpp"
#include "Model/plasticity.hpp"
#include <Equations/datastructures.hpp>
#include <unordered_map>
#include <vector>
#include <array>
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "easi/Query.h"
#include "easi/Component.h"
#include "SeisSol.h"
#include "Physics/Attenuation.hpp"

namespace {
struct ModelQuadratureRule {
  ::std::vector<::std::array<double, 3>> points;
  ::std::vector<double> weights;

  ModelQuadratureRule(int order) {
    if (order == 0) {
      // explicitly barycentric (regardless of quadrature rule)
      points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
      weights = {1};
    } else {
      // TODO(David): convert to quadrature methods to array etc.
      auto length = ::seissol::initializers::numQuadPoints(order);
      auto prePoints = ::std::vector<double[3]>(length);
      weights.resize(length);
      ::seissol::quadrature::TetrahedronQuadrature(prePoints.data(), weights.data(), order);
      points.resize(length);
      for (::std::size_t i = 0; i < length; ++i) {
        points[i] = {prePoints[i][0], prePoints[i][1], prePoints[i][2]};
      }
    }
  }

  ::std::vector<::std::array<double, 3>>
      map(const ::std::array<Eigen::Vector3d, 4>& vertices) const {
    ::std::vector<::std::array<double, 3>> mappedPoints(points.size());
    ::std::array<::std::array<double, 3>, 4> verticesConv{
        ::std::array<double, 3>{vertices[0](0), vertices[0](1), vertices[0](2)},
        ::std::array<double, 3>{vertices[1](0), vertices[1](1), vertices[1](2)},
        ::std::array<double, 3>{vertices[2](0), vertices[2](1), vertices[2](2)},
        ::std::array<double, 3>{vertices[3](0), vertices[3](1), vertices[3](2)},
    };
    for (size_t i = 0; i < points.size(); ++i) {
      ::seissol::transformations::tetrahedronReferenceToGlobal(verticesConv[0].data(),
                                                               verticesConv[1].data(),
                                                               verticesConv[2].data(),
                                                               verticesConv[3].data(),
                                                               points[i].data(),
                                                               mappedPoints[i].data());
    }
  }
};

class ModelQuadratureRuleCollection {
  public:
  ModelQuadratureRuleCollection(int maxOrder) {
    for (int i = 0; i < rules.size(); ++i) {
      rules.emplace_back(ModelQuadratureRule(i));
    }
  }

  const ModelQuadratureRule& getRule(int order) const { return rules[order]; }

  private:
  ::std::vector<ModelQuadratureRule> rules;
};

// TODO(David): maybe merge with the Material class
template <typename MaterialT>
struct MaterialAverager {
  constexpr static bool Supported = false;
  static MaterialT compute(const ::std::vector<MaterialT>& materials,
                           size_t start,
                           const ModelQuadratureRule& rule) {
    if (rule.weights.size() > 1) {
      logError() << "Material averaging is not supported for" << MaterialT::Text;
    }

    return materials[start];
  }
};

template <>
struct MaterialAverager<::seissol::model::SolidMaterial> {
  constexpr static bool Supported = true;
  using MaterialT = ::seissol::model::SolidMaterial;
  static MaterialT compute(const ::std::vector<MaterialT>& materials,
                           size_t start,
                           const ModelQuadratureRule& rule) {
    double rhoMean = 0.0;

    for (unsigned quadPointIdx = 0; quadPointIdx < rule.weights.size(); ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * rule.weights[quadPointIdx];
      const size_t globalPointIdx = start + quadPointIdx;
      const auto& elementMaterial = materials[globalPointIdx];
      rhoMean += elementMaterial.rho * quadWeight;
    }

    MaterialT result{};
    result.rho = rhoMean;

    return result;
  }
};

template <>
struct MaterialAverager<::seissol::model::AcousticMaterial> {
  constexpr static bool Supported = true;
  using MaterialT = ::seissol::model::AcousticMaterial;
  static MaterialT compute(const ::std::vector<MaterialT>& materials,
                           size_t start,
                           const ModelQuadratureRule& rule) {
    double rhoMean = 0.0;

    // Average of the bulk modulus, used for acoustic material
    double kMeanInv = 0.0;

    for (unsigned quadPointIdx = 0; quadPointIdx < rule.weights.size(); ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * rule.weights[quadPointIdx];
      const size_t globalPointIdx = start + quadPointIdx;
      const auto& elementMaterial = materials[globalPointIdx];
      rhoMean += elementMaterial.rho * quadWeight;
      kMeanInv += 1.0 / elementMaterial.lambda * quadWeight;
    }

    MaterialT result{};
    result.rho = rhoMean;

    // Harmonic average is used for mu/K, so take the reciprocal
    result.lambda = 1.0 / kMeanInv;
    result.mu = 0.0;

    return result;
  }
};

template <>
struct MaterialAverager<::seissol::model::ElasticMaterial> {
  constexpr static bool Supported = true;
  using MaterialT =
      ::seissol::model::ElasticMaterial; // TODO(David): differentiate to acoustic properly
  static MaterialT compute(const ::std::vector<MaterialT>& materials,
                           size_t start,
                           const ModelQuadratureRule& rule) {
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

    for (unsigned quadPointIdx = 0; quadPointIdx < rule.weights.size(); ++quadPointIdx) {
      // Divide by volume of reference tetrahedron (1/6)
      const double quadWeight = 6.0 * rule.weights[quadPointIdx];
      const size_t globalPointIdx = start + quadPointIdx;
      const auto& elementMaterial = materials[globalPointIdx];
      isAcoustic |= elementMaterial.mu == 0.0;
      // TODO(David or Lukas): is this ifdef really correct? Wouldn't we need to iterate twice over
      // everything first?
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

    MaterialT result{};
    result.rho = rhoMean;

    // Harmonic average is used for mu/K, so take the reciprocal
    if (isAcoustic) {
      result.lambda = 1.0 / kMeanInv;
      result.mu = 0.0;
    } else {
      const auto muMean = 1.0 / muMeanInv;
      // Derive lambda from averaged mu and (Poisson ratio / elastic modulus)
      result.lambda =
          (4.0 * ::std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);
      result.mu = muMean;
    }

    return result;
  }
};

template <::std::size_t Mechanisms>
struct MaterialAverager<::seissol::model::ViscoElasticMaterial<Mechanisms>> {
  constexpr static bool Supported = true;
  using MaterialT = ::seissol::model::ViscoElasticMaterial<Mechanisms>;
  static MaterialT compute(const ::std::vector<MaterialT>& materials,
                           size_t start,
                           const ModelQuadratureRule& rule) {
    double muMeanInv = 0.0;
    double rhoMean = 0.0;
    double vERatioMean = 0.0;
    double QpMean = 0.0;
    double QsMean = 0.0;

#pragma omp for simd
    for (unsigned quadPointIdx = 0; quadPointIdx < rule.weights.size(); ++quadPointIdx) {
      const double quadWeight = 6.0 * rule.weights[quadPointIdx];
      const size_t globalPointIdx = start + quadPointIdx;
      const auto& elementMaterial = materials[globalPointIdx];
      muMeanInv += 1.0 / elementMaterial.mu * quadWeight;
      rhoMean += elementMaterial.rho * quadWeight;
      vERatioMean +=
          elementMaterial.lambda /
          (2.0 * elementMaterial.mu * (3.0 * elementMaterial.lambda + 2.0 * elementMaterial.mu)) *
          quadWeight;
      QpMean += elementMaterial.Qp * quadWeight;
      QsMean += elementMaterial.Qs * quadWeight;
    }

    // Harmonic average is used for mu, so take the reciprocal
    double muMean = 1.0 / muMeanInv;
    // Derive lambda from averaged mu and (Poisson ratio / elastic modulus)
    double lambdaMean =
        (4.0 * ::std::pow(muMean, 2) * vERatioMean) / (1.0 - 6.0 * muMean * vERatioMean);

    MaterialT result{};
    result.rho = rhoMean;
    result.mu = muMean;
    result.lambda = lambdaMean;
    result.Qp = QpMean;
    result.Qs = QsMean;

    return result;
  }
};

// TODO(David):: move both MaterialBinder and MaterialPostprocessor to the material classes
// themselves, maybe make inheritable

template <typename MaterialT>
struct MaterialBinder {
  static void addBindingPoints(easi::ArrayOfStructsAdapter<MaterialT>& adapter) {}
};

template <>
struct MaterialBinder<::seissol::model::ElasticMaterial> {
  using MaterialT = ::seissol::model::ElasticMaterial;
  static void addBindingPoints(easi::ArrayOfStructsAdapter<MaterialT>& adapter) {
    adapter.addBindingPoint("rho", &MaterialT::rho);
    adapter.addBindingPoint("mu", &MaterialT::mu);
    adapter.addBindingPoint("lambda", &MaterialT::lambda);
  }
};

template <::std::size_t Mechanisms>
struct MaterialBinder<::seissol::model::ViscoElasticMaterial<Mechanisms>> {
  using MaterialT = ::seissol::model::ViscoElasticMaterial<Mechanisms>;
  static void addBindingPoints(easi::ArrayOfStructsAdapter<MaterialT>& adapter) {
    adapter.addBindingPoint("rho", &MaterialT::rho);
    adapter.addBindingPoint("mu", &MaterialT::mu);
    adapter.addBindingPoint("lambda", &MaterialT::lambda);
    adapter.addBindingPoint("Qp", &MaterialT::Qp);
    adapter.addBindingPoint("Qs", &MaterialT::Qs);
  }
};

template <>
struct MaterialBinder<::seissol::model::AnisotropicMaterial> {
  using MaterialT = ::seissol::model::AnisotropicMaterial;
  static void addBindingPoints(easi::ArrayOfStructsAdapter<MaterialT>& adapter) {
    adapter.addBindingPoint("rho", &MaterialT::rho);
    adapter.addBindingPoint("c11", &MaterialT::c11);
    adapter.addBindingPoint("c12", &MaterialT::c12);
    adapter.addBindingPoint("c13", &MaterialT::c13);
    adapter.addBindingPoint("c14", &MaterialT::c14);
    adapter.addBindingPoint("c15", &MaterialT::c15);
    adapter.addBindingPoint("c16", &MaterialT::c16);
    adapter.addBindingPoint("c22", &MaterialT::c22);
    adapter.addBindingPoint("c23", &MaterialT::c23);
    adapter.addBindingPoint("c24", &MaterialT::c24);
    adapter.addBindingPoint("c25", &MaterialT::c25);
    adapter.addBindingPoint("c26", &MaterialT::c26);
    adapter.addBindingPoint("c33", &MaterialT::c33);
    adapter.addBindingPoint("c34", &MaterialT::c34);
    adapter.addBindingPoint("c35", &MaterialT::c35);
    adapter.addBindingPoint("c36", &MaterialT::c36);
    adapter.addBindingPoint("c44", &MaterialT::c44);
    adapter.addBindingPoint("c45", &MaterialT::c45);
    adapter.addBindingPoint("c46", &MaterialT::c46);
    adapter.addBindingPoint("c55", &MaterialT::c55);
    adapter.addBindingPoint("c56", &MaterialT::c56);
    adapter.addBindingPoint("c66", &MaterialT::c66);
  }
};

template <>
struct MaterialBinder<::seissol::model::PoroElasticMaterial> {
  using MaterialT = ::seissol::model::PoroElasticMaterial;
  static void addBindingPoints(easi::ArrayOfStructsAdapter<MaterialT>& adapter) {
    adapter.addBindingPoint("bulk_solid", &MaterialT::bulkSolid);
    adapter.addBindingPoint("rho", &MaterialT::rho);
    adapter.addBindingPoint("lambda", &MaterialT::lambda);
    adapter.addBindingPoint("mu", &MaterialT::mu);
    adapter.addBindingPoint("porosity", &MaterialT::porosity);
    adapter.addBindingPoint("permeability", &MaterialT::permeability);
    adapter.addBindingPoint("tortuosity", &MaterialT::tortuosity);
    adapter.addBindingPoint("bulk_fluid", &MaterialT::bulkFluid);
    adapter.addBindingPoint("rho_fluid", &MaterialT::rhoFluid);
    adapter.addBindingPoint("viscosity", &MaterialT::viscosity);
  }
};

template <>
struct MaterialBinder<::seissol::model::Plasticity> {
  using MaterialT = ::seissol::model::Plasticity;
  static void addBindingPoints(easi::ArrayOfStructsAdapter<MaterialT>& adapter) {
    adapter.addBindingPoint("bulkFriction", &MaterialT::bulkFriction);
    adapter.addBindingPoint("plastCo", &MaterialT::plastCo);
    adapter.addBindingPoint("s_xx", &MaterialT::s_xx);
    adapter.addBindingPoint("s_yy", &MaterialT::s_yy);
    adapter.addBindingPoint("s_zz", &MaterialT::s_zz);
    adapter.addBindingPoint("s_xy", &MaterialT::s_xy);
    adapter.addBindingPoint("s_yz", &MaterialT::s_yz);
    adapter.addBindingPoint("s_xz", &MaterialT::s_xz);
  }
};

template <typename MaterialT>
struct MaterialPostprocess {
  static void postprocess(MaterialT& material) {}
};

template <size_t Mechanisms>
struct MaterialPostprocess<::seissol::model::ViscoElasticMaterial<Mechanisms>> {
  static void postprocess(::seissol::model::ViscoElasticMaterial<Mechanisms>& material) {
    const auto& parameters = ::seissol::SeisSol::main.getSeisSolParameters();
    ::seissol::physics::fitAttenuation(
        material, parameters.model.freqCentral, parameters.model.freqRatio);
  }
};

/*template<typename QueryMaterialT, typename OutputMaterialT>
OutputMaterialT convertMaterial(QueryMaterialT&& material);

template<typename MaterialT>
MaterialT convertMaterial(MaterialT&& material) {
    return ::std::forward<MaterialT>(material);
}*/

template <typename QueryMaterialT, typename OutputMaterialT>
auto convertMaterial(QueryMaterialT&& source)
    -> decltype(OutputMaterialT(::std::forward<QueryMaterialT>(source))) {
  return ::std::move(OutputMaterialT(::std::forward<QueryMaterialT>(source)));
}

// TODO(David): add up/downcasts here (partially done). E.g. Elastic <=> Acoustic, Elastic <=>
// Anisotropic etc.

/*template<> // NOT YET!
model::PlasticityData<> convertMaterial<model::Plasticity,
model::PlasticityData<>>(model::Plasticity&& material)
{
    return model::PlasticityData<>(::std::forward<model::Plasticity>(material));
}*/

template <typename Config>
constexpr bool supportsAveraging() {
  return !Config::Plasticity && MaterialAverager<typename Config::MaterialT>::Supported;
}

struct MaterialInfo {
  int order;
  ::std::size_t index;
  ::std::size_t uncompressedIndex;
  ::std::size_t uncompressedSize;
};

template <typename QueryMaterialT, typename OutputMaterialT, typename OutputArrayT>
void runMaterialQuery(::std::vector<OutputArrayT>& output,
                      const ::std::vector<MaterialInfo>& materialInfo,
                      const ::std::string& modelFile,
                      const initializers::CellToVertexArray& ctov,
                      const ModelQuadratureRuleCollection& rules,
                      bool postprocess) {
  const auto& last = materialInfo[materialInfo.size() - 1];
  auto uncompressedCount = last.uncompressedIndex + last.uncompressedSize;

  auto count = materialInfo.size();

  if constexpr (::std::is_same_v<QueryMaterialT, OutputMaterialT>) {
    logInfo() << "Material query for" << QueryMaterialT::Text << "(" << count << "cells,"
              << uncompressedCount << "quadrature points)";
  } else {
    logInfo() << "Material query for" << QueryMaterialT::Text << "=>" << OutputMaterialT::Text
              << "(" << count << "cells," << uncompressedCount << "quadrature points)";
  }

  easi::Query query(uncompressedCount, 3);
#pragma omp parallel for
  for (size_t i = 0; i < count; ++i) {
    auto mapped =
        rules.getRule(materialInfo[i].order).map(ctov.elementCoordinates(materialInfo[i].index));
    auto group = ctov.elementGroups(materialInfo[i].index);
    for (size_t j = 0; j < materialInfo[i].uncompressedSize; ++j) {
      query.x(j + materialInfo[i].uncompressedSize, 0) = mapped[j][0];
      query.x(j + materialInfo[i].uncompressedSize, 1) = mapped[j][1];
      query.x(j + materialInfo[i].uncompressedSize, 2) = mapped[j][2];
      query.group(j + materialInfo[i].uncompressedSize) = group;
    }
  }

  ::std::vector<QueryMaterialT> tempOutput(uncompressedCount);
  easi::ArrayOfStructsAdapter<QueryMaterialT> adapter(tempOutput.data());
  MaterialBinder<QueryMaterialT>::addBindingPoints(adapter);
  easi::Component* model = ::seissol::initializers::loadEasiModel(modelFile);
  model->evaluate(query, adapter);
  delete model;

#pragma omp parallel for
  for (size_t i = 0; i < count; ++i) {
    auto averagedMaterial =
        materialInfo[i].order == 0
            ? tempOutput[materialInfo[i].uncompressedIndex]
            : MaterialAverager<QueryMaterialT>::compute(tempOutput,
                                                        materialInfo[i].uncompressedIndex,
                                                        rules.getRule(materialInfo[i].order));
    auto outputMaterial = ::std::move(
        convertMaterial<QueryMaterialT, OutputMaterialT>(::std::move(averagedMaterial)));
    if (postprocess) {
      MaterialPostprocess<OutputMaterialT>::postprocess(outputMaterial);
    }
    output[materialInfo[i].index] = ::std::move(outputMaterial);
  }
}
} // namespace

namespace seissol::model {

::std::pair<::std::vector<SupportedMaterials>, ::std::vector<::seissol::model::Plasticity>>
    queryMaterial(const initializer::CellConfigInfoMap& infoMap,
                  const initializers::CellToVertexArray& ctov,
                  bool fullCompute) {
  logInfo(MPI::mpi.rank()) << "Beginning material queries...";

  const auto& parameters = ::seissol::SeisSol::main.getSeisSolParameters();
  const bool averaging = !parameters.model.useCellHomogenizedMaterial;
  size_t count = 0;
  size_t matcount = 0;

  ::std::unordered_map<int, ::std::size_t> compressedInfoMap;
  ::std::unordered_map<int, int> configOrder;
  ::std::unordered_map<int, bool> configPlasticity;
  ::std::size_t infoMapIndex = 0;

  for (const auto& [groupArg, info] : infoMap) {
    // TODO(C++20): assignment needed for C++17 compatibility (once C++20 hits, it can be removed)
    int group = groupArg;
    ::std::visit(
        [&](const auto& config) {
          using Config = ::std::decay_t<decltype(config)>;
          if (averaging && supportsAveraging<Config>()) {
            // averaging case
            configOrder[group] = Config::ConvergenceOrder;
          } else {
            // barycentric case
            configOrder[group] = 0;
          }

          configPlasticity[group] = Config::Plasticity;
        },
        info.config);
    // (for now, no compression; i.e. different groups with the same material get computed in two
    // passes)
    compressedInfoMap[group] = infoMapIndex;
    ++infoMapIndex;
  }

  int maxOrder = 0;
  ::std::vector<::std::vector<MaterialInfo>> materialInfo(compressedInfoMap.size());
  ::std::vector<::std::vector<MaterialInfo>> plasticityInfo(compressedInfoMap.size());
  for (size_t i = 0; i < ctov.size; ++i) {
    auto group = ctov.elementGroups(i);
    auto index = compressedInfoMap[group];
    int order = configOrder[group];
    {
      const auto& last = materialInfo[index].at(materialInfo[index].size() - 1);
      MaterialInfo info{order,
                        i,
                        last.uncompressedIndex + last.uncompressedSize,
                        order > 0 ? initializers::numQuadPoints(order) : 1};
      materialInfo[index].push_back(info);
    }
    if (configPlasticity[group] && fullCompute) {
      const auto& last = plasticityInfo[index].at(plasticityInfo[index].size() - 1);
      MaterialInfo info{order,
                        i,
                        last.uncompressedIndex + last.uncompressedSize,
                        order > 0 ? initializers::numQuadPoints(order) : 1};
      plasticityInfo[index].push_back(info);
    }
    maxOrder = ::std::max(order, maxOrder);
  }

  ModelQuadratureRuleCollection rules(maxOrder);
  ::std::vector<SupportedMaterials> materialOutput(ctov.size);
  ::std::vector<Plasticity> plasticityOutput(fullCompute ? ctov.size : 0);

  for (::std::size_t i = 0; i < compressedInfoMap.size(); ++i) {
    const auto& materials = materialInfo[i];
    const auto& plasticities = plasticityInfo[i];
    if (!materials.empty()) {
      int exampleGroup = ctov.elementGroups(materials[0].index);
      const auto& exampleInfo = infoMap.at(exampleGroup);
      ::std::visit(
          [&](const auto& config) {
            using Config = ::std::decay_t<decltype(config)>;
            using MaterialT = typename Config::MaterialT;
            // TODO(David): add functionality to use a different easi file material than is
            // outputted (by casting up/down)
            runMaterialQuery<MaterialT, MaterialT>(
                materialOutput, materials, exampleInfo.model, ctov, rules, fullCompute);
            if (fullCompute) {
              runMaterialQuery<Plasticity, Plasticity>(
                  plasticityOutput, plasticities, exampleInfo.model, ctov, rules, fullCompute);
            }
          },
          exampleInfo.config);
    }
  }

  return {materialOutput, plasticityOutput};
}

} // namespace seissol::model
