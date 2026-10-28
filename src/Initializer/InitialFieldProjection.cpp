// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "InitialFieldProjection.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/MeshReader.h"
#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Numerical/Quadrature.h"
#include "ParameterDB.h"
#include "Physics/InitialField.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cstddef>
#include <easi/Query.h>
#include <easi/ResultAdapter.h>
#include <easi/YAMLParser.h>
#include <exception>
#include <memory>
#include <string>
#include <utils/logger.h>
#include <vector>

#ifdef USE_ASAGI
#include "Common/Real.h"
#include "Reader/AsagiReader.h"

#include <easi/util/AsagiReader.h>
#endif

// time-dependent conditions require easi version 1.5.0 or higher
#ifdef EASI_VERSION_MAJOR
#if EASI_VERSION_MINOR >= 5 || EASI_VERSION_MAJOR > 1
#define SUPPORTS_EASI_TIME
#endif
#endif

#ifdef SUPPORTS_EASI_TIME
#include <set>
#endif

GENERATE_HAS_MEMBER(Qane)

#ifndef USE_ASAGI
namespace easi {
class AsagiReader {};
} // namespace easi
#endif

namespace {
struct EasiLoader {
  bool hasTime;
  std::vector<std::unique_ptr<easi::Component>> components;
  std::unique_ptr<easi::AsagiReader> asagiReader;
  std::unique_ptr<easi::YAMLParser> parser;
  EasiLoader(bool hasTime, const std::vector<std::string>& files) : hasTime(hasTime) {
#ifdef USE_ASAGI
    asagiReader = std::make_unique<seissol::asagi::AsagiReader>();
#else
    asagiReader.reset();
#endif

    // NOTE: easi currently sorts the dimension names lexicographically (due to using std::set)
    // hence: if we have time as a dimension, it will come first
#ifdef SUPPORTS_EASI_TIME
    const auto dimensionNames =
        hasTime ? std::set<std::string>{"t", "x", "y", "z"} : std::set<std::string>{"x", "y", "z"};
    parser = std::make_unique<easi::YAMLParser>(dimensionNames, asagiReader.get());
#else
    // ignore time
    if (hasTime) {
      logError() << "easi is too old for time-dependent initial conditions. You need at least "
                    "version 1.5.0.";
    }
    parser = std::make_unique<easi::YAMLParser>(3, asagiReader.get(), 'x');
#endif
    components.resize(files.size());
    for (std::size_t i = 0; i < files.size(); ++i) {
      try {
        components[i] = std::unique_ptr<easi::Component>(parser->parse(files.at(i)));
      } catch (const std::exception& error) {
        logError() << "Error while parsing easi file" << files.at(i) << ":"
                   << std::string(error.what());
      }
    }
  }
};
} // namespace

namespace seissol::initializer {

namespace {

/// Projects the initial fields onto the cells of `layer`, which compute in the configuration `Cfg`.
template <typename Cfg>
void projectInitialFieldOnLayer(
    const std::vector<std::unique_ptr<physics::InitialField>>& iniFields,
    const seissol::geometry::MeshReader& meshReader,
    LTS::Layer& layer) {
  using real = Real<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // Looked up rather than named: a configuration without anelastic unknowns has no Qane.
  const auto* anelasticLayout = runtime::tensorTable(Variant).find("Qane", {});

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(QuadPolyDegree);
  const auto& quadraturePoints = rule.first;

#if !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
  {
    alignas(Alignment) real iniCondData[tensor::iniCond<Cfg>::size()] = {};
    auto iniCond = init::iniCond<Cfg>::view::create(iniCondData);

    std::vector<std::array<double, Cell::Dim>> quadraturePointsXyz;
    quadraturePointsXyz.resize(NumQuadPoints);

    runtime::kernel::projectIniCond krnl;
    krnl.iniCond = runtime::init::iniCond::view(Variant, iniCondData);

    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    const auto* material = layer.var<LTS::Material>();
    auto* dofs = layer.var<LTS::Dofs>(Cfg());
    auto* dofsAne = layer.var<LTS::DofsAne>(Cfg());

#if !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const auto meshId = secondaryInformation[cell].meshId;
      const auto transform = seissol::geometry::AffineTransform::fromMeshCell(meshId, meshReader);
      transform.refToSpace(
          quadraturePoints.data(), quadraturePointsXyz.data(), quadraturePoints.size());

      const CellMaterialData& materialData = material[cell];
      for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
        auto sub = multisim::simtensor<Cfg>(iniCond, s);
        iniFields[s % iniFields.size()]->evaluate(
            0.0, quadraturePointsXyz.data(), quadraturePointsXyz.size(), materialData, sub);
      }

      krnl.Q = runtime::init::Q::view(Variant, dofs[cell]);
      if constexpr (kernels::HasSize<tensor::Qane<Cfg>>::Value) {
        set_Qane(krnl, yateto::viewOf(anelasticLayout, dofsAne[cell]));
      }
      krnl.execute(Variant);
    }
  }
}

/// Projects the values `data` of the easi fields, as `projectEasiFields<Cfg>` gives them, onto the
/// cells of `layer`, which compute in the configuration `Cfg`.
template <typename Cfg>
void projectEasiFieldsOnLayer(const std::vector<double>& data,
                              std::size_t fieldCount,
                              LTS::Layer& layer) {
  using real = Real<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // Looked up rather than named: a configuration without anelastic unknowns has no Qane.
  const auto* anelasticLayout = runtime::tensorTable(Variant).find("Qane", {});
  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto quantityCount = model::MaterialOf<Cfg>::Quantities.size();
  const auto dataStride = NumQuadPoints * fieldCount * quantityCount;

#if !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
  {
    alignas(Alignment) real iniCondData[tensor::iniCond<Cfg>::size()] = {};
    auto iniCond = init::iniCond<Cfg>::view::create(iniCondData);

    runtime::kernel::projectIniCond krnl;
    krnl.iniCond = runtime::init::iniCond::view(Variant, iniCondData);

    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* dofs = layer.var<LTS::Dofs>(Cfg());
    auto* dofsAne = layer.var<LTS::DofsAne>(Cfg());

#if !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const auto meshId = secondaryInformation[cell].meshId;
      // TODO: multisim loop

      for (std::size_t s = 0; s < Cfg::NumSimulations; s++) {
        auto sub = multisim::simtensor<Cfg>(iniCond, s);
        for (std::size_t i = 0; i < NumQuadPoints; ++i) {
          for (std::size_t j = 0; j < quantityCount; ++j) {
            sub(i, j) = data.at(meshId * dataStride + quantityCount * i + j);
          }
        }
      }

      krnl.Q = runtime::init::Q::view(Variant, dofs[cell]);
      if constexpr (kernels::HasSize<tensor::Qane<Cfg>>::Value) {
        set_Qane(krnl, yateto::viewOf(anelasticLayout, dofsAne[cell]));
      }
      krnl.execute(Variant);
    }
  }
}

} // namespace

void projectInitialField(
    const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& iniFields,
    const seissol::geometry::MeshReader& meshReader,
    LTS::Storage& storage) {
  for (auto& layer : storage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      projectInitialFieldOnLayer<decltype(cfg)>(
          iniFields.at(layer.getIdentifier().config), meshReader, layer);
    });
  }
}

template <typename Cfg>
std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime) {
  using MaterialT = model::MaterialOf<Cfg>;
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const int dimensions = needsTime ? (Cell::Dim + 1) : Cell::Dim;
  const int spaceStart = needsTime ? 1 : 0;
  easi::Query query(elements.size() * NumQuadPoints, dimensions);

  {
    const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(QuadPolyDegree);
    const auto& quadraturePoints = rule.first;

#pragma omp parallel for schedule(static)
    for (std::size_t elem = 0; elem < elements.size(); ++elem) {
      const auto transform = seissol::geometry::AffineTransform::fromMeshCell(elem, meshReader);
      for (size_t i = 0; i < NumQuadPoints; ++i) {
        const auto transformed = transform.refToSpace(quadraturePoints[i]);
        for (std::size_t d = 0; d < Cell::Dim; ++d) {
          query.x(elem * NumQuadPoints + i, spaceStart + d) = transformed[d];
        }
        if (needsTime) {
          query.x(elem * NumQuadPoints + i, 0) = time;
        }
        query.group(elem * NumQuadPoints + i) = elements[elem].group;
      }
    }
  }

  std::vector<double> data(NumQuadPoints * iniFields.size() * MaterialT::Quantities.size() *
                           elements.size());
  const auto dataPointStride = iniFields.size() * MaterialT::Quantities.size();
  {
    auto models = EasiLoader(needsTime, iniFields);
    for (std::size_t i = 0; i < iniFields.size(); ++i) {
      auto adapter = easi::ArraysAdapter();
      for (std::size_t j = 0; j < MaterialT::Quantities.size(); ++j) {
        const auto& quantity = MaterialT::Quantities.at(j);
        const std::size_t bindOffset = i + j * iniFields.size();
        adapter.addBindingPoint(quantity, data.data() + bindOffset, dataPointStride);
      }
      try {
        models.components.at(i)->evaluate(query, adapter);
      } catch (const std::exception& error) {
        logError() << "Error while applying easi file" << iniFields.at(i) << ":"
                   << std::string(error.what());
      }
    }
  }

  return data;
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template std::vector<double> projectEasiFields<Cfg>(                                             \
      const std::vector<std::string>&, double, const seissol::geometry::MeshReader&, bool);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

void projectEasiInitialField(const std::vector<std::string>& iniFields,
                             const seissol::geometry::MeshReader& meshReader,
                             LTS::Storage& storage,
                             bool needsTime) {
  // the fields are sampled at the points of each configuration, once for all of its layers
  for (const auto config : storage.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      const auto data = projectEasiFields<Cfg>(iniFields, 0, meshReader, needsTime);
      for (auto& layer : storage.leaves(Ghost)) {
        if (layer.getIdentifier().config == config) {
          projectEasiFieldsOnLayer<Cfg>(data, iniFields.size(), layer);
        }
      }
    });
  }
}

} // namespace seissol::initializer
