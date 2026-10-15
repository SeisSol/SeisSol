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
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Expr/Backend.h"
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
#include "Model/CommonDatastructures.h"
#include "Model/ScriptParameters.h"
#include "Model/StateNames.h"
#include "Numerical/Quadrature.h"
#include "Physics/InitialField.h"
#include "Physics/ScriptScope.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

GENERATE_HAS_MEMBER(Qane)

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

/// Projects the values `data` of a scripted field, as `projectScriptFields<Cfg>` gives them, onto
/// the cells of `layer`, which compute in the configuration `Cfg`.
template <typename Cfg>
void projectScriptFieldsOnLayer(const std::vector<double>& data, LTS::Layer& layer) {
  using real = Real<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // Looked up rather than named: a configuration without anelastic unknowns has no Qane.
  const auto* anelasticLayout = runtime::tensorTable(Variant).find("Qane", {});
  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto componentCount = model::stateNames<model::MaterialOf<Cfg>>().size();

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

      for (std::size_t s = 0; s < Cfg::NumSimulations; s++) {
        auto sub = multisim::simtensor<Cfg>(iniCond, s);
        // the state of the solver may leave out memory variables that it keeps apart
        const auto columns = std::min<std::size_t>(componentCount, sub.shape(1));
        for (std::size_t i = 0; i < NumQuadPoints; ++i) {
          const auto* point =
              data.data() +
              ((meshId * NumQuadPoints + i) * Cfg::NumSimulations + s) * componentCount;
          for (std::size_t j = 0; j < columns; ++j) {
            sub(i, j) = point[j];
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
std::vector<double> projectScriptFields(const std::string& path,
                                        double time,
                                        const seissol::geometry::MeshReader& meshReader,
                                        const LTS::Storage& storage,
                                        bool needsTime) {
  const auto& elements = meshReader.getElements();
  const auto names = model::stateNames<model::MaterialOf<Cfg>>();

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const std::size_t numPoints = elements.size() * NumQuadPoints;

  // The point set is materialised once and bound as plain strided views, so a compiled program
  // can read it through raw pointers instead of a per-point callback.
  std::vector<std::array<double, Cell::Dim>> points(numPoints);
  std::vector<std::int32_t> groups(numPoints);
  {
    const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(QuadPolyDegree);
    const auto& quadraturePoints = rule.first;

#pragma omp parallel for schedule(static)
    for (std::size_t elem = 0; elem < elements.size(); ++elem) {
      const auto transform = seissol::geometry::AffineTransform::fromMeshCell(elem, meshReader);
      for (size_t i = 0; i < NumQuadPoints; ++i) {
        const auto transformed = transform.refToSpace(quadraturePoints[i]);
        for (std::size_t d = 0; d < Cell::Dim; ++d) {
          points[elem * NumQuadPoints + i][d] = transformed[d];
        }
        groups[elem * NumQuadPoints + i] = elements[elem].group;
      }
    }
  }

  const auto inVars = needsTime ? std::vector<std::string>{"t", "x", "y", "z"}
                                : std::vector<std::string>{"x", "y", "z"};

  // a script that compiles runs compiled, with what it computes in this configuration of the run;
  // easi, and Lua that does not trace, interpreted
  std::unique_ptr<reader::scripting::DataReader> reader;
  if (auto program = reader::scripting::buildProgram(path)) {
    reader::scripting::CompiledReaderOptions options;
    options.backend.preferred = expr::BackendKind::RtcCpu;
    reader = std::make_unique<reader::scripting::CompiledReader>(
        physics::fieldScope(storage.configs()).restrict(*program, configIdOf<Cfg>(), path, false),
        nullptr,
        options);
  } else {
    reader = reader::scripting::buildReader(path, inVars);
  }

  // the parameters of the material the script reads, by their position among those of the material
  const auto parameterNames = model::scriptParameterNames<model::MaterialOf<Cfg>>();
  std::vector<std::size_t> readParameters;
  for (const auto& input : reader->inputVars()) {
    if (input == "x" || input == "y" || input == "z" || input == "t" || input == "sim" ||
        input == "group") {
      continue;
    }
    const auto parameter = model::parameterIndexOf(parameterNames, input);
    if (!parameter.has_value()) {
      logError() << "The scripted field" << path << "reads" << input
                 << ", which a field of the initial condition does not have; it reads x, y, z, t, "
                    "sim, group and the parameters of the material:"
                 << model::parameterList(parameterNames);
    }
    readParameters.push_back(parameter.value_or(0));
  }

  // per parameter read, its value at every point: that of the material of the cell of the element
  std::vector<std::vector<double>> parameterColumns(
      readParameters.size(),
      std::vector<double>(numPoints, std::numeric_limits<double>::quiet_NaN()));
  if (!readParameters.empty()) {
    for (const auto& layer : storage.leaves(Ghost)) {
      if (layer.getIdentifier().config != configIdOf<Cfg>()) {
        continue;
      }
      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      const auto* material = layer.var<LTS::Material>();
#pragma omp parallel
      {
        std::vector<double> values(parameterNames.size());
#pragma omp for schedule(static)
        for (std::size_t cell = 0; cell < layer.size(); ++cell) {
          // a cell kept twice has the material of the first
          if (secondaryInformation[cell].duplicate > 0) {
            continue;
          }
          material[cell].local->scriptParameters(values.data());
          const auto first = secondaryInformation[cell].meshId * NumQuadPoints;
          for (std::size_t k = 0; k < readParameters.size(); ++k) {
            std::fill_n(parameterColumns[k].begin() + static_cast<std::ptrdiff_t>(first),
                        NumQuadPoints,
                        values[readParameters[k]]);
          }
        }
      }
    }
  }

  constexpr auto Simulations = Cfg::NumSimulations;
  std::vector<double> data(numPoints * Simulations * names.size());
  const auto pointStride = Simulations * names.size();

  for (std::size_t simulation = 0; simulation < Simulations; ++simulation) {
    reader::scripting::DataTable table(numPoints);
    table.bindViewConst("x", reader::scripting::Direction::In, points.data()->data(), Cell::Dim, 0);
    table.bindViewConst("y", reader::scripting::Direction::In, points.data()->data(), Cell::Dim, 1);
    table.bindViewConst("z", reader::scripting::Direction::In, points.data()->data(), Cell::Dim, 2);
    table.bindViewConst("group", reader::scripting::Direction::In, groups.data());
    table.bindConstant("t", time);
    table.bindConstant("sim", static_cast<std::int32_t>(simulation));
    for (std::size_t k = 0; k < readParameters.size(); ++k) {
      table.bindViewConst(parameterNames[readParameters[k]],
                          reader::scripting::Direction::In,
                          parameterColumns[k].data());
    }
    for (std::size_t j = 0; j < names.size(); ++j) {
      table.bindView(names[j],
                     reader::scripting::Direction::Out,
                     data.data(),
                     pointStride,
                     simulation * names.size() + j);
    }
    reader->call(table);
  }

  return data;
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template std::vector<double> projectScriptFields<Cfg>(const std::string&,                        \
                                                        double,                                    \
                                                        const seissol::geometry::MeshReader&,      \
                                                        const LTS::Storage&,                       \
                                                        bool);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

void projectScriptInitialField(const std::string& path,
                               const seissol::geometry::MeshReader& meshReader,
                               LTS::Storage& storage,
                               bool needsTime) {
  // the field is sampled at the points of each configuration, once for all of its layers
  for (const auto config : storage.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      const auto data = projectScriptFields<Cfg>(path, 0, meshReader, storage, needsTime);
      for (auto& layer : storage.leaves(Ghost)) {
        if (layer.getIdentifier().config == config) {
          projectScriptFieldsOnLayer<Cfg>(data, layer);
        }
      }
    });
  }
}

} // namespace seissol::initializer
