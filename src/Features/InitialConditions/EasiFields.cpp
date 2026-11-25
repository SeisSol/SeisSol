// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EasiFields.h"

#include "Config.h"
#include "Core/Common/Constants.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Easi/EasiContext.h"
#include "Core/Numerical/Quadrature.h"
#include "WavePropagation/Physics/Datastructures.h"

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
  EasiLoader(bool hasTime,
             const std::vector<std::string>& files,
             const seissol::data::EasiContext& easi)
      : hasTime(hasTime), asagiReader(seissol::data::makeAsagiReader(easi)) {
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

template <typename Cfg>
std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime,
                                      const data::EasiContext& easi) {
  using MaterialT = model::MaterialOf<Cfg>;
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const int dimensions = needsTime ? (Cell::Dim + 1) : Cell::Dim;
  const int spaceStart = needsTime ? 1 : 0;
  const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(QuadPolyDegree);
  const auto& quadraturePoints = rule.first;

  // a query of its own for every file: easi maps the points of a query in place as it evaluates it
  const auto makeQuery = [&] {
    easi::Query query(elements.size() * NumQuadPoints, dimensions);
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
    return query;
  };

  std::vector<double> data(NumQuadPoints * iniFields.size() * MaterialT::Quantities.size() *
                           elements.size());
  const auto dataPointStride = iniFields.size() * MaterialT::Quantities.size();
  {
    auto models = EasiLoader(needsTime, iniFields, easi);
    for (std::size_t i = 0; i < iniFields.size(); ++i) {
      auto adapter = easi::ArraysAdapter();
      for (std::size_t j = 0; j < MaterialT::Quantities.size(); ++j) {
        const auto& quantity = MaterialT::Quantities.at(j);
        const std::size_t bindOffset = i + j * iniFields.size();
        adapter.addBindingPoint(quantity, data.data() + bindOffset, dataPointStride);
      }
      try {
        auto query = makeQuery();
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
  template std::vector<double> projectEasiFields<Cfg>(const std::vector<std::string>&,             \
                                                      double,                                      \
                                                      const seissol::geometry::MeshReader&,        \
                                                      bool,                                        \
                                                      const data::EasiContext&);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

} // namespace seissol::initializer
