// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ImposedSlipRatesInitializer.h"

#include "Common/ConfigDispatch.h"
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/SlipRateScript.h"
#include "DynamicRupture/Misc.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/ParameterDB.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Reader/Scripting/DataTable.h"
#include "SeisSol.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::initializer {

namespace {

/// The rotation from strike and dip into the coordinate system of the face: slip along strike s
/// and dip d is cos * s + sin * d along the first direction of the face, cos * d - sin * s along
/// the second.
std::pair<double, double> slipRotation(const Fault& fault) {
  CoordinateT strikeVector{};
  CoordinateT dipVector{};
  misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

  // cos^2 can be greater than 1 because of rounding errors
  const double cos = std::clamp(MeshTools::dot(strikeVector, fault.tangent1), -1.0, 1.0);
  CoordinateT crossProduct{};
  MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
  const double scalarProduct = MeshTools::dot(crossProduct, fault.normal);
  const double sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
  return {cos, sin};
}

/// The name of the coordinate a row of the script reads.
std::string coordinateName(friction_law::SlipRateScript::Row kind) {
  switch (kind) {
  case friction_law::SlipRateScript::Row::X:
    return "x";
  case friction_law::SlipRateScript::Row::Y:
    return "y";
  default:
    return "z";
  }
}

/// No initial stress and no nucleation: the slip is imposed.
template <typename Cfg>
void clearStressSources(DynamicRupture::Layer& layer, std::uint32_t sourceCount) {
  auto* stressInFaultCS = layer.var<DynamicRupture::StressSourceInFaultCS>(Cfg());
  auto* pressure = layer.var<DynamicRupture::StressSourcePressure>(Cfg());
  auto* onset = layer.var<DynamicRupture::StressSourceOnset>(Cfg());
  auto* riseTime = layer.var<DynamicRupture::StressSourceRiseTime>(Cfg());
  for (std::uint32_t source = 0; source < sourceCount; ++source) {
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
        for (std::uint32_t dim = 0; dim < 6; ++dim) {
          stressInFaultCS[ltsFace * sourceCount + source][dim][pointIndex] = 0;
        }
        pressure[ltsFace * sourceCount + source][pointIndex] = 0;
        onset[ltsFace * sourceCount + source][pointIndex] = 0;
        riseTime[ltsFace * sourceCount + source][pointIndex] = 0;
      }
    }
  }
}

} // namespace

void ImposedSlipRatesInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  logQuadratureRules(drStorage);
  for (auto& layer : drStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
      // parameters to be read from fault parameters yaml file
      std::unordered_map<std::string, void*> parameterToStorageMap;

      auto* imposedSlipDirection1 = layer.var<LTSImposedSlipRates::ImposedSlipDirection1>(Cfg());
      auto* imposedSlipDirection2 = layer.var<LTSImposedSlipRates::ImposedSlipDirection2>(Cfg());
      auto* onsetTime = layer.var<LTSImposedSlipRates::OnsetTime>(Cfg());

      // First read slip in strike/dip direction. Later we will rotate this to the face aligned
      // coordinate system.
      using VectorOfArraysT = std::vector<std::array<real, misc::NumPaddedPoints<Cfg>>>;
      VectorOfArraysT strikeSlip(layer.size());
      VectorOfArraysT dipSlip(layer.size());
      parameterToStorageMap.insert({"strike_slip", strikeSlip.data()->data()});
      parameterToStorageMap.insert({"dip_slip", dipSlip.data()->data()});
      parameterToStorageMap.insert({"rupture_onset", reinterpret_cast<real*>(onsetTime)});

      // get additional parameters (for derived friction laws)
      addAdditionalParameters(parameterToStorageMap, layer);

      for (std::size_t i = 0; i < Cfg::NumSimulations; ++i) {
        seissol::initializer::FaultParameterDB<real> faultParameterDB(i, Cfg::NumSimulations);
        // read parameters from yaml file
        for (const auto& parameterStoragePair : parameterToStorageMap) {
          faultParameterDB.addParameter(parameterStoragePair.first,
                                        static_cast<real*>(parameterStoragePair.second));
        }
        const auto faceIDs = getFaceIDsInIterator(layer);
        queryModel<Cfg>(faultParameterDB, faceIDs, i);
      }

      rotateSlipToFaultCS<Cfg>(
          layer, strikeSlip, dipSlip, imposedSlipDirection1, imposedSlipDirection2);

      clearStressSources<Cfg>(layer, stressSourceCount(*drParameters_));

      // Set initial and nucleation stress to zero, these are not needed for this FL

      fixInterpolatedSTFParameters(layer);

      initializeOtherVariables(layer);
    });
  }
}

template <typename Cfg>
void ImposedSlipRatesInitializer::rotateSlipToFaultCS(
    DynamicRupture::Layer& layer,
    const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& strikeSlip,
    const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& dipSlip,
    Real<Cfg> (*imposedSlipDirection1)[misc::NumPaddedPoints<Cfg>],
    Real<Cfg> (*imposedSlipDirection2)[misc::NumPaddedPoints<Cfg>]) {
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const auto meshFace = drFaceInformation[ltsFace].meshFace;
    const Fault& fault = seissolInstance_.meshReader().getFault().at(meshFace);
    const auto [cos, sin] = slipRotation(fault);

    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      imposedSlipDirection1[ltsFace][pointIndex] =
          cos * strikeSlip[ltsFace][pointIndex] + sin * dipSlip[ltsFace][pointIndex];
      imposedSlipDirection2[ltsFace][pointIndex] =
          -sin * strikeSlip[ltsFace][pointIndex] + cos * dipSlip[ltsFace][pointIndex];
    }
  }
}

void ImposedSlipRatesInitializer::fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) {
  // do nothing
}

void ImposedSlipRatesYoffeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*tauS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauS>(Cfg());
    real(*tauR)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauR>(Cfg());
    parameterToStorageMap.insert({"tau_S", reinterpret_cast<real*>(tauS)});
    parameterToStorageMap.insert({"tau_R", reinterpret_cast<real*>(tauR)});
  });
}

void ImposedSlipRatesYoffeInitializer::fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*tauS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauS>(Cfg());
    real(*tauR)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSImposedSlipRatesYoffe::TauR>(Cfg());
    // ensure that tauR is larger than tauS and that tauS and tauR are greater than 0 (the contrary
    // can happen due to ASAGI interpolation)
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
        tauS[ltsFace][pointIndex] = std::max(static_cast<real>(0.0), tauS[ltsFace][pointIndex]);
        tauR[ltsFace][pointIndex] = std::max(tauR[ltsFace][pointIndex], tauS[ltsFace][pointIndex]);
      }
    }
  });
}

void ImposedSlipRatesGaussianInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*riseTime)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSImposedSlipRatesGaussian::RiseTime>(Cfg());
    parameterToStorageMap.insert({"rupture_rise_time", reinterpret_cast<real*>(riseTime)});
  });
}

void ImposedSlipRatesDeltaInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {}

void ImposedSlipRatesScriptInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  using Row = friction_law::SlipRateScript::Row;
  logQuadratureRules(drStorage);
  const auto& rowNames = script_->rowNames();
  const auto& rowKinds = script_->rowKinds();
  for (std::size_t row = 0; row < rowNames.size(); ++row) {
    if (rowKinds[row] == Row::Parameter && faultParameterNames_.count(rowNames[row]) == 0) {
      logError() << "imposed slip rates: the script" << script_->path() << "reads"
                 << rowNames[row].c_str()
                 << ", which is neither t, dt, x, y, z nor sim, nor given by the fault "
                    "parameter file"
                 << drParameters_->faultFileName;
    }
  }

  for (auto& layer : drStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      constexpr std::size_t PointsPerFace = misc::NumPaddedPoints<Cfg>;
      constexpr std::size_t PointsPerSimulation = misc::NumPaddedPointsSingleSim<Cfg>;
      constexpr std::size_t Simulations = Cfg::NumSimulations;
      const std::size_t faces = layer.size();
      const std::size_t points = faces * PointsPerFace;
      // the rows, one after the other over all points of the layer (see SlipRateEvaluator)
      auto* rows =
          reinterpret_cast<double*>(layer.var<LTSImposedSlipRatesScript::ScriptParameters>(Cfg()));
      std::fill_n(rows, rowNames.size() * points, 0.0);
      std::fill_n(reinterpret_cast<Real<Cfg>*>(
                      layer.var<LTSImposedSlipRatesScript::ScriptSlipRates>(Cfg())),
                  2 * misc::TimeSteps<Cfg> * points,
                  static_cast<Real<Cfg>>(0));

      if (faces > 0) {
        const auto faceIDs = getFaceIDsInIterator(layer);
        const bool readsParameters =
            std::find(rowKinds.begin(), rowKinds.end(), Row::Parameter) != rowKinds.end();
        for (std::size_t sim = 0; readsParameters && sim < Simulations; ++sim) {
          seissol::initializer::FaultParameterDB<double> faultParameterDB(sim, Simulations);
          for (std::size_t row = 0; row < rowNames.size(); ++row) {
            if (rowKinds[row] == Row::Parameter) {
              faultParameterDB.addParameter(rowNames[row], rows + row * points);
            }
          }
          const auto& fileName = drParameters_->faultFileNames[sim].has_value()
                                     ? drParameters_->faultFileNames[sim]
                                     : drParameters_->faultFileNames[0];
          const seissol::initializer::FaultGPGenerator<Cfg> queryGen(seissolInstance_.meshReader(),
                                                                     faceIDs);
          faultParameterDB.evaluateModel(fileName.value(), queryGen);
        }

        // the points as the fault parameters are queried; all simulations of a point share its
        // position
        const reader::scripting::DataTable table =
            seissol::initializer::FaultGPGenerator<Cfg>(seissolInstance_.meshReader(), faceIDs)
                .generate();
        std::vector<double> coordinate(faces * PointsPerSimulation);
        const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
        for (std::size_t row = 0; row < rowNames.size(); ++row) {
          const auto kind = rowKinds[row];
          double* values = rows + row * points;
          if (kind == Row::X || kind == Row::Y || kind == Row::Z) {
            const std::string name = coordinateName(kind);
            for (const auto& entry : table.dataEntries()) {
              if (entry.name == name) {
                entry.getValues<double>(0, coordinate.size(), coordinate.data());
              }
            }
            for (std::size_t q = 0; q < coordinate.size(); ++q) {
              const std::size_t face = q / PointsPerSimulation;
              const std::size_t point = q % PointsPerSimulation;
              for (std::size_t sim = 0; sim < Simulations; ++sim) {
                values[face * PointsPerFace + point * Simulations + sim] = coordinate[q];
              }
            }
          } else if (kind == Row::Simulation) {
            for (std::size_t p = 0; p < points; ++p) {
              values[p] = static_cast<double>(p % Simulations);
            }
          } else if (kind == Row::Cos || kind == Row::Sin) {
            for (std::size_t face = 0; face < faces; ++face) {
              const auto& fault =
                  seissolInstance_.meshReader().getFault().at(faceInformation[face].meshFace);
              const auto [cos, sin] = slipRotation(fault);
              std::fill_n(
                  values + face * PointsPerFace, PointsPerFace, kind == Row::Cos ? cos : sin);
            }
          }
        }
      }

      clearStressSources<Cfg>(layer, stressSourceCount(*drParameters_));
      initializeOtherVariables(layer);
    });
  }
}
} // namespace seissol::dr::initializer
