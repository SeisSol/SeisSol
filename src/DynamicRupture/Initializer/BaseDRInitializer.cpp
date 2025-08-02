// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BaseDRInitializer.h"

#include "DynamicRupture/Misc.h"
#include "GeneratedCode/kernel.h"
#include "Geometry/MeshDefinition.h"
#include "Initializer/ParameterDB.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include <Eigen/Dense>
#include <Equations/Datastructures.h>
#include <GeneratedCode/init.h>
#include <Model/CommonDatastructures.h>
#include <Solver/MultipleSimulations.h>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifndef NDEBUG
#include <cmath>
#endif

namespace seissol::dr::initializer {
void BaseDRInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  logInfo() << "Initializing Fault, using a quadrature rule with " << misc::NumBoundaryGaussPoints<Cfg>
            << " points.";
  for (auto& layer : drStorage.leaves(Ghost)) {
    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    // read initial stress and nucleation stress
    auto addStressesToStorageMap =
        [&parameterToStorageMap, &layer, this](StressTensor& initialStress, int readNucleation) {
          // return pointer to first element
          auto getRawData = [](StressTensor::VectorOfArraysT& vectorOfArrays) {
            return vectorOfArrays.data()->data();
          };
          // fault can be either initialized by traction or by cartesian stress
          // this method reads either the nucleation stress or the initial stress
          auto [identifiers, parametrization] = this->stressIdentifiers(readNucleation);
          const bool isFaultParameterizedByTraction = parametrization == Parametrization::Traction;
          if (isFaultParameterizedByTraction) {
            // only read traction in normal, strike and dip direction
            parameterToStorageMap.insert({identifiers[0], getRawData(initialStress.xx)});
            parameterToStorageMap.insert({identifiers[1], getRawData(initialStress.xy)});
            parameterToStorageMap.insert({identifiers[2], getRawData(initialStress.xz)});
            // set the rest to zero
            for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
              for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
                initialStress.yy[ltsFace][pointIndex] = 0.0;
                initialStress.zz[ltsFace][pointIndex] = 0.0;
                initialStress.yz[ltsFace][pointIndex] = 0.0;
              }
            }
          } else { // read all stress components from the parameter file
            parameterToStorageMap.insert({identifiers[0], getRawData(initialStress.xx)});
            parameterToStorageMap.insert({identifiers[1], getRawData(initialStress.yy)});
            parameterToStorageMap.insert({identifiers[2], getRawData(initialStress.zz)});
            parameterToStorageMap.insert({identifiers[3], getRawData(initialStress.xy)});
            parameterToStorageMap.insert({identifiers[4], getRawData(initialStress.yz)});
            parameterToStorageMap.insert({identifiers[5], getRawData(initialStress.xz)});
          }
          if constexpr (model::MaterialT::Type == model::MaterialType::Poroelastic) {
            if (isFaultParameterizedByTraction) {
              parameterToStorageMap.insert({identifiers[3], getRawData(initialStress.p)});
            } else {
              parameterToStorageMap.insert({identifiers[6], getRawData(initialStress.p)});
            }
          } else {
            for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
              for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
                initialStress.p[ltsFace][pointIndex] = 0.0;
              }
            }
          }

          return isFaultParameterizedByTraction;
        };

    StressTensor initialStress(layer.size());
    const bool initialStressParameterizedByTraction = addStressesToStorageMap(initialStress, 0);

    std::vector<bool> nucleationStressParameterizedByTraction(drParameters->nucleationCount);
    std::vector<StressTensor> nucleationStresses;
    nucleationStresses.reserve(drParameters->nucleationCount);
    for (unsigned i = 0; i < drParameters->nucleationCount; ++i) {
      nucleationStresses.emplace_back(layer.size());
      nucleationStressParameterizedByTraction[i] =
          addStressesToStorageMap(nucleationStresses[i], i + 1);
    }

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, layer);

    for (std::size_t i = 0; i < multisim::NumSimulations; ++i) {
      seissol::initializer::FaultParameterDB faultParameterDB(i);
      // read parameters from yaml file
      for (const auto& parameterStoragePair : parameterToStorageMap) {
        faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
      }
      const auto faceIDs = getFaceIDsInIterator(layer);
      queryModel(faultParameterDB, faceIDs, i);
    }

    // rotate initial stress to fault coordinate system
    if (initialStressParameterizedByTraction) {
      rotateTractionToCartesianStress(layer, initialStress);
    }

    auto* initialStressInFaultCS = layer.var<DynamicRupture::InitialStressInFaultCS>();
    rotateStressToFaultCS(layer, initialStressInFaultCS, 0, 1, initialStress);
    // rotate nucleation stress to fault coordinate system
    for (unsigned i = 0; i < drParameters->nucleationCount; ++i) {
      if (nucleationStressParameterizedByTraction[i]) {
        rotateTractionToCartesianStress(layer, nucleationStresses[i]);
      }
      auto* nucleationStressInFaultCS = layer.var<DynamicRupture::NucleationStressInFaultCS>();
      rotateStressToFaultCS(layer,
                            nucleationStressInFaultCS,
                            i,
                            drParameters->nucleationCount,
                            nucleationStresses[i]);
    }

    auto* initialPressure = layer.var<DynamicRupture::InitialPressure>();
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
        initialPressure[ltsFace][pointIndex] = initialStress.p[ltsFace][pointIndex];
        for (unsigned i = 0; i < drParameters->nucleationCount; ++i) {
          auto* nucleationPressure = layer.var<DynamicRupture::NucleationPressure>();
          nucleationPressure[ltsFace * drParameters->nucleationCount + i][pointIndex] =
              nucleationStresses[i].p[ltsFace][pointIndex];
        }
      }
    }

    initializeOtherVariables(layer);
  }
}

std::vector<unsigned> BaseDRInitializer::getFaceIDsInIterator(DynamicRupture::Layer& layer) {
  const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
  std::vector<unsigned> faceIDs;
  faceIDs.reserve(layer.size());
  // collect all face IDs within this lts leaf
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    faceIDs.push_back(drFaceInformation[ltsFace].meshFace);
  }
  return faceIDs;
}

void BaseDRInitializer::queryModel(seissol::initializer::FaultParameterDB& faultParameterDB,
                                   const std::vector<unsigned>& faceIDs,
                                   std::size_t simid) {
  // create a query and evaluate the model
  if (!drParameters->faultFileNames[simid].has_value()) {
    simid = 0;
  }
  {
    const seissol::initializer::FaultGPGenerator queryGen(seissolInstance.meshReader(), faceIDs);
    faultParameterDB.evaluateModel(drParameters->faultFileNames[simid].value(), queryGen);
  }
}

void BaseDRInitializer::rotateTractionToCartesianStress(DynamicRupture::Layer& layer,
                                                        StressTensor& stress) {
  // create rotation kernel
  real faultTractionToCartesianMatrixValues[init::stressRotationMatrix<Cfg>::size()];
  auto faultTractionToCartesianMatrixView =
      init::stressRotationMatrix<Cfg>::view::create(faultTractionToCartesianMatrixValues);
  dynamicRupture::kernel::rotateStress<Cfg> faultTractionToCartesianRotationKernel;
  faultTractionToCartesianRotationKernel.stressRotationMatrix =
      faultTractionToCartesianMatrixValues;

  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissolInstance.meshReader().getFault().at(meshFace);

    // if we read the traction in strike, dip and normal direction, we first transform it to stress
    // in cartesian coordinates
    VrtxCoords strike{};
    VrtxCoords dip{};
    misc::computeStrikeAndDipVectors(fault.normal, strike, dip);
    seissol::transformations::symmetricTensor2RotationMatrix(
        fault.normal, strike, dip, faultTractionToCartesianMatrixView, 0, 0);

    using namespace dr::misc::quantity_indices;
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      const real initialTraction[init::initialStress<Cfg>::size()] = {stress.xx[ltsFace][pointIndex],
                                                                 stress.yy[ltsFace][pointIndex],
                                                                 stress.zz[ltsFace][pointIndex],
                                                                 stress.xy[ltsFace][pointIndex],
                                                                 stress.yz[ltsFace][pointIndex],
                                                                 stress.xz[ltsFace][pointIndex]};
      assert(std::abs(initialTraction[YY]) < 1e-15);
      assert(std::abs(initialTraction[ZZ]) < 1e-15);
      assert(std::abs(initialTraction[YZ]) < 1e-15);

      real cartesianStress[init::initialStress<Cfg>::size()]{};
      faultTractionToCartesianRotationKernel.initialStress = initialTraction;
      faultTractionToCartesianRotationKernel.rotatedStress = cartesianStress;
      faultTractionToCartesianRotationKernel.execute();
      stress.xx[ltsFace][pointIndex] = cartesianStress[XX];
      stress.yy[ltsFace][pointIndex] = cartesianStress[YY];
      stress.zz[ltsFace][pointIndex] = cartesianStress[ZZ];
      stress.xy[ltsFace][pointIndex] = cartesianStress[XY];
      stress.yz[ltsFace][pointIndex] = cartesianStress[YZ];
      stress.xz[ltsFace][pointIndex] = cartesianStress[XZ];
    }
  }
}

void BaseDRInitializer::rotateStressToFaultCS(DynamicRupture::Layer& layer,
                                              real (*stressInFaultCS)[6][misc::NumPaddedPoints<Cfg>],
                                              std::size_t index,
                                              std::size_t count,
                                              const StressTensor& stress) {
  // create rotation kernel
  real cartesianToFaultCSMatrixValues[init::stressRotationMatrix<Cfg>::size()];
  auto cartesianToFaultCSMatrixView =
      init::stressRotationMatrix<Cfg>::view::create(cartesianToFaultCSMatrixValues);
  dynamicRupture::kernel::rotateStress<Cfg> cartesianToFaultCSRotationKernel;
  cartesianToFaultCSRotationKernel.stressRotationMatrix = cartesianToFaultCSMatrixValues;

  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    constexpr unsigned int NumStressComponents = model::MaterialT::TractionQuantities;
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissolInstance.meshReader().getFault().at(meshFace);

    // now rotate the stress in cartesian coordinates to the element aligned coordinate system.
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, cartesianToFaultCSMatrixView, 0, 0);

    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      const real initialStress[init::initialStress<Cfg>::size()] = {stress.xx[ltsFace][pointIndex],
                                                               stress.yy[ltsFace][pointIndex],
                                                               stress.zz[ltsFace][pointIndex],
                                                               stress.xy[ltsFace][pointIndex],
                                                               stress.yz[ltsFace][pointIndex],
                                                               stress.xz[ltsFace][pointIndex]};
      real rotatedStress[init::initialStress<Cfg>::size()]{};
      cartesianToFaultCSRotationKernel.initialStress = initialStress;
      cartesianToFaultCSRotationKernel.rotatedStress = rotatedStress;
      cartesianToFaultCSRotationKernel.execute();
      for (std::size_t stressIndex = 0; stressIndex < NumStressComponents; ++stressIndex) {
        stressInFaultCS[ltsFace * count + index][stressIndex][pointIndex] =
            rotatedStress[stressIndex];
      }
    }
  }
}

void BaseDRInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  // do nothing for base friction law
}

void BaseDRInitializer::initializeOtherVariables(DynamicRupture::Layer& layer) {
  // initialize rupture front flag
  bool (*ruptureTimePending)[misc::NumPaddedPoints<Cfg>] =
      layer.var<DynamicRupture::RuptureTimePending>();
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      ruptureTimePending[ltsFace][pointIndex] = true;
    }
  }

  // initialize all other variables to zero
  real(*peakSlipRate)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::PeakSlipRate>();
  real(*ruptureTime)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::RuptureTime>();
  real(*dynStressTime)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::DynStressTime>();
  real(*accumulatedSlipMagnitude)[misc::NumPaddedPoints<Cfg>] =
      layer.var<DynamicRupture::AccumulatedSlipMagnitude>();
  real(*slip1)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::Slip1>();
  real(*slip2)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::Slip2>();
  real(*slipRateMagnitude)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::SlipRateMagnitude>();
  real(*traction1)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::Traction1>();
  real(*traction2)[misc::NumPaddedPoints<Cfg>] = layer.var<DynamicRupture::Traction2>();

  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      peakSlipRate[ltsFace][pointIndex] = 0;
      ruptureTime[ltsFace][pointIndex] = 0;
      dynStressTime[ltsFace][pointIndex] = 0;
      accumulatedSlipMagnitude[ltsFace][pointIndex] = 0;
      slip1[ltsFace][pointIndex] = 0;
      slip2[ltsFace][pointIndex] = 0;
      slipRateMagnitude[ltsFace][pointIndex] = 0;
      traction1[ltsFace][pointIndex] = 0;
      traction2[ltsFace][pointIndex] = 0;
    }
  }
}

bool BaseDRInitializer::faultProvides(const std::string& parameter) {
  // TODO: Use C++20 contains
  return faultParameterNames.count(parameter) > 0;
}

std::pair<std::vector<std::string>, BaseDRInitializer::Parametrization>
    BaseDRInitializer::stressIdentifiers(int readNucleation) {
  std::vector<std::string> tractionNames;
  std::vector<std::string> cartesianNames;

  const std::string index = readNucleation > 1 ? std::to_string(readNucleation) : "";

  const auto insertIndex = [&](const std::string& name, const std::string& subscript) {
    return name + index + "_" + subscript;
  };

  if (readNucleation > 0) {
    tractionNames = {insertIndex("Tnuc", "n"), insertIndex("Tnuc", "s"), insertIndex("Tnuc", "d")};
    cartesianNames = {insertIndex("nuc", "xx"),
                      insertIndex("nuc", "yy"),
                      insertIndex("nuc", "zz"),
                      insertIndex("nuc", "xy"),
                      insertIndex("nuc", "yz"),
                      insertIndex("nuc", "xz")};
  } else {
    tractionNames = {"T_n", "T_s", "T_d"};
    cartesianNames = {"s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"};
  }
  if (model::MaterialT::Type == model::MaterialType::Poroelastic) {
    if (readNucleation > 0) {
      tractionNames.emplace_back(insertIndex("nuc", "p"));
      cartesianNames.emplace_back(insertIndex("nuc", "p"));
    } else {
      tractionNames.emplace_back("p");
      cartesianNames.emplace_back("p");
    }
  }

  bool allTractionParametersSupplied = true;
  bool allCartesianParametersSupplied = true;
  bool anyTractionParametersSupplied = false;
  bool anyCartesianParametersSupplied = false;
  for (size_t i = 0; i < 3; ++i) {
    const auto b = faultProvides(tractionNames[i]);
    allTractionParametersSupplied &= b;
    anyTractionParametersSupplied |= b;
  }
  for (size_t i = 0; i < 6; ++i) {
    const auto b = faultProvides(cartesianNames[i]);
    allCartesianParametersSupplied &= b;
    anyCartesianParametersSupplied |= b;
  }

  if (allCartesianParametersSupplied && !anyTractionParametersSupplied) {
    return {cartesianNames, Parametrization::Cartesian};
  } else if (allTractionParametersSupplied && !anyCartesianParametersSupplied) {
    return {tractionNames, Parametrization::Traction};
  } else {
    logError() << "Please specify a correct parametrization of the "
               << (readNucleation > 0
                       ? "nucleation stress " + std::to_string(readNucleation + 1) + "."
                       : "initial stress.")
               << "You have either not specified all parameters or an uncommom mixture of "
                  "parameters. Give either all of "
               << cartesianNames << "or all of" << tractionNames << ", but not a mixture of them.";
    return {};
  }
}

} // namespace seissol::dr::initializer
