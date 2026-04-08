// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BaseDRInitializer.h"

#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Initializer/ParameterDB.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <array>
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

namespace {

/**
 * Stores the initialStresses.
 */
struct StressTensor {
  explicit StressTensor(size_t size) {
    xx.resize(size);
    yy.resize(size);
    zz.resize(size);
    xy.resize(size);
    yz.resize(size);
    xz.resize(size);
    p.resize(size);
  }
  using VectorOfArraysT = std::vector<std::array<real, misc::NumPaddedPoints>>;
  VectorOfArraysT xx;
  VectorOfArraysT yy;
  VectorOfArraysT zz;
  VectorOfArraysT xy;
  VectorOfArraysT yz;
  VectorOfArraysT xz;
  VectorOfArraysT p;
};

/**
 * Rotates the fault-aligned traction to cartesian stress coordinates
 * @param layer reference to a Storage layer
 * @param stress reference to a StressTensor
 * IN: stores traction in fault strike/dip coordinate system OUT: stores the the stress in
 * cartesian coordinates
 */
void rotateTractionToCartesianStress(DynamicRupture::Layer& layer,
                                     StressTensor& stress,
                                     const geometry::MeshReader& mesh) {
  // create rotation kernel
  std::array<double, seissol::general::init::stressRotationMatrix::size()>
      faultTractionToCartesianMatrixValues{};
  auto faultTractionToCartesianMatrixView =
      seissol::general::init::stressRotationMatrix::view::create(
          faultTractionToCartesianMatrixValues.data());
  seissol::general::dynamicRupture::kernel::rotateStress faultTractionToCartesianRotationKernel;
  faultTractionToCartesianRotationKernel.stressRotationMatrix =
      faultTractionToCartesianMatrixValues.data();

  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const auto meshFace = drFaceInformation[ltsFace].meshFace;
    const Fault& fault = mesh.getFault().at(meshFace);

    // if we read the traction in strike, dip and normal direction, we first transform it to stress
    // in cartesian coordinates
    VrtxCoords strike{};
    VrtxCoords dip{};
    misc::computeStrikeAndDipVectors(fault.normal, strike, dip);
    seissol::transformations::symmetricTensor2RotationMatrix(
        fault.normal, strike, dip, faultTractionToCartesianMatrixView, 0, 0);

    using namespace dr::misc::quantity_indices;
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const std::array<double, seissol::general::init::initialStress::size()> initialTraction{
          stress.xx[ltsFace][pointIndex],
          stress.yy[ltsFace][pointIndex],
          stress.zz[ltsFace][pointIndex],
          stress.xy[ltsFace][pointIndex],
          stress.yz[ltsFace][pointIndex],
          stress.xz[ltsFace][pointIndex]};
      assert(std::abs(initialTraction[YY]) < 1e-15);
      assert(std::abs(initialTraction[ZZ]) < 1e-15);
      assert(std::abs(initialTraction[YZ]) < 1e-15);

      std::array<double, seissol::general::init::initialStress::size()> cartesianStress{};
      faultTractionToCartesianRotationKernel.initialStress = initialTraction.data();
      faultTractionToCartesianRotationKernel.rotatedStress = cartesianStress.data();
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

/**
 * Rotates the stress tensor to a fault aligned coordinate system and stores it in stressInFaultCS
 * @param layer reference to a Storage layer
 * @param stressInFaultCS pointer to array of size [numCells][6][numPaddedPoints], stores rotated
 * stress
 * @param index stress index per cell (set to 0, unless initializing multi-nucleation)
 * @param count stress count per cell (set to 1, unless initializing multi-nucleation)
 * @param stress reference to a StressTensor, stores the stress in cartesian coordinates
 */
void rotateStressToFaultCS(DynamicRupture::Layer& layer,
                           real (*stressInFaultCS)[6][misc::NumPaddedPoints],
                           std::size_t index,
                           std::size_t count,
                           const StressTensor& stress,
                           const geometry::MeshReader& mesh) {
  // create rotation kernel
  std::array<double, seissol::general::init::stressRotationMatrix::size()>
      cartesianToFaultCSMatrixValues{};
  auto cartesianToFaultCSMatrixView = seissol::general::init::stressRotationMatrix::view::create(
      cartesianToFaultCSMatrixValues.data());
  seissol::general::dynamicRupture::kernel::rotateStress cartesianToFaultCSRotationKernel;
  cartesianToFaultCSRotationKernel.stressRotationMatrix = cartesianToFaultCSMatrixValues.data();

  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    constexpr auto NumStressComponents = model::MaterialT::TractionQuantities;
    const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
    const auto meshFace = drFaceInformation[ltsFace].meshFace;
    const Fault& fault = mesh.getFault().at(meshFace);

    // now rotate the stress in cartesian coordinates to the element aligned coordinate system.
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, cartesianToFaultCSMatrixView, 0, 0);

    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const std::array<double, seissol::general::init::initialStress::size()> initialStress{
          stress.xx[ltsFace][pointIndex],
          stress.yy[ltsFace][pointIndex],
          stress.zz[ltsFace][pointIndex],
          stress.xy[ltsFace][pointIndex],
          stress.yz[ltsFace][pointIndex],
          stress.xz[ltsFace][pointIndex]};
      std::array<double, seissol::general::init::initialStress::size()> rotatedStress{};
      cartesianToFaultCSRotationKernel.initialStress = initialStress.data();
      cartesianToFaultCSRotationKernel.rotatedStress = rotatedStress.data();
      cartesianToFaultCSRotationKernel.execute();
      for (std::size_t stressIndex = 0; stressIndex < NumStressComponents; ++stressIndex) {
        stressInFaultCS[ltsFace * count + index][stressIndex][pointIndex] =
            rotatedStress[stressIndex];
      }
    }
  }
}
} // namespace

void BaseDRInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  logInfo() << "Initializing Fault, using a quadrature rule with " << misc::NumBoundaryGaussPoints
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
              for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
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
              for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
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
    for (std::uint32_t i = 0; i < drParameters->nucleationCount; ++i) {
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
      rotateTractionToCartesianStress(layer, initialStress, seissolInstance.meshReader());
    }

    auto* initialStressInFaultCS = layer.var<DynamicRupture::InitialStressInFaultCS>();
    rotateStressToFaultCS(
        layer, initialStressInFaultCS, 0, 1, initialStress, seissolInstance.meshReader());
    // rotate nucleation stress to fault coordinate system
    for (std::uint32_t i = 0; i < drParameters->nucleationCount; ++i) {
      if (nucleationStressParameterizedByTraction[i]) {
        rotateTractionToCartesianStress(layer, nucleationStresses[i], seissolInstance.meshReader());
      }
      auto* nucleationStressInFaultCS = layer.var<DynamicRupture::NucleationStressInFaultCS>();
      rotateStressToFaultCS(layer,
                            nucleationStressInFaultCS,
                            i,
                            drParameters->nucleationCount,
                            nucleationStresses[i],
                            seissolInstance.meshReader());
    }

    auto* initialPressure = layer.var<DynamicRupture::InitialPressure>();
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        initialPressure[ltsFace][pointIndex] = initialStress.p[ltsFace][pointIndex];
        for (std::uint32_t i = 0; i < drParameters->nucleationCount; ++i) {
          auto* nucleationPressure = layer.var<DynamicRupture::NucleationPressure>();
          nucleationPressure[ltsFace * drParameters->nucleationCount + i][pointIndex] =
              nucleationStresses[i].p[ltsFace][pointIndex];
        }
      }
    }

    initializeOtherVariables(layer);
  }
}

std::vector<std::size_t> BaseDRInitializer::getFaceIDsInIterator(DynamicRupture::Layer& layer) {
  const auto& drFaceInformation = layer.var<DynamicRupture::FaceInformation>();
  std::vector<std::size_t> faceIDs;
  faceIDs.reserve(layer.size());
  // collect all face IDs within this lts leaf
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    faceIDs.push_back(drFaceInformation[ltsFace].meshFace);
  }
  return faceIDs;
}

void BaseDRInitializer::queryModel(seissol::initializer::FaultParameterDB& faultParameterDB,
                                   const std::vector<std::size_t>& faceIDs,
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

void BaseDRInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  // do nothing for base friction law
}

void BaseDRInitializer::initializeOtherVariables(DynamicRupture::Layer& layer) {
  // initialize rupture front flag
  bool (*ruptureTimePending)[misc::NumPaddedPoints] =
      layer.var<DynamicRupture::RuptureTimePending>();
  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      ruptureTimePending[ltsFace][pointIndex] = true;
    }
  }

  // initialize all other variables to zero
  real(*peakSlipRate)[misc::NumPaddedPoints] = layer.var<DynamicRupture::PeakSlipRate>();
  real(*ruptureTime)[misc::NumPaddedPoints] = layer.var<DynamicRupture::RuptureTime>();
  real(*dynStressTime)[misc::NumPaddedPoints] = layer.var<DynamicRupture::DynStressTime>();
  real(*accumulatedSlipMagnitude)[misc::NumPaddedPoints] =
      layer.var<DynamicRupture::AccumulatedSlipMagnitude>();
  real(*slip1)[misc::NumPaddedPoints] = layer.var<DynamicRupture::Slip1>();
  real(*slip2)[misc::NumPaddedPoints] = layer.var<DynamicRupture::Slip2>();
  real(*slipRateMagnitude)[misc::NumPaddedPoints] = layer.var<DynamicRupture::SlipRateMagnitude>();
  real(*traction1)[misc::NumPaddedPoints] = layer.var<DynamicRupture::Traction1>();
  real(*traction2)[misc::NumPaddedPoints] = layer.var<DynamicRupture::Traction2>();

  for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
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
  std::vector<std::string> commonNames;

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
    const auto tnName = faultNameAlternatives({"T_n", "Pn0"});
    const auto tsName = faultNameAlternatives({"T_s", "Ts0"});
    const auto tdName = faultNameAlternatives({"T_d", "Td0"});

    tractionNames = {tnName, tsName, tdName};
    cartesianNames = {"s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"};
  }
  if (model::MaterialT::Type == model::MaterialType::Poroelastic) {
    if (readNucleation > 0) {
      commonNames.emplace_back(insertIndex("nuc", "p"));
    } else {
      commonNames.emplace_back("p");
    }
  }

  bool allTractionParametersSupplied = true;
  bool allCartesianParametersSupplied = true;
  bool anyTractionParametersSupplied = false;
  bool anyCartesianParametersSupplied = false;
  for (const auto& name : tractionNames) {
    const auto b = faultProvides(name);
    allTractionParametersSupplied &= b;
    anyTractionParametersSupplied |= b;
  }
  for (const auto& name : cartesianNames) {
    const auto b = faultProvides(name);
    allCartesianParametersSupplied &= b;
    anyCartesianParametersSupplied |= b;
  }
  for (const auto& name : commonNames) {
    const auto b = faultProvides(name);
    allCartesianParametersSupplied &= b;
    allTractionParametersSupplied &= b;

    // only insert common names after we checked their existence
    cartesianNames.push_back(name);
    tractionNames.push_back(name);
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

std::string BaseDRInitializer::faultNameAlternatives(const std::vector<std::string>& parameter) {
  for (const auto& name : parameter) {
    if (faultProvides(name)) {
      if (name != parameter[0]) {
        logWarning() << "You are using the deprecated fault parameter name" << name << "for"
                     << parameter[0];
      }
      return name;
    }
  }
  return parameter[0];
}

} // namespace seissol::dr::initializer
