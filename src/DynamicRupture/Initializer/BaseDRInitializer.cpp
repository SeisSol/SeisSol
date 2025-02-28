// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BaseDRInitializer.h"

#include "DynamicRupture/Misc.h"
#include "Geometry/MeshDefinition.h"
#include "Initializer/ParameterDB.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Numerical/Transformation.h"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <init.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::initializer {
void BaseDRInitializer::initializeFault(const seissol::initializer::DynamicRupture* const dynRup,
                                        seissol::initializer::LTSTree* const dynRupTree) {
  logInfo() << "Initializing Fault, using a quadrature rule with " << misc::NumBoundaryGaussPoints
            << " points.";
  seissol::initializer::FaultParameterDB faultParameterDB;
  for (auto& layer : dynRupTree->leaves(Ghost)) {
    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    // read initial stress and nucleation stress
    auto addStressesToStorageMap = [&parameterToStorageMap, &layer, this](
                                       StressTensor& initialStress, bool readNucleation) {
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
        for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
          for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
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
#ifdef USE_POROELASTIC
      if (isFaultParameterizedByTraction) {
        parameterToStorageMap.insert({identifiers[3], getRawData(initialStress.p)});
      } else {
        parameterToStorageMap.insert({identifiers[6], getRawData(initialStress.p)});
      }
#else
      for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
        for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
          initialStress.p[ltsFace][pointIndex] = 0.0;
        }
      }
#endif

      return isFaultParameterizedByTraction;
    };

    StressTensor initialStress(layer.getNumberOfCells());
    const bool initialStressParameterizedByTraction = addStressesToStorageMap(initialStress, false);
    StressTensor nucleationStress(layer.getNumberOfCells());
    const bool nucleationStressParameterizedByTraction =
        addStressesToStorageMap(nucleationStress, true);

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, layer);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, layer);
    queryModel(faultParameterDB, faceIDs);

    // rotate initial stress to fault coordinate system
    if (initialStressParameterizedByTraction) {
      rotateTractionToCartesianStress(dynRup, layer, initialStress);
    }

    auto* initialStressInFaultCS = layer.var(dynRup->initialStressInFaultCS);
    rotateStressToFaultCS(dynRup, layer, initialStressInFaultCS, initialStress);
    // rotate nucleation stress to fault coordinate system
    if (nucleationStressParameterizedByTraction) {
      rotateTractionToCartesianStress(dynRup, layer, nucleationStress);
    }
    real(*nucleationStressInFaultCS)[misc::NumPaddedPoints][6] =
        layer.var(dynRup->nucleationStressInFaultCS);
    rotateStressToFaultCS(dynRup, layer, nucleationStressInFaultCS, nucleationStress);

    auto* initialPressure = layer.var(dynRup->initialPressure);
    auto* nucleationPressure = layer.var(dynRup->nucleationPressure);
    for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        initialPressure[ltsFace][pointIndex] = initialStress.p[ltsFace][pointIndex];
        nucleationPressure[ltsFace][pointIndex] = nucleationStress.p[ltsFace][pointIndex];
      }
    }
    initializeOtherVariables(dynRup, layer);
  }
}

std::vector<unsigned> BaseDRInitializer::getFaceIDsInIterator(
    const seissol::initializer::DynamicRupture* const dynRup, seissol::initializer::Layer& layer) {
  const auto& drFaceInformation = layer.var(dynRup->faceInformation);
  std::vector<unsigned> faceIDs;
  faceIDs.reserve(layer.getNumberOfCells());
  // collect all face IDs within this lts leaf
  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    faceIDs.push_back(drFaceInformation[ltsFace].meshFace);
  }
  return faceIDs;
}

void BaseDRInitializer::queryModel(seissol::initializer::FaultParameterDB& faultParameterDB,
                                   const std::vector<unsigned>& faceIDs) {
  // create a query and evaluate the model
  const seissol::initializer::FaultGPGenerator queryGen(seissolInstance.meshReader(), faceIDs);
  faultParameterDB.evaluateModel(drParameters->faultFileName, queryGen);
}

void BaseDRInitializer::rotateTractionToCartesianStress(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer,
    StressTensor& stress) {
  // create rotation kernel
  real faultTractionToCartesianMatrixValues[init::stressRotationMatrix::size()];
  auto faultTractionToCartesianMatrixView =
      init::stressRotationMatrix::view::create(faultTractionToCartesianMatrixValues);
  dynamicRupture::kernel::rotateStress faultTractionToCartesianRotationKernel;
  faultTractionToCartesianRotationKernel.stressRotationMatrix =
      faultTractionToCartesianMatrixValues;

  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = layer.var(dynRup->faceInformation);
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
    for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const real initialTraction[init::initialStress::size()] = {stress.xx[ltsFace][pointIndex],
                                                                 stress.yy[ltsFace][pointIndex],
                                                                 stress.zz[ltsFace][pointIndex],
                                                                 stress.xy[ltsFace][pointIndex],
                                                                 stress.yz[ltsFace][pointIndex],
                                                                 stress.xz[ltsFace][pointIndex]};
      assert(std::abs(initialTraction[YY]) < 1e-15);
      assert(std::abs(initialTraction[ZZ]) < 1e-15);
      assert(std::abs(initialTraction[YZ]) < 1e-15);

      real cartesianStress[init::initialStress::size()]{};
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

void BaseDRInitializer::rotateStressToFaultCS(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer,
    real (*stressInFaultCS)[misc::NumPaddedPoints][6],
    const StressTensor& stress) {
  // create rotation kernel
  real cartesianToFaultCSMatrixValues[init::stressRotationMatrix::size()];
  auto cartesianToFaultCSMatrixView =
      init::stressRotationMatrix::view::create(cartesianToFaultCSMatrixValues);
  dynamicRupture::kernel::rotateStress cartesianToFaultCSRotationKernel;
  cartesianToFaultCSRotationKernel.stressRotationMatrix = cartesianToFaultCSMatrixValues;

  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    constexpr unsigned int NumStressComponents = 6;
    const auto& drFaceInformation = layer.var(dynRup->faceInformation);
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissolInstance.meshReader().getFault().at(meshFace);

    // now rotate the stress in cartesian coordinates to the element aligned coordinate system.
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, cartesianToFaultCSMatrixView, 0, 0);

    for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const real initialStress[init::initialStress::size()] = {stress.xx[ltsFace][pointIndex],
                                                               stress.yy[ltsFace][pointIndex],
                                                               stress.zz[ltsFace][pointIndex],
                                                               stress.xy[ltsFace][pointIndex],
                                                               stress.yz[ltsFace][pointIndex],
                                                               stress.xz[ltsFace][pointIndex]};
      real rotatedStress[init::initialStress::size()]{};
      cartesianToFaultCSRotationKernel.initialStress = initialStress;
      cartesianToFaultCSRotationKernel.rotatedStress = rotatedStress;
      cartesianToFaultCSRotationKernel.execute();
      for (unsigned int stressIndex = 0; stressIndex < NumStressComponents; ++stressIndex) {
        stressInFaultCS[ltsFace][pointIndex][stressIndex] = rotatedStress[stressIndex];
      }
    }
  }
}

void BaseDRInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  // do nothing for base friction law
}

void BaseDRInitializer::initializeOtherVariables(
    const seissol::initializer::DynamicRupture* const dynRup, seissol::initializer::Layer& layer) {
  // initialize rupture front flag
  bool(*ruptureTimePending)[misc::NumPaddedPoints] = layer.var(dynRup->ruptureTimePending);
  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      ruptureTimePending[ltsFace][pointIndex] = true;
    }
  }

  // initialize all other variables to zero
  real(*peakSlipRate)[misc::NumPaddedPoints] = layer.var(dynRup->peakSlipRate);
  real(*ruptureTime)[misc::NumPaddedPoints] = layer.var(dynRup->ruptureTime);
  real(*dynStressTime)[misc::NumPaddedPoints] = layer.var(dynRup->dynStressTime);
  real(*accumulatedSlipMagnitude)[misc::NumPaddedPoints] =
      layer.var(dynRup->accumulatedSlipMagnitude);
  real(*slip1)[misc::NumPaddedPoints] = layer.var(dynRup->slip1);
  real(*slip2)[misc::NumPaddedPoints] = layer.var(dynRup->slip2);
  real(*slipRateMagnitude)[misc::NumPaddedPoints] = layer.var(dynRup->slipRateMagnitude);
  real(*traction1)[misc::NumPaddedPoints] = layer.var(dynRup->traction1);
  real(*traction2)[misc::NumPaddedPoints] = layer.var(dynRup->traction2);

  for (unsigned int ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
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
    BaseDRInitializer::stressIdentifiers(bool readNucleation) {
  std::vector<std::string> tractionNames;
  std::vector<std::string> cartesianNames;
  if (readNucleation) {
    tractionNames = {"Tnuc_n", "Tnuc_s", "Tnuc_d"};
    cartesianNames = {"nuc_xx", "nuc_yy", "nuc_zz", "nuc_xy", "nuc_yz", "nuc_xz"};
  } else {
    tractionNames = {"T_n", "T_s", "T_d"};
    cartesianNames = {"s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"};
  }
#ifdef USE_POROELASTIC
  if (readNucleation) {
    tractionNames.push_back("nuc_p");
    cartesianNames.push_back("nuc_p");
  } else {
    tractionNames.push_back("p");
    cartesianNames.push_back("p");
  }
#endif

  bool allTractionParametersSupplied = true;
  bool allCartesianParametersSupplied = true;
  bool anyTractionParametersSupplied = false;
  bool anyCartesianParametersSupplied = false;
  for (size_t i = 0; i < 3; i++) {
    const auto b = faultProvides(tractionNames[i]);
    allTractionParametersSupplied &= b;
    anyTractionParametersSupplied |= b;
  }
  for (size_t i = 0; i < 6; i++) {
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
               << (readNucleation ? "nucleation stress." : "initial stress.")
               << "You have either not specified all parameters or an uncommom mixture of "
                  "parameters. Give either all of "
               << (readNucleation
                       ? "(nuc_xx, nuc_yy, nuc_zz, nuc_xy, nuc_yz, nuc_xz) or all of (Tnuc_n, "
                         "Tnuc_s, Tnuc_d)"
                       : "(s_xx, s_yy, s_zz, s_xy, s_yz, s_xz) or all of (T_n, T_s, T_d)")
               << ", but not a mixture";
    return {};
  }
}

} // namespace seissol::dr::initializer
