#include "BaseDRInitializer.h"

#include "Initializer/ParameterDB.h"
#include "Model/common.hpp"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <Eigen/Dense>

namespace seissol::dr::initializers {
void BaseDRInitializer::initializeFault(seissol::initializers::DynamicRupture const* const dynRup,
                                        seissol::initializers::LTSTree* const dynRupTree) {
  const int rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Initializing Fault, using a quadrature rule with "
                << misc::numberOfBoundaryGaussPoints << " points.";
  seissol::initializers::FaultParameterDB faultParameterDB;
  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    // read initial stress and nucleation stress
    using VectorOfArraysT = std::vector<std::array<real, misc::numPaddedPoints>>;

    auto addStressesToStorageMap = [&parameterToStorageMap, &it, this](VectorOfArraysT& stressXX,
                                                                       VectorOfArraysT& stressYY,
                                                                       VectorOfArraysT& stressZZ,
                                                                       VectorOfArraysT& stressXY,
                                                                       VectorOfArraysT& stressYZ,
                                                                       VectorOfArraysT& stressXZ,
                                                                       bool readNucleation) {
      // return pointer to first element
      auto getRawData = [](VectorOfArraysT& vectorOfArrays) {
        return vectorOfArrays.data()->data();
      };
      // fault can be either initialized by traction or by cartesian stress
      // this method reads either the nucleation stress or the initial stress
      std::vector<std::string> identifiers = this->stressIdentifiers(readNucleation);
      bool isFaultParameterizedByTraction = identifiers.size() == 3;
      if (isFaultParameterizedByTraction) {
        // only read traction in normal, strike and dip direction
        parameterToStorageMap.insert({identifiers[0], getRawData(stressXX)});
        parameterToStorageMap.insert({identifiers[1], getRawData(stressXY)});
        parameterToStorageMap.insert({identifiers[2], getRawData(stressXZ)});
        // set the rest to zero
        for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
          for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
            stressYY[ltsFace][pointIndex] = 0.0;
            stressZZ[ltsFace][pointIndex] = 0.0;
            stressYZ[ltsFace][pointIndex] = 0.0;
          }
        }
      } else { // read all stress components from the parameter file
        parameterToStorageMap.insert({identifiers[0], getRawData(stressXX)});
        parameterToStorageMap.insert({identifiers[1], getRawData(stressYY)});
        parameterToStorageMap.insert({identifiers[2], getRawData(stressZZ)});
        parameterToStorageMap.insert({identifiers[3], getRawData(stressXY)});
        parameterToStorageMap.insert({identifiers[4], getRawData(stressYZ)});
        parameterToStorageMap.insert({identifiers[5], getRawData(stressXZ)});
      }
      return isFaultParameterizedByTraction;
    };

    VectorOfArraysT initialStressXX(it->getNumberOfCells());
    VectorOfArraysT initialStressYY(it->getNumberOfCells());
    VectorOfArraysT initialStressZZ(it->getNumberOfCells());
    VectorOfArraysT initialStressXY(it->getNumberOfCells());
    VectorOfArraysT initialStressXZ(it->getNumberOfCells());
    VectorOfArraysT initialStressYZ(it->getNumberOfCells());
    const bool initialStressParameterizedByTraction = addStressesToStorageMap(initialStressXX,
                                                                              initialStressYY,
                                                                              initialStressZZ,
                                                                              initialStressXY,
                                                                              initialStressYZ,
                                                                              initialStressXZ,
                                                                              false);
    VectorOfArraysT nucleationStressXX(it->getNumberOfCells());
    VectorOfArraysT nucleationStressYY(it->getNumberOfCells());
    VectorOfArraysT nucleationStressZZ(it->getNumberOfCells());
    VectorOfArraysT nucleationStressXY(it->getNumberOfCells());
    VectorOfArraysT nucleationStressXZ(it->getNumberOfCells());
    VectorOfArraysT nucleationStressYZ(it->getNumberOfCells());
    const bool nucleationStressParameterizedByTraction = addStressesToStorageMap(nucleationStressXX,
                                                                                 nucleationStressYY,
                                                                                 nucleationStressZZ,
                                                                                 nucleationStressXY,
                                                                                 nucleationStressYZ,
                                                                                 nucleationStressXZ,
                                                                                 true);

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, it);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, it);
    queryModel(faultParameterDB, faceIDs);

    // rotate initial stress to fault coordinate system
    if (initialStressParameterizedByTraction) {
      rotateTractionToCartesianStress(dynRup,
                                      it,
                                      initialStressXX,
                                      initialStressYY,
                                      initialStressZZ,
                                      initialStressXY,
                                      initialStressYZ,
                                      initialStressXZ);
    }
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->initialStressInFaultCS);
    rotateStressToFaultCS(dynRup,
                          it,
                          initialStressInFaultCS,
                          initialStressXX,
                          initialStressYY,
                          initialStressZZ,
                          initialStressXY,
                          initialStressYZ,
                          initialStressXZ);
    // rotate nucleation stress to fault coordinate system
    if (nucleationStressParameterizedByTraction) {
      rotateTractionToCartesianStress(dynRup,
                                      it,
                                      nucleationStressXX,
                                      nucleationStressYY,
                                      nucleationStressZZ,
                                      nucleationStressXY,
                                      nucleationStressYZ,
                                      nucleationStressXZ);
    }
    real(*nucleationStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    rotateStressToFaultCS(dynRup,
                          it,
                          nucleationStressInFaultCS,
                          nucleationStressXX,
                          nucleationStressYY,
                          nucleationStressZZ,
                          nucleationStressXY,
                          nucleationStressYZ,
                          nucleationStressXZ);

    initializeOtherVariables(dynRup, it);
  }
}

std::vector<unsigned> BaseDRInitializer::getFaceIDsInIterator(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it) {
  const auto& drFaceInformation = it->var(dynRup->faceInformation);
  std::vector<unsigned> faceIDs;
  faceIDs.reserve(it->getNumberOfCells());
  // collect all face IDs within this lts leaf
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    faceIDs.push_back(drFaceInformation[ltsFace].meshFace);
  }
  return faceIDs;
}

void BaseDRInitializer::queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                                   std::vector<unsigned> const& faceIDs) {
  // create a query and evaluate the model
  seissol::initializers::FaultGPGenerator queryGen(seissol::SeisSol::main.meshReader(), faceIDs);
  faultParameterDB.evaluateModel(drParameters->faultFileName, queryGen);
}

void BaseDRInitializer::rotateTractionToCartesianStress(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressXX,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressYY,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressZZ,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressXY,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressYZ,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressXZ) {
  // create rotation kernel
  real faultTractionToCartesianMatrixValues[init::stressRotationMatrix::size()];
  auto faultTractionToCartesianMatrixView =
      init::stressRotationMatrix::view::create(faultTractionToCartesianMatrixValues);
  dynamicRupture::kernel::rotateStress faultTractionToCartesianRotationKernel;
  faultTractionToCartesianRotationKernel.stressRotationMatrix =
      faultTractionToCartesianMatrixValues;

  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);

    // if we read the traction in strike, dip and normal direction, we first transform it to stress
    // in cartesian coordinates
    VrtxCoords strike{};
    VrtxCoords dip{};
    misc::computeStrikeAndDipVectors(fault.normal, strike, dip);
    seissol::transformations::symmetricTensor2RotationMatrix(
        fault.normal, strike, dip, faultTractionToCartesianMatrixView, 0, 0);

    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      const real initialTraction[init::initialStress::size()] = {stressXX[ltsFace][pointIndex],
                                                                 stressYY[ltsFace][pointIndex],
                                                                 stressZZ[ltsFace][pointIndex],
                                                                 stressXY[ltsFace][pointIndex],
                                                                 stressYZ[ltsFace][pointIndex],
                                                                 stressXZ[ltsFace][pointIndex]};
      assert(std::abs(initialTraction[1]) < 1e-15);
      assert(std::abs(initialTraction[2]) < 1e-15);
      assert(std::abs(initialTraction[4]) < 1e-15);

      real cartesianStress[init::initialStress::size()]{};
      faultTractionToCartesianRotationKernel.initialStress = initialTraction;
      faultTractionToCartesianRotationKernel.rotatedStress = cartesianStress;
      faultTractionToCartesianRotationKernel.execute();
      stressXX[ltsFace][pointIndex] = cartesianStress[0];
      stressYY[ltsFace][pointIndex] = cartesianStress[1];
      stressZZ[ltsFace][pointIndex] = cartesianStress[2];
      stressXY[ltsFace][pointIndex] = cartesianStress[3];
      stressYZ[ltsFace][pointIndex] = cartesianStress[4];
      stressXZ[ltsFace][pointIndex] = cartesianStress[5];
    }
  }
}

void BaseDRInitializer::rotateStressToFaultCS(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    real (*stressInFaultCS)[misc::numPaddedPoints][6],
    std::vector<std::array<real, misc::numPaddedPoints>> const& stressXX,
    std::vector<std::array<real, misc::numPaddedPoints>> const& stressYY,
    std::vector<std::array<real, misc::numPaddedPoints>> const& stressZZ,
    std::vector<std::array<real, misc::numPaddedPoints>> const& stressXY,
    std::vector<std::array<real, misc::numPaddedPoints>> const& stressYZ,
    std::vector<std::array<real, misc::numPaddedPoints>> const& stressXZ) {
  // create rotation kernel
  real cartesianToFaultCSMatrixValues[init::stressRotationMatrix::size()];
  auto cartesianToFaultCSMatrixView =
      init::stressRotationMatrix::view::create(cartesianToFaultCSMatrixValues);
  dynamicRupture::kernel::rotateStress cartesianToFaultCSRotationKernel;
  cartesianToFaultCSRotationKernel.stressRotationMatrix = cartesianToFaultCSMatrixValues;

  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    constexpr unsigned int numberOfStressComponents = 6;
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);

    // now rotate the stress in cartesian coordinates to the element aligned coordinate system.
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, cartesianToFaultCSMatrixView, 0, 0);

    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      const real initialStress[init::initialStress::size()] = {stressXX[ltsFace][pointIndex],
                                                               stressYY[ltsFace][pointIndex],
                                                               stressZZ[ltsFace][pointIndex],
                                                               stressXY[ltsFace][pointIndex],
                                                               stressYZ[ltsFace][pointIndex],
                                                               stressXZ[ltsFace][pointIndex]};
      real rotatedStress[init::initialStress::size()]{};
      cartesianToFaultCSRotationKernel.initialStress = initialStress;
      cartesianToFaultCSRotationKernel.rotatedStress = rotatedStress;
      cartesianToFaultCSRotationKernel.execute();
      for (unsigned int stressIndex = 0; stressIndex < numberOfStressComponents; ++stressIndex) {
        stressInFaultCS[ltsFace][pointIndex][stressIndex] = rotatedStress[stressIndex];
      }
    }
  }
}

void BaseDRInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing for base friction law
}

void BaseDRInitializer::initializeOtherVariables(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // initialize rupture front flag
  bool(*ruptureTimePending)[misc::numPaddedPoints] = it->var(dynRup->ruptureTimePending);
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      ruptureTimePending[ltsFace][pointIndex] = drParameters->isRfOutputOn;
    }
  }

  // initialize all other variables to zero
  real(*peakSlipRate)[misc::numPaddedPoints] = it->var(dynRup->peakSlipRate);
  real(*ruptureTime)[misc::numPaddedPoints] = it->var(dynRup->ruptureTime);
  real(*dynStressTime)[misc::numPaddedPoints] = it->var(dynRup->dynStressTime);
  real(*accumulatedSlipMagnitude)[misc::numPaddedPoints] =
      it->var(dynRup->accumulatedSlipMagnitude);
  real(*slip1)[misc::numPaddedPoints] = it->var(dynRup->slip1);
  real(*slip2)[misc::numPaddedPoints] = it->var(dynRup->slip2);
  real(*slipRateMagnitude)[misc::numPaddedPoints] = it->var(dynRup->slipRateMagnitude);
  real(*traction1)[misc::numPaddedPoints] = it->var(dynRup->traction1);
  real(*traction2)[misc::numPaddedPoints] = it->var(dynRup->traction2);

  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
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

bool BaseDRInitializer::faultProvides(std::string&& parameter) {
  return seissol::initializers::FaultParameterDB::faultProvides(parameter,
                                                                drParameters->faultFileName);
}

std::vector<std::string> BaseDRInitializer::stressIdentifiers(bool readNucleation) {
  std::vector<std::string> tractionNames;
  std::vector<std::string> cartesianNames;
  if (readNucleation) {
    tractionNames = {"Tnuc_n", "Tnuc_s", "Tnuc_d"};
    cartesianNames = {"nuc_xx", "nuc_yy", "nuc_zz", "nuc_xy", "nuc_yz", "nuc_xz"};
  } else {
    tractionNames = {"T_n", "T_s", "T_d"};
    cartesianNames = {"s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"};
  }

  bool allTractionParametersSupplied = true;
  bool allCartesianParametersSupplied = true;
  bool anyTractionParametersSupplied = false;
  bool anyCartesianParametersSupplied = false;
  for (size_t i = 0; i < 3; i++) {
    auto b = seissol::initializers::FaultParameterDB::faultProvides(tractionNames[i],
                                                                    drParameters->faultFileName);
    allTractionParametersSupplied &= b;
    anyTractionParametersSupplied |= b;
  }
  for (size_t i = 0; i < 6; i++) {
    const auto b = seissol::initializers::FaultParameterDB::faultProvides(
        cartesianNames[i], drParameters->faultFileName);
    allCartesianParametersSupplied &= b;
    anyCartesianParametersSupplied |= b;
  }

  if (allCartesianParametersSupplied && !anyTractionParametersSupplied) {
    return cartesianNames;
  } else if (allTractionParametersSupplied && !anyCartesianParametersSupplied) {
    return tractionNames;
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

} // namespace seissol::dr::initializers
