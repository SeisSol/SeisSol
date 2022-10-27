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
    auto addStressesToStorageMap = [&parameterToStorageMap, &it, this](StressTensor& initialStress,
                                                                       bool readNucleation) {
      // return pointer to first element
      auto getRawData = [](StressTensor::VectorOfArrays_t& vectorOfArrays) {
        return vectorOfArrays.data()->data();
      };
      // fault can be either initialized by traction or by cartesian stress
      // this method reads either the nucleation stress or the initial stress
      std::vector<std::string> identifiers;
      Parametrization parametrization;
      std::tie(identifiers, parametrization) = this->stressIdentifiers(readNucleation);
      bool isFaultParameterizedByTraction = parametrization == Parametrization::Traction;
      if (isFaultParameterizedByTraction) {
        // only read traction in normal, strike and dip direction
        parameterToStorageMap.insert({identifiers[0], getRawData(initialStress.xx)});
        parameterToStorageMap.insert({identifiers[1], getRawData(initialStress.xy)});
        parameterToStorageMap.insert({identifiers[2], getRawData(initialStress.xz)});
        // set the rest to zero
        for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
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
      return isFaultParameterizedByTraction;
    };

    StressTensor initialStress(it->getNumberOfCells());
    const bool initialStressParameterizedByTraction = addStressesToStorageMap(initialStress, false);
    StressTensor nucleationStress(it->getNumberOfCells());
    const bool nucleationStressParameterizedByTraction =
        addStressesToStorageMap(nucleationStress, true);

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
      rotateTractionToCartesianStress(dynRup, it, initialStress);
    }

    auto* initialStressInFaultCS = it->var(dynRup->initialStressInFaultCS);
    rotateStressToFaultCS(dynRup, it, initialStressInFaultCS, initialStress);
    // rotate nucleation stress to fault coordinate system
    if (nucleationStressParameterizedByTraction) {
      rotateTractionToCartesianStress(dynRup, it, nucleationStress);
    }
    real(*nucleationStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    rotateStressToFaultCS(dynRup, it, nucleationStressInFaultCS, nucleationStress);

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
    StressTensor& stress) {
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

    using namespace dr::misc::quantity_indices;
    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
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
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    real (*stressInFaultCS)[misc::numPaddedPoints][6],
    StressTensor const& stress) {
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

} // namespace seissol::dr::initializers
