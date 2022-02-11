#include "BaseDRInitializer.h"

#include "Initializer/ParameterDB.h"
#include "Model/common.hpp"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <Eigen/Dense>

namespace seissol::dr::initializers {
void BaseDRInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                        seissol::initializers::LTSTree* dynRupTree,
                                        seissol::Interoperability* eInteroperability) {
  logInfo() << "Initialize Fault, using a quadrature rule with "
            << misc::numberOfBoundaryGaussPoints << " points.";
  seissol::initializers::FaultParameterDB faultParameterDB;
  dynRup->isFaultParameterizedByTraction =
      faultParameterDB.faultParameterizedByTraction(drParameters.faultFileName);
  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    // read initial stress and nucleation stress
    using VectorOfArrays = std::vector<std::array<real, misc::numPaddedPoints>>;

    auto addStressesToStorageMap = [&dynRup, &parameterToStorageMap, &it](VectorOfArrays& stressXX,
                                                                          VectorOfArrays& stressYY,
                                                                          VectorOfArrays& stressZZ,
                                                                          VectorOfArrays& stressXY,
                                                                          VectorOfArrays& stressYZ,
                                                                          VectorOfArrays& stressXZ,
                                                                          bool readNucleation) {
      // return pointer to first element
      auto getRawData = [](VectorOfArrays& vectorOfArrays) {
        return vectorOfArrays.data()->data();
      };
      // fault can be either initialized by traction or by cartesian stress
      // this method reads either the nucleation stress or the initial stress
      std::array<std::string, 3> tractionNames;
      std::array<std::string, 6> cartesianNames;
      if (readNucleation) {
        tractionNames = {"Tnuc_n", "Tnuc_d", "Tnuc_s"};
        cartesianNames = {"nuc_xx", "nuc_yy", "nuc_zz", "nuc_xy", "nuc_yz", "nuc_xz"};
      } else {
        tractionNames = {"T_n", "T_s", "T_d"};
        cartesianNames = {"s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"};
      }
      if (dynRup->isFaultParameterizedByTraction) {
        // only read traction in normal, strike and dip direction
        parameterToStorageMap.insert({tractionNames[0], getRawData(stressXX)});
        parameterToStorageMap.insert({tractionNames[1], getRawData(stressXY)});
        parameterToStorageMap.insert({tractionNames[2], getRawData(stressXZ)});
        // set the rest to zero
        for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
          for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
            stressYY[ltsFace][pointIndex] = 0.0;
            stressZZ[ltsFace][pointIndex] = 0.0;
            stressYZ[ltsFace][pointIndex] = 0.0;
          }
        }
      } else {
        // read all stress components from the parameter file
        parameterToStorageMap.insert({cartesianNames[0], getRawData(stressXX)});
        parameterToStorageMap.insert({cartesianNames[1], getRawData(stressYY)});
        parameterToStorageMap.insert({cartesianNames[2], getRawData(stressZZ)});
        parameterToStorageMap.insert({cartesianNames[3], getRawData(stressXY)});
        parameterToStorageMap.insert({cartesianNames[4], getRawData(stressYZ)});
        parameterToStorageMap.insert({cartesianNames[5], getRawData(stressXZ)});
      }
    };

    VectorOfArrays initialStressXX(it->getNumberOfCells());
    VectorOfArrays initialStressYY(it->getNumberOfCells());
    VectorOfArrays initialStressZZ(it->getNumberOfCells());
    VectorOfArrays initialStressXY(it->getNumberOfCells());
    VectorOfArrays initialStressXZ(it->getNumberOfCells());
    VectorOfArrays initialStressYZ(it->getNumberOfCells());
    addStressesToStorageMap(initialStressXX,
                            initialStressYY,
                            initialStressZZ,
                            initialStressXY,
                            initialStressYZ,
                            initialStressXZ,
                            false);
    VectorOfArrays nucleationStressXX(it->getNumberOfCells());
    VectorOfArrays nucleationStressYY(it->getNumberOfCells());
    VectorOfArrays nucleationStressZZ(it->getNumberOfCells());
    VectorOfArrays nucleationStressXY(it->getNumberOfCells());
    VectorOfArrays nucleationStressXZ(it->getNumberOfCells());
    VectorOfArrays nucleationStressYZ(it->getNumberOfCells());
    addStressesToStorageMap(nucleationStressXX,
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

    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      eInteroperability->copyFrictionOutputToFortranInitialStressInFaultCS(ltsFace,
                                                                           meshFace,
                                                                           initialStressInFaultCS,
                                                                           initialStressXX,
                                                                           initialStressYY,
                                                                           initialStressZZ,
                                                                           initialStressXY,
                                                                           initialStressYZ,
                                                                           initialStressXZ);
    }

    // initialize rupture front flag
    bool(*ruptureTimePending)[misc::numPaddedPoints] = it->var(dynRup->ruptureTimePending);
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        ruptureTimePending[ltsFace][pointIndex] = drParameters.isRfOutputOn;
      }
    }

    // initialize all other variables to zero
    real(*peakSlipRate)[misc::numPaddedPoints] = it->var(dynRup->peakSlipRate);
    real(*ruptureTime)[misc::numPaddedPoints] = it->var(dynRup->ruptureTime);
    real(*dynStressTime)[misc::numPaddedPoints] = it->var(dynRup->dynStressTime);
    real(*accumulatedSlipMagnitude)[misc::numPaddedPoints] =
        it->var(dynRup->accumulatedSlipMagnitude);
    real(*slip1)[misc::numPaddedPoints] = it->var(dynRup->slip2);
    real(*slip2)[misc::numPaddedPoints] = it->var(dynRup->slip1);
    real(*slipRateMagnitude)[misc::numPaddedPoints] = it->var(dynRup->slipRateMagnitude);
    real(*tractionXY)[misc::numPaddedPoints] = it->var(dynRup->tractionXY);
    real(*tractionXZ)[misc::numPaddedPoints] = it->var(dynRup->tractionXZ);

    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        peakSlipRate[ltsFace][pointIndex] = 0;
        ruptureTime[ltsFace][pointIndex] = 0;
        dynStressTime[ltsFace][pointIndex] = 0;
        accumulatedSlipMagnitude[ltsFace][pointIndex] = 0;
        slip1[ltsFace][pointIndex] = 0;
        slip2[ltsFace][pointIndex] = 0;
        slipRateMagnitude[ltsFace][pointIndex] = 0;
        tractionXY[ltsFace][pointIndex] = 0;
        tractionXZ[ltsFace][pointIndex] = 0;
      }
    }
    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      eInteroperability->copyFrictionOutputToFortranGeneral(ltsFace,
                                                            meshFace,
                                                            accumulatedSlipMagnitude,
                                                            slip2,
                                                            slip1,
                                                            ruptureTime,
                                                            dynStressTime,
                                                            peakSlipRate,
                                                            tractionXY,
                                                            tractionXZ);
    }
  }
}

std::vector<unsigned>
    BaseDRInitializer::getFaceIDsInIterator(seissol::initializers::DynamicRupture* dynRup,
                                            seissol::initializers::LTSTree::leaf_iterator& it) {
  const auto& drFaceInformation = it->var(dynRup->faceInformation);
  std::vector<unsigned> faceIDs;
  // collect all face IDs within this lts leaf
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    faceIDs.push_back(drFaceInformation[ltsFace].meshFace);
  }
  return faceIDs;
}

void BaseDRInitializer::queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                                   std::vector<unsigned> faceIDs) {
  // create a query and evaluate the model
  seissol::initializers::FaultGPGenerator queryGen(seissol::SeisSol::main.meshReader(), faceIDs);
  faultParameterDB.evaluateModel(drParameters.faultFileName, queryGen);
}

void BaseDRInitializer::rotateStressToFaultCS(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    real (*stressInFaultCS)[misc::numPaddedPoints][6],
    std::vector<std::array<real, misc::numPaddedPoints>>& stressXX,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressYY,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressZZ,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressXY,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressYZ,
    std::vector<std::array<real, misc::numPaddedPoints>>& stressXZ) {
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    constexpr unsigned int numberOfStressComponents = 6;
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);
    real stressRotationMatrixValues[init::stressRotationMatrix::size()];
    auto stressRotationMatrixView =
        init::stressRotationMatrix::view::create(stressRotationMatrixValues);
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, stressRotationMatrixView, 0, 0);

    dynamicRupture::kernel::rotateStressToFaultCS rotationKernel;
    rotationKernel.stressRotationMatrix = stressRotationMatrixValues;

    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      real initialStress[init::initialStress::size()] = {stressXX[ltsFace][pointIndex],
                                                         stressYY[ltsFace][pointIndex],
                                                         stressZZ[ltsFace][pointIndex],
                                                         stressXY[ltsFace][pointIndex],
                                                         stressYZ[ltsFace][pointIndex],
                                                         stressXZ[ltsFace][pointIndex]};
      real rotatedStress[init::initialStress::size()]{};
      rotationKernel.initialStress = initialStress;
      rotationKernel.rotatedStress = rotatedStress;
      rotationKernel.execute();
      for (unsigned int stressIndex = 0; stressIndex < numberOfStressComponents; ++stressIndex) {
        stressInFaultCS[ltsFace][pointIndex][stressIndex] = rotatedStress[stressIndex];
      }
    }
  }
}

void BaseDRInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing for base friction law
}

} // namespace seissol::dr::initializers
