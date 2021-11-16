#include "BaseDRInitializer.h"

#include <Eigen/Dense>
#include "Initializer/ParameterDB.h"
#include "SeisSol.h"
#include "Numerical_aux/Quadrature.h"
#include "Model/common.hpp"
#include "generated_code/kernel.h"

namespace seissol::dr::initializers {
void BaseDRInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                        seissol::initializers::LTSTree* dynRupTree,
                                        seissol::Interoperability* e_interoperability) {
  if (drParameters.faultFileName == "") {
    // Todo find better way to enable/disable dynamic rupture
    logInfo() << "No fault filename set. Assume, we don't compute dynamic rupture";
    return;
  }
  seissol::initializers::FaultParameterDB faultParameterDB;
  dynRup->isFaultParameterizedByTraction =
      faultParameterDB.faultParameterizedByTraction(drParameters.faultFileName);
  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    // read initial stress and nucleation stress
    using VectorOfArrays = std::vector<std::array<real, numPaddedPoints>>;

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
        // only read traction in normal strike and slip direction
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

    VectorOfArrays iniXX(it->getNumberOfCells());
    VectorOfArrays iniYY(it->getNumberOfCells());
    VectorOfArrays iniZZ(it->getNumberOfCells());
    VectorOfArrays iniXY(it->getNumberOfCells());
    VectorOfArrays iniXZ(it->getNumberOfCells());
    VectorOfArrays iniYZ(it->getNumberOfCells());
    addStressesToStorageMap(iniXX, iniYY, iniZZ, iniXY, iniYZ, iniXZ, false);
    VectorOfArrays nucXX(it->getNumberOfCells());
    VectorOfArrays nucYY(it->getNumberOfCells());
    VectorOfArrays nucZZ(it->getNumberOfCells());
    VectorOfArrays nucXY(it->getNumberOfCells());
    VectorOfArrays nucXZ(it->getNumberOfCells());
    VectorOfArrays nucYZ(it->getNumberOfCells());
    addStressesToStorageMap(nucXX, nucYY, nucZZ, nucXY, nucYZ, nucXZ, true);

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, it);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, it);
    queryModel(faultParameterDB, faceIDs);

    // rotate initial stress to fault coordinate system
    real(*initialStressInFaultCS)[numPaddedPoints][6] = it->var(dynRup->initialStressInFaultCS);
    rotateStressToFaultCS(
        dynRup, it, initialStressInFaultCS, iniXX, iniYY, iniZZ, iniXY, iniYZ, iniXZ);
    // rotate nucleation stress to fault coordinate system
    real(*nucleationStressInFaultCS)[numPaddedPoints][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    rotateStressToFaultCS(
        dynRup, it, nucleationStressInFaultCS, nucXX, nucYY, nucZZ, nucXY, nucYZ, nucXZ);

    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      e_interoperability->copyFrictionOutputToFortranInitialStressInFaultCS(
          ltsFace, meshFace, initialStressInFaultCS, iniXX, iniYY, iniZZ, iniXY, iniYZ, iniXZ);
    }

    // initialize rupture front flag
    bool(*ruptureFront)[numPaddedPoints] = it->var(dynRup->ruptureFront);
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        ruptureFront[ltsFace][pointIndex] = drParameters.isRfOutputOn;
      }
    }

    // initialize all other variables to zero
    real(*peakSlipRate)[numPaddedPoints] = it->var(dynRup->peakSlipRate);
    real(*ruptureTime)[numPaddedPoints] = it->var(dynRup->ruptureTime);
    real(*dynStressTime)[numPaddedPoints] = it->var(dynRup->dynStressTime);
    real(*slip)[numPaddedPoints] = it->var(dynRup->slip);
    real(*slipDip)[numPaddedPoints] = it->var(dynRup->slipDip);
    real(*slipStrike)[numPaddedPoints] = it->var(dynRup->slipStrike);
    real(*slipRateMagnitude)[numPaddedPoints] = it->var(dynRup->slipRateMagnitude);
    real(*tractionXY)[numPaddedPoints] = it->var(dynRup->tractionXY);
    real(*tractionXZ)[numPaddedPoints] = it->var(dynRup->tractionXZ);

    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        peakSlipRate[ltsFace][pointIndex] = 0;
        ruptureTime[ltsFace][pointIndex] = 0;
        dynStressTime[ltsFace][pointIndex] = 0;
        slip[ltsFace][pointIndex] = 0;
        slipDip[ltsFace][pointIndex] = 0;
        slipStrike[ltsFace][pointIndex] = 0;
        slipRateMagnitude[ltsFace][pointIndex] = 0;
        tractionXY[ltsFace][pointIndex] = 0;
        tractionXZ[ltsFace][pointIndex] = 0;
      }
    }
    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      e_interoperability->copyFrictionOutputToFortranGeneral(ltsFace,
                                                             meshFace,
                                                             slip,
                                                             slipStrike,
                                                             slipDip,
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
  double boundaryGaussPoints[numPaddedPoints][2] = {0};
  double weights[numPaddedPoints] = {0};
  quadrature::TriangleQuadrature(boundaryGaussPoints, weights, CONVERGENCE_ORDER + 1);
  seissol::initializers::FaultGPGenerator queryGen(
      seissol::SeisSol::main.meshReader(),
      reinterpret_cast<double(*)[2]>(boundaryGaussPoints),
      numPaddedPoints,
      faceIDs);
  faultParameterDB.evaluateModel(drParameters.faultFileName, queryGen);
}

void BaseDRInitializer::rotateStressToFaultCS(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    real (*stressInFaultCS)[numPaddedPoints][6],
    std::vector<std::array<real, numPaddedPoints>>& stressXX,
    std::vector<std::array<real, numPaddedPoints>>& stressYY,
    std::vector<std::array<real, numPaddedPoints>>& stressZZ,
    std::vector<std::array<real, numPaddedPoints>>& stressXY,
    std::vector<std::array<real, numPaddedPoints>>& stressYZ,
    std::vector<std::array<real, numPaddedPoints>>& stressXZ) {
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

    for (unsigned int pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
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
