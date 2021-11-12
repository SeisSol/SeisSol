#include "BaseDRInitializer.h"

#include <Eigen/Dense>
#include "Initializer/ParameterDB.h"
#include "SeisSol.h"
#include "Numerical_aux/Quadrature.h"
#include "Model/common.hpp"

namespace seissol::dr::initializers {
void BaseDRInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                        seissol::initializers::LTSTree* dynRupTree,
                                        seissol::Interoperability* e_interoperability) {
  if (!drParameters.isGpWiseInitialization) {
    logError() << "Dynamic Rupture with cell average not supported any more";
  }

  seissol::initializers::FaultParameterDB faultParameterDB;
  dynRup->isFaultParameterizedByTraction =
      faultParameterDB.faultParameterizedByTraction(drParameters.faultFileName);
  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::map<std::string, real*> parameterToStorageMap;

    // read initial stress
    real(*iniBulkXX)[numPaddedPoints] = it->var(dynRup->iniBulkXX);
    real(*iniBulkYY)[numPaddedPoints] = it->var(dynRup->iniBulkYY);
    real(*iniBulkZZ)[numPaddedPoints] = it->var(dynRup->iniBulkZZ);
    real(*iniShearXY)[numPaddedPoints] = it->var(dynRup->iniShearXY);
    real(*iniShearXZ)[numPaddedPoints] = it->var(dynRup->iniShearXZ);
    real(*iniShearYZ)[numPaddedPoints] = it->var(dynRup->iniShearYZ);
    if (dynRup->isFaultParameterizedByTraction) {
      // only read traction in normal strike and slip direction
      parameterToStorageMap.insert({"T_n", (real*)iniBulkXX});
      parameterToStorageMap.insert({"T_s", (real*)iniShearXY});
      parameterToStorageMap.insert({"T_d", (real*)iniShearXZ});
      // set the rest to zero
      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
          iniBulkYY[ltsFace][pointIndex] = 0.0;
          iniBulkZZ[ltsFace][pointIndex] = 0.0;
          iniShearYZ[ltsFace][pointIndex] = 0.0;
        }
      }
    } else {
      // read all stress components from the parameter file
      parameterToStorageMap.insert({"s_xx", (real*)iniBulkXX});
      parameterToStorageMap.insert({"s_yy", (real*)iniBulkYY});
      parameterToStorageMap.insert({"s_zz", (real*)iniBulkZZ});
      parameterToStorageMap.insert({"s_xy", (real*)iniShearXY});
      parameterToStorageMap.insert({"s_yz", (real*)iniShearYZ});
      parameterToStorageMap.insert({"s_xz", (real*)iniShearXZ});
    }

    // do the same for nucleation stress
    real(*nucXX)[numPaddedPoints] = new real[it->getNumberOfCells()][numPaddedPoints];
    real(*nucYY)[numPaddedPoints] = new real[it->getNumberOfCells()][numPaddedPoints];
    real(*nucZZ)[numPaddedPoints] = new real[it->getNumberOfCells()][numPaddedPoints];
    real(*nucXY)[numPaddedPoints] = new real[it->getNumberOfCells()][numPaddedPoints];
    real(*nucXZ)[numPaddedPoints] = new real[it->getNumberOfCells()][numPaddedPoints];
    real(*nucYZ)[numPaddedPoints] = new real[it->getNumberOfCells()][numPaddedPoints];
    if (dynRup->isFaultParameterizedByTraction) {
      parameterToStorageMap.insert({"Tnuc_n", (real*)nucXX});
      parameterToStorageMap.insert({"Tnuc_s", (real*)nucXY});
      parameterToStorageMap.insert({"Tnuc_d", (real*)nucXZ});
      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
          nucYY[ltsFace][pointIndex] = 0.0;
          nucZZ[ltsFace][pointIndex] = 0.0;
          nucYZ[ltsFace][pointIndex] = 0.0;
        }
      }
    } else {
      parameterToStorageMap.insert({"nuc_xx", (real*)nucXX});
      parameterToStorageMap.insert({"nuc_yy", (real*)nucYY});
      parameterToStorageMap.insert({"nuc_zz", (real*)nucZZ});
      parameterToStorageMap.insert({"nuc_xy", (real*)nucXY});
      parameterToStorageMap.insert({"nuc_yz", (real*)nucYZ});
      parameterToStorageMap.insert({"nuc_xz", (real*)nucXZ});
    }

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
    rotateStressToFaultCS(dynRup,
                          it,
                          initialStressInFaultCS,
                          iniBulkXX,
                          iniBulkYY,
                          iniBulkZZ,
                          iniShearXY,
                          iniShearYZ,
                          iniShearXZ);
    // rotate nucleation stress to fault coordinate system
    real(*nucleationStressInFaultCS)[numPaddedPoints][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    rotateStressToFaultCS(
        dynRup, it, nucleationStressInFaultCS, nucXX, nucYY, nucZZ, nucXY, nucYZ, nucXZ);
    delete[] nucXX;
    delete[] nucYY;
    delete[] nucZZ;
    delete[] nucXY;
    delete[] nucYZ;
    delete[] nucXZ;

    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      e_interoperability->copyFrictionOutputToFortranInitialStressInFaultCS(ltsFace,
                                                                            meshFace,
                                                                            initialStressInFaultCS,
                                                                            iniBulkXX,
                                                                            iniBulkYY,
                                                                            iniBulkZZ,
                                                                            iniShearXY,
                                                                            iniShearYZ,
                                                                            iniShearXZ);
    }

    // initialize rupture front flag
    bool(*ruptureFront)[numPaddedPoints] = it->var(dynRup->ruptureFront);
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int j = 0; j < numPaddedPoints; ++j) {
        ruptureFront[ltsFace][j] = drParameters.isRfOutputOn;
      }
    }

    // initialize all other variables to zero
    real(*peakSlipRate)[numPaddedPoints] = it->var(dynRup->peakSlipRate);
    real(*ruptureTime)[numPaddedPoints] = it->var(dynRup->ruptureTime);
    real(*slip)[numPaddedPoints] = it->var(dynRup->slip);
    real(*slipDip)[numPaddedPoints] = it->var(dynRup->slipDip);
    real(*slipStrike)[numPaddedPoints] = it->var(dynRup->slipStrike);
    real(*slipRateMagnitude)[numPaddedPoints] = it->var(dynRup->slipRateMagnitude);
    real(*tractionXY)[numPaddedPoints] = it->var(dynRup->tractionXY);
    real(*tractionXZ)[numPaddedPoints] = it->var(dynRup->tractionXZ);

    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int j = 0; j < numPaddedPoints; ++j) {
        peakSlipRate[ltsFace][j] = 0;
        ruptureTime[ltsFace][j] = 0;
        slip[ltsFace][j] = 0;
        slipDip[ltsFace][j] = 0;
        slipStrike[ltsFace][j] = 0;
        slipRateMagnitude[ltsFace][j] = 0;
        tractionXY[ltsFace][j] = 0;
        tractionXZ[ltsFace][j] = 0;
      }
    }
    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      e_interoperability->copyFrictionOutputToFortran(ltsFace,
                                                      meshFace,
                                                      slip,
                                                      slipStrike,
                                                      slipDip,
                                                      ruptureTime,
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

void BaseDRInitializer::rotateStressToFaultCS(seissol::initializers::DynamicRupture* dynRup,
                                              seissol::initializers::LTSTree::leaf_iterator& it,
                                              real (*stressInFaultCS)[numPaddedPoints][6],
                                              real (*stressXX)[numPaddedPoints],
                                              real (*stressYY)[numPaddedPoints],
                                              real (*stressZZ)[numPaddedPoints],
                                              real (*stressXY)[numPaddedPoints],
                                              real (*stressYZ)[numPaddedPoints],
                                              real (*stressXZ)[numPaddedPoints]) {
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);
    Eigen::Matrix<double, 6, 6> rotationMatrix;
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, rotationMatrix, 0, 0);

    for (unsigned int j = 0; j < numPaddedPoints; ++j) {
      Eigen::Vector<double, 6> stress;
      stress << stressXX[ltsFace][j], stressYY[ltsFace][j], stressZZ[ltsFace][j],
          stressXY[ltsFace][j], stressYZ[ltsFace][j], stressXZ[ltsFace][j];
      Eigen::Vector<double, 6> rotatedStress = rotationMatrix * stress;
      for (int k = 0; k < 6; ++k) {
        stressInFaultCS[ltsFace][j][k] = rotatedStress(k);
      }
    }
  }
}

void BaseDRInitializer::addAdditionalParameters(
    std::map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing for base friction law
}

} // namespace seissol::dr::initializers
