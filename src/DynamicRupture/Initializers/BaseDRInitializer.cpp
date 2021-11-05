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
    std::map<std::string, double*> parameterToStorageMap;

    real(*cohesion)[numPaddedPoints] = it->var(dynRup->cohesion);
    parameterToStorageMap.insert({"cohesion", (double*)cohesion});

    real(*iniBulkXX)[numPaddedPoints] = it->var(dynRup->iniBulkXX);
    real(*iniBulkYY)[numPaddedPoints] = it->var(dynRup->iniBulkYY);
    real(*iniBulkZZ)[numPaddedPoints] = it->var(dynRup->iniBulkZZ);
    real(*iniShearXY)[numPaddedPoints] = it->var(dynRup->iniShearXY);
    real(*iniShearXZ)[numPaddedPoints] = it->var(dynRup->iniShearXZ);
    real(*iniShearYZ)[numPaddedPoints] = it->var(dynRup->iniShearYZ);
    if (dynRup->isFaultParameterizedByTraction) {
      parameterToStorageMap.insert({"T_n", (double*)iniBulkXX});
      parameterToStorageMap.insert({"T_s", (double*)iniShearXY});
      parameterToStorageMap.insert({"T_d", (double*)iniShearXZ});
      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        for (unsigned pointIndex = 0; pointIndex < init::QInterpolated::Stop[0]; ++pointIndex) {
          iniBulkYY[ltsFace][pointIndex] = 0.0;
          iniBulkZZ[ltsFace][pointIndex] = 0.0;
          iniShearYZ[ltsFace][pointIndex] = 0.0;
        }
      }
    } else {
      parameterToStorageMap.insert({"s_xx", (double*)iniBulkXX});
      parameterToStorageMap.insert({"s_yy", (double*)iniBulkYY});
      parameterToStorageMap.insert({"s_zz", (double*)iniBulkZZ});
      parameterToStorageMap.insert({"s_xy", (double*)iniShearXY});
      parameterToStorageMap.insert({"s_yz", (double*)iniShearYZ});
      parameterToStorageMap.insert({"s_xz", (double*)iniShearXZ});
    }

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, it);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(it, dynRup);
    queryModel(faultParameterDB, faceIDs);

    // rotate initial stress to fault coordinate system
    real(*initialStressInFaultCS)[numPaddedPoints][6] = it->var(dynRup->initialStressInFaultCS);
    rotateStressToFaultCS(it,
                          dynRup,
                          initialStressInFaultCS,
                          iniBulkXX,
                          iniBulkYY,
                          iniBulkZZ,
                          iniShearXY,
                          iniShearYZ,
                          iniShearXZ);

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
    for (int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
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
    BaseDRInitializer::getFaceIDsInIterator(seissol::initializers::LTSTree::leaf_iterator& it,
                                            seissol::initializers::DynamicRupture* dynRup) {
  const auto& meshReader = seissol::SeisSol::main.meshReader();
  const auto& drFaceInformation = it->var(dynRup->faceInformation);
  std::vector<unsigned> faceIDs;
  // collect all face IDs within this lts leaf
  for (int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
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

void BaseDRInitializer::rotateStressToFaultCS(seissol::initializers::LTSTree::leaf_iterator& it,
                                              seissol::initializers::DynamicRupture* dynRup,
                                              real (*initialStressInFaultCS)[52][6],
                                              real (*iniBulkXX)[52],
                                              real (*iniBulkYY)[52],
                                              real (*iniBulkZZ)[52],
                                              real (*iniShearXY)[52],
                                              real (*iniShearYZ)[52],
                                              real (*iniShearXZ)[52]) {
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);
    Eigen::Matrix<double, 6, 6> rotationMatrix;
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, rotationMatrix, 0, 0);

    for (unsigned int j = 0; j < numPaddedPoints; ++j) {
      Eigen::Vector<double, 6> stress;
      stress << iniBulkXX[ltsFace][j], iniBulkYY[ltsFace][j], iniBulkZZ[ltsFace][j],
          iniShearXY[ltsFace][j], iniShearYZ[ltsFace][j], iniShearXZ[ltsFace][j];
      Eigen::Vector<double, 6> rotatedStress = rotationMatrix * stress;
      for (int k = 0; k < 6; ++k) {
        initialStressInFaultCS[ltsFace][j][k] = rotatedStress(k);
      }
    }
  }
}

void BaseDRInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    unsigned* ltsFaceToMeshFace) {}

void BaseDRInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned int* ltsFaceToMeshFace,
    Interoperability& e_interoperability) {}

void BaseDRInitializer::addAdditionalParameters(
    std::map<std::string, double*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing for base friction law
}

} // namespace seissol::dr::initializers