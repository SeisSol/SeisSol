#ifndef SEISSOL_BASEDRINITIALIZER_H
#define SEISSOL_BASEDRINITIALIZER_H

#include <Solver/Interoperability.h>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/FrictionLaws/FrictionLaws.h"
#include "Initializer/InputAux.hpp"
#include "Initializer/ParameterDB.h"
#include "DynamicRupture/Parameters.h"

namespace seissol::dr::initializers {
class BaseDRInitializer; // general parameters initialized that are required by all friction laws
}

/*
 * initial values are obtained
 * from m_Params (direct in read from parameter.par file)
 * from easi, accessed through FaultParameterDB
 */
class seissol::dr::initializers::BaseDRInitializer {
  protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  static constexpr int numPaddedPoints = init::QInterpolated::Stop[0];

  DRParameters& drParameters;

  public:
  BaseDRInitializer(DRParameters& drParameters) : drParameters(drParameters){};
  virtual ~BaseDRInitializer() {}

  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability);

  protected:
  virtual void addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                                       seissol::initializers::DynamicRupture* dynRup,
                                       seissol::initializers::LTSInternalNode::leaf_iterator& it);

  private:
  std::vector<unsigned>
      getFaceIDsInIterator(seissol::initializers::LTSInternalNode::leaf_iterator& it,
                           seissol::initializers::DynamicRupture* dynRup);

  void queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                  std::vector<unsigned> faceIDs);

  void rotateStressToFaultCS(seissol::initializers::LTSTree::leaf_iterator& it,
                             seissol::initializers::DynamicRupture* dynRup,
                             real (*initialStressInFaultCS)[numPaddedPoints][6],
                             real (*iniBulkXX)[numPaddedPoints],
                             real (*iniBulkYY)[numPaddedPoints],
                             real (*iniBulkZZ)[numPaddedPoints],
                             real (*iniShearXY)[numPaddedPoints],
                             real (*iniShearYZ)[numPaddedPoints],
                             real (*iniShearXZ)[numPaddedPoints]);
};

#endif // SEISSOL_BASEDRINITIALIZER_H
