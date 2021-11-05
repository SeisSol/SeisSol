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
 * from std::unordered_map<std::string, double*> faultParameters -> allocated matrices in Fortran
 * with C++ pointer from Easi inread. from interoperabilty functions that access Fortran values and
 * copy them in C++ memory (copy by value) (e.g. e_interoperability.getDynRupParameters())
 */
class seissol::dr::initializers::BaseDRInitializer {
  protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  static constexpr int numPaddedPoints = init::QInterpolated::Stop[0];
  // YAML::Node m_InputParam;

  DRParameters& drParameters;

  public:
  BaseDRInitializer(DRParameters& drParameters) : drParameters(drParameters){};
  virtual ~BaseDRInitializer() {}

  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability);

  std::vector<unsigned>
      getFaceIDsInIterator(seissol::initializers::LTSInternalNode::leaf_iterator& it,
                           seissol::initializers::DynamicRupture* dynRup);

  virtual void addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                                       seissol::initializers::DynamicRupture* dynRup,
                                       seissol::initializers::LTSInternalNode::leaf_iterator& it);

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

  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture* dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          unsigned* ltsFaceToMeshFace);

  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture* dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability& e_interoperability);
};

#endif // SEISSOL_BASEDRINITIALIZER_H
