#ifndef SEISSOL_BASEDRINITIALIZER_H
#define SEISSOL_BASEDRINITIALIZER_H

#include <Solver/Interoperability.h>
#include <yaml-cpp/yaml.h>
#include "Initializer/InputAux.hpp"
#include "DynamicRupture/FrictionLaws/FrictionLaws.h"
#include <DynamicRupture/DR_Parameters.h>

namespace seissol::dr::initializers {
  class BaseDRInitializer;            //general parameters initialized that are required by all friction laws
}

/*
 * initial values are obtained
 * from m_Params (direct in read from parameter.par file)
 * from std::unordered_map<std::string, double*> faultParameters -> allocated matrices in Fortran with C++ pointer from Easi inread.
 * from interoperabilty functions that access Fortran values and copy them in C++ memory (copy by value) (e.g. e_interoperability.getDynRupParameters())
 */
class seissol::dr::initializers::BaseDRInitializer {
protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
  //YAML::Node m_InputParam;
  dr::DRParameters *m_Params;

public:
  virtual ~BaseDRInitializer() {}

  //set the parameters from .par file with yaml to this class attributes.
  void setInputParam(dr::DRParameters *DynRupParameter);

  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          seissol::initializers::LTSTree* dynRupTree,
                                          seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
                                          std::unordered_map<std::string, double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability);
};


#endif //SEISSOL_BASEDRINITIALIZER_H
