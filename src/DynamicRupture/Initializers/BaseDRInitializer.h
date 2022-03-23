#ifndef SEISSOL_BASEDRINITIALIZER_H
#define SEISSOL_BASEDRINITIALIZER_H

#include <Solver/Interoperability.h>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/InputAux.hpp"
#include "Initializer/ParameterDB.h"

namespace seissol::dr::initializers {
/**
 * Base class for dynamic rupture initializers
 * This class reads space dependent parameters from an easi file through a FaultParameterDB and
 * global parameters from the parameters.par file Furthermore derived quantities (such as e.g.
 * intial friction) are computed.
 */
class BaseDRInitializer {
  protected:
  /**
   * reference to the dynamic rupture parameters, which describe the global behaviour
   */
  DRParameters& drParameters;

  public:
  /**
   * @param drParameters reference to the DRParameters, which contain all information from the
   * DynamicRupture namelist in the parameters.par file
   */
  BaseDRInitializer(DRParameters& drParameters) : drParameters(drParameters){};

  virtual ~BaseDRInitializer() {}

  /**
   * Main function to initialize all fault dependent parameters.
   * In particular this initializes the initial and nucleation stress and rotates them to the fault
   * aligned coordinate system Furthermore data is copied to the fortran part
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param dynRupTree pointer to the dynamic rupture lts tree
   * @param e_interoperability pointer to the interoperability instance, can be removed once we do
   * not need to store values in the Fortran parts
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* eInteroperability);

  protected:
  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  virtual void
      addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it);

  /**
   * Finds all faceIDs in one iterator. This is the mapping idInLTSTree -> idInMesh
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @return vector containing all faceIDs which are stored in the leaf_iterator
   */
  std::vector<unsigned>
      getFaceIDsInIterator(seissol::initializers::DynamicRupture* dynRup,
                           seissol::initializers::LTSInternalNode::leaf_iterator& it);

  /**
   * Initialize all other variables:
   * ruptureTimePending
   * peakSlipRate
   * ruptureTime
   * dynStressTime
   * accumulatedSlipMagnitude
   * slip1
   * slip2
   * slipRateMagnitude
   * tractionXY
   * tractionXZ
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @param e_interoperability pointer to the interoperability instance, can be removed once we do
   * not need to store values in the Fortran parts
   */
  void initializeOtherVariables(seissol::initializers::DynamicRupture* dynRup,
                                seissol::initializers::LTSInternalNode::leaf_iterator& it,
                                Interoperability* eInteroperability);

  /**
   * Reads the parameters from the easi file
   * @param faultParameterDB reference to a FaultParameterDB, which manages easi
   * @param faceIDs faceIDs of the cells which are to be read
   */
  void queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                  std::vector<unsigned> faceIDs);

  private:
  /**
   * Rotates the stress tensor a fault aligned coordinate system
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @param stressInFaultCS pointer to array of size [numCells][numPaddedPoints][6], stores rotated
   * stress
   * @param stressXX reference to std::vector of std::array<real, numPaddedPoints>, stores the XX
   * component of the stress in cartesian coordinates
   * @param stressYY reference to std::vector of std::array<real, numPaddedPoints>, stores the YY
   * component of the stress in cartesian coordinates
   * @param stressZZ reference to std::vector of std::array<real, numPaddedPoints>, stores the ZZ
   * component of the stress in cartesian coordinates
   * @param stressXX reference to std::vector of std::array<real, numPaddedPoints>, stores the XY
   * component of the stress in cartesian coordinates
   * @param stressYZ reference to std::vector of std::array<real, numPaddedPoints>, stores the YZ
   * component of the stress in cartesian coordinates
   * @param stressXZ reference to std::vector of std::array<real, numPaddedPoints>, stores the XZ
   * component of the stress in cartesian coordinates
   */
  void rotateStressToFaultCS(seissol::initializers::DynamicRupture* dynRup,
                             seissol::initializers::LTSTree::leaf_iterator& it,
                             real (*stressInFaultCS)[misc::numPaddedPoints][6],
                             std::vector<std::array<real, misc::numPaddedPoints>>& stressXX,
                             std::vector<std::array<real, misc::numPaddedPoints>>& stressYY,
                             std::vector<std::array<real, misc::numPaddedPoints>>& stressZZ,
                             std::vector<std::array<real, misc::numPaddedPoints>>& stressXY,
                             std::vector<std::array<real, misc::numPaddedPoints>>& stressYZ,
                             std::vector<std::array<real, misc::numPaddedPoints>>& stressXZ);
};

} // namespace seissol::dr::initializers
#endif // SEISSOL_BASEDRINITIALIZER_H
