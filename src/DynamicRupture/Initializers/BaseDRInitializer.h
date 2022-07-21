#ifndef SEISSOL_BASEDRINITIALIZER_H
#define SEISSOL_BASEDRINITIALIZER_H

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
   * not need to store values in the Fortran parts
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree);

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
   * traction1
   * traction2
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  void initializeOtherVariables(seissol::initializers::DynamicRupture* dynRup,
                                seissol::initializers::LTSInternalNode::leaf_iterator& it);

  /**
   * Reads the parameters from the easi file
   * @param faultParameterDB reference to a FaultParameterDB, which manages easi
   * @param faceIDs faceIDs of the cells which are to be read
   */
  void queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                  std::vector<unsigned> faceIDs);

  /**
   * Evaluates, whether the FaultParameterDB provides a certain parameter.
   * @param parameter
   * @return true if the FaultParameterDB provides parameter.
   */
  bool faultProvides(std::string&& parameter);

  private:
  /**
   * Rotates the fault-aligned traction to cartesian stress coordinates
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @param stressXX reference to std::vector of std::array<real, numPaddedPoints>
   * IN: stores traction in normal direction, OUT: stores the XX component of the stress in
   * cartesian coordinates
   * @param stressYY reference to std::vector of std::array<real, numPaddedPoints>
   * IN: zero, OUT: stores the YY component of the stress in cartesian coordinates
   * @param stressZZ reference to std::vector of std::array<real, numPaddedPoints>
   * IN: zero, OUT: stores the ZZ component of the stress in cartesian coordinates
   * @param stressXX reference to std::vector of std::array<real, numPaddedPoints>
   * IN: stores traction in strike direction, OUT: stores the XY component of the stress in
   * cartesian coordinates
   * @param stressYZ reference to std::vector of std::array<real, numPaddedPoints>
   * IN: zero, OUT: stores the YZ component of the stress in cartesian coordinates
   * @param stressXZ reference to std::vector of std::array<real, numPaddedPoints>
   * IN: stroes traction in dip direction, OUT: stores the XZ component of the stress in cartesian
   * coordinates
   */
  void rotateTractionToCartesianStress(
      seissol::initializers::DynamicRupture* dynRup,
      seissol::initializers::LTSTree::leaf_iterator& it,
      std::vector<std::array<real, misc::numPaddedPoints>>& stressXX,
      std::vector<std::array<real, misc::numPaddedPoints>>& stressYY,
      std::vector<std::array<real, misc::numPaddedPoints>>& stressZZ,
      std::vector<std::array<real, misc::numPaddedPoints>>& stressXY,
      std::vector<std::array<real, misc::numPaddedPoints>>& stressYZ,
      std::vector<std::array<real, misc::numPaddedPoints>>& stressXZ);

  /**
   * Rotates the stress tensor to a fault aligned coordinate system and stores it in stressInFaultCS
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

  /**
   * Checks how the initial stress (nucleation stress) is characterized. The user can either provide
   * the traction "T_n", "T_s, "T_d" ("Tnuc_n", "Tnuc_s", "Tnuc_d") or the full stress tensor
   * "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz" ("nuc_xx", "nuc_yy", "nuc_zz", "nuc_xy",
   * "nuc_yz", "nuc_xz"). The user either has to provide all three traction components or all six
   * stress components, but no mixture.
   * @param readNucleation if set to true, check the identifiers for the nucleation stress. If set
   * to false, check identifiers for the initial stress
   * @return vector of strings, with the identifiers for the initial stress.
   */
  std::vector<std::string> stressIdentifiers(bool readNucleation);
};

} // namespace seissol::dr::initializers
#endif // SEISSOL_BASEDRINITIALIZER_H
