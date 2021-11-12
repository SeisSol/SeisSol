#ifndef SEISSOL_BASEDRINITIALIZER_H
#define SEISSOL_BASEDRINITIALIZER_H

#include <Solver/Interoperability.h>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/FrictionLaws/FrictionLaws.h"
#include "Initializer/InputAux.hpp"
#include "Initializer/ParameterDB.h"
#include "DynamicRupture/Parameters.h"

namespace seissol::dr::initializers {
class BaseDRInitializer;
}

/**
 * Base class for dynamic rupture initializers
 * This class reads space dependent parameters from an easi file through a FaultParameterDB and
 * global parameters from the parameters.par file Furthermore derived quantities (such as e.g.
 * intial friction) are computed.
 */
class seissol::dr::initializers::BaseDRInitializer {
  protected:
  /**
   * Number of quadrature points on an surface element
   */
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  /**
   * Number of quadrature points on an surface element padded for vectorization
   */
  static constexpr int numPaddedPoints = init::QInterpolated::Stop[0];

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
                               seissol::Interoperability* e_interoperability);

  protected:
  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::map<std::string, double*>, which maps the
   * parameter name, to the address in memory, where the parameter shall be stored
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  virtual void addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                                       seissol::initializers::DynamicRupture* dynRup,
                                       seissol::initializers::LTSInternalNode::leaf_iterator& it);

  private:
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
   * Reads the parameters from the easi file
   * @param faultParameterDB reference to a FaultParameterDB, which manages easi
   * @param faceIDs faceIDs of the cells which are to be read
   */
  void queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                  std::vector<unsigned> faceIDs);

  /**
   * Rotates the stress tensor a fault aligned coordinate system
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @param stressInFaultCS pointer to array of size [numCells][numPaddedPoints][6], stores rotated
   * stress
   * @param stressXX pointer to array of size[numCells][numPaddedPoints], stores the XX component of
   * the stress in cartesian coordinates
   * @param stressYY pointer to array of size[numCells][numPaddedPoints], stores the YY component of
   * the stress in cartesian coordinates
   * @param stressZZ pointer to array of size[numCells][numPaddedPoints], stores the ZZ component of
   * the stress in cartesian coordinates
   * @param stressXX pointer to array of size[numCells][numPaddedPoints], stores the XY component of
   * the stress in cartesian coordinates
   * @param stressYZ pointer to array of size[numCells][numPaddedPoints], stores the YZ component of
   * the stress in cartesian coordinates
   * @param stressXZ pointer to array of size[numCells][numPaddedPoints], stores the XZ component of
   * the stress in cartesian coordinates
   */
  void rotateStressToFaultCS(seissol::initializers::DynamicRupture* dynRup,
                             seissol::initializers::LTSTree::leaf_iterator& it,
                             real (*stressInFaultCS)[numPaddedPoints][6],
                             real (*stressXX)[numPaddedPoints],
                             real (*stressYY)[numPaddedPoints],
                             real (*stressZZ)[numPaddedPoints],
                             real (*stressXY)[numPaddedPoints],
                             real (*stressYZ)[numPaddedPoints],
                             real (*stressXZ)[numPaddedPoints]);
};

#endif // SEISSOL_BASEDRINITIALIZER_H
