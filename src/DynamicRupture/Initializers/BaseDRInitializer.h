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
 * initial friction) are computed.
 */
class BaseDRInitializer {
  protected:
  /**
   * reference to the dynamic rupture parameters, which describe the global behaviour
   */
  std::shared_ptr<DRParameters> drParameters;

  /**
   * Set which contains all parameteres, which are provided by easi for the fault
   */
  std::set<std::string> faultParameterNames;

  /**
   * Characterizes how the initial stress on the fault is provided.
   */
  enum class Parametrization { Traction, Cartesian };

  /**
   * Stores the initialStresses.
   */
  struct StressTensor {
    StressTensor(size_t size) {
      xx.reserve(size);
      yy.reserve(size);
      zz.reserve(size);
      xy.reserve(size);
      yz.reserve(size);
      xz.reserve(size);
    }
    using VectorOfArrays_t = std::vector<std::array<real, misc::numPaddedPoints>>;
    VectorOfArrays_t xx;
    VectorOfArrays_t yy;
    VectorOfArrays_t zz;
    VectorOfArrays_t xy;
    VectorOfArrays_t yz;
    VectorOfArrays_t xz;
  };

  public:
  /**
   * @param drParameters reference to the DRParameters, which contain all information from the
   * DynamicRupture namelist in the parameters.par file
   */
  BaseDRInitializer(std::shared_ptr<DRParameters> drParameters)
      : drParameters(drParameters),
        faultParameterNames(
            seissol::initializers::FaultParameterDB::faultProvides(drParameters->faultFileName)){};

  virtual ~BaseDRInitializer() = default;

  /**
   * Main function to initialize all fault dependent parameters.
   * In particular this initializes the initial and nucleation stress and rotates them to the fault
   * aligned coordinate system Furthermore data is copied to the fortran part
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param dynRupTree pointer to the dynamic rupture lts tree
   * not need to store values in the Fortran parts
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture const* const dynRup,
                               seissol::initializers::LTSTree* const dynRupTree);

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
                              seissol::initializers::DynamicRupture const* const dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it);

  /**
   * Finds all faceIDs in one iterator. This is the mapping idInLTSTree -> idInMesh
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @return vector containing all faceIDs which are stored in the leaf_iterator
   */
  std::vector<unsigned>
      getFaceIDsInIterator(seissol::initializers::DynamicRupture const* const dynRup,
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
  void initializeOtherVariables(seissol::initializers::DynamicRupture const* const dynRup,
                                seissol::initializers::LTSInternalNode::leaf_iterator& it);

  /**
   * Reads the parameters from the easi file
   * @param faultParameterDB reference to a FaultParameterDB, which manages easi
   * @param faceIDs faceIDs of the cells which are to be read
   */
  void queryModel(seissol::initializers::FaultParameterDB& faultParameterDB,
                  std::vector<unsigned> const& faceIDs);

  /**
   * Evaluates, whether the FaultParameterDB provides a certain parameter.
   * @param parameter
   * @return true if the FaultParameterDB provides parameter.
   */
  bool faultProvides(const std::string& parameter);

  private:
  /**
   * Rotates the fault-aligned traction to cartesian stress coordinates
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @param stress reference to a StressTensor
   * IN: stores traction in fault strike/dip coordinate system OUT: stores the the stress in
   * cartesian coordinates
   */
  void rotateTractionToCartesianStress(seissol::initializers::DynamicRupture const* const dynRup,
                                       seissol::initializers::LTSTree::leaf_iterator& it,
                                       StressTensor& stress);

  /**
   * Rotates the stress tensor to a fault aligned coordinate system and stores it in stressInFaultCS
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   * @param stressInFaultCS pointer to array of size [numCells][numPaddedPoints][6], stores rotated
   * stress
   * @param stress reference to a StressTensor, stores the stress in cartesian coordinates
   */
  void rotateStressToFaultCS(seissol::initializers::DynamicRupture const* const dynRup,
                             seissol::initializers::LTSTree::leaf_iterator& it,
                             real (*stressInFaultCS)[misc::numPaddedPoints][6],
                             StressTensor const& stress);

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
  std::pair<std::vector<std::string>, Parametrization> stressIdentifiers(bool readNucleation);
};

} // namespace seissol::dr::initializers
#endif // SEISSOL_BASEDRINITIALIZER_H
