// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_BASEDRINITIALIZER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_BASEDRINITIALIZER_H_

#include <Model/CommonDatastructures.h>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "Initializer/InputAux.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/SeisSolParameters.h"

namespace seissol {
class SeisSol;
namespace dr::initializer {
/**
 * Base class for dynamic rupture initializers
 * This class reads space dependent parameters from an easi file through a FaultParameterDB and
 * global parameters from the parameters.par file Furthermore derived quantities (such as e.g.
 * initial friction) are computed.
 */
class BaseDRInitializer {
  protected:
  seissol::SeisSol& seissolInstance;
  /**
   * reference to the dynamic rupture parameters, which describe the global behaviour
   */
  std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters;

  /**
   * Set which contains all parameteres, which are provided by easi for the fault
   */
  std::set<std::string> faultParameterNames;

  /**
   * Characterizes how the initial stress on the fault is provided.
   */
  enum class Parametrization { Traction, Cartesian };

  public:
  /**
   * @param drParameters reference to the DRParameters, which contain all information from the
   * DynamicRupture namelist in the parameters.par file
   */
  BaseDRInitializer(
      const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
      seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance), drParameters(drParameters),
        faultParameterNames(seissol::initializer::FaultParameterDB<double>::faultProvides(
            drParameters->faultFileName)) {};
  // (NOTE: the FaultParameterDB type didn't matter in the last line)

  virtual ~BaseDRInitializer() = default;

  /**
   * Main function to initialize all fault dependent parameters.
   * In particular this initializes the initial and nucleation stress and rotates them to the fault
   * aligned coordinate system Furthermore data is copied to the fortran part
   * @param drStorage pointer to the dynamic rupture storage
   * not need to store values in the Fortran parts
   */
  virtual void initializeFault(DynamicRupture::Storage& drStorage);

  protected:
  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param layer reference to a Storage layer
   */
  virtual void
      addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                              DynamicRupture::Layer& layer);

  /**
   * Finds all faceIDs in one iterator. This is the mapping idInStorage -> idInMesh
   * @param layer reference to a Storage layer
   * @return vector containing all faceIDs which are stored in the leaf_iterator
   */
  static std::vector<std::size_t> getFaceIDsInIterator(DynamicRupture::Layer& layer);

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
   * @param layer reference to a Storage layer
   */
  static void initializeOtherVariables(DynamicRupture::Layer& layer);

  /**
   * Reads the parameters from the easi file
   * @param faultParameterDB reference to a FaultParameterDB, which manages easi
   * @param faceIDs faceIDs of the cells which are to be read
   */
  template <typename T>
  void queryModel(seissol::initializer::FaultParameterDB<T>& faultParameterDB,
                  const std::vector<std::size_t>& faceIDs,
                  std::size_t simid,
                  std::size_t configId);

  /**
   * Evaluates, whether the FaultParameterDB provides a certain parameter.
   * @param parameter
   * @return true if the FaultParameterDB provides parameter.
   */
  bool faultProvides(const std::string& parameter);

  private:
  /**
   * Checks how the initial stress (nucleation stress) is characterized. The user can either provide
   * the traction "T_n", "T_s, "T_d" ("Tnuc_n", "Tnuc_s", "Tnuc_d") or the full stress tensor
   * "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz" ("nuc_xx", "nuc_yy", "nuc_zz", "nuc_xy",
   * "nuc_yz", "nuc_xz"). The user either has to provide all three traction components or all six
   * stress components, but no mixture.
   * @param readNucleation if set to > 0, check the identifiers for the nucleation stress. If set
   * to 0, check identifiers for the initial stress
   * @return vector of strings, with the identifiers for the initial stress.
   */
  std::pair<std::vector<std::string>, Parametrization>
      stressIdentifiers(int readNucleation, model::MaterialType materialType);
};

} // namespace dr::initializer
} // namespace seissol

#endif // SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_BASEDRINITIALIZER_H_
