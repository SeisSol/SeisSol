// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_TYPEDEFS_H_
#define SEISSOL_SRC_INITIALIZER_TYPEDEFS_H_

#include "BasicTypedefs.h"
#include "CellLocalInformation.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include <Eigen/Dense>
#include <Solver/MultipleSimulations.h>
#include <complex>
#include <cstddef>
#include <vector>

namespace seissol {

namespace kernels {
template <typename Config>
constexpr std::size_t NumSpaceQuadraturePoints =
    (Config::ConvergenceOrder + 1) * (Config::ConvergenceOrder + 1);
} // namespace kernels

template <typename Cfg>
struct GlobalDataCfg {
  using real = Real<Cfg>;
  /**
   * Addresses of the global change of basis matrices (multiplied by the inverse diagonal mass
   *matrix):
   *
   *    0: \f$ M^{-1} R^1 \f$
   *    1: \f$ M^{-1} R^2 \f$
   *    2: \f$ M^{-1} R^3 \f$
   *    3: \f$ M^{-1} R^4 \f$
   **/
  typename seissol::tensor::rDivM<Cfg>::template Container<const real*> changeOfBasisMatrices;

  /**
   * Addresses of the transposed global change of basis matrices left-multiplied with the local flux
   *matrix:
   *
   *    0: \f$ F^- ( R^1 )^T \f$
   *    1: \f$ F^- ( R^2 )^T \f$
   *    2: \f$ F^- ( R^3 )^T \f$
   *    3: \f$ F^- ( R^4 )^T \f$
   **/
  typename seissol::tensor::fMrT<Cfg>::template Container<const real*>
      localChangeOfBasisMatricesTransposed;

  /**
   * Addresses of the transposed global change of basis matrices:
   *
   *    0: \f$ ( R^1 )^T \f$
   *    1: \f$ ( R^2 )^T \f$
   *    2: \f$ ( R^3 )^T \f$
   *    3: \f$ ( R^4 )^T \f$
   **/
  typename seissol::tensor::rT<Cfg>::template Container<const real*>
      neighborChangeOfBasisMatricesTransposed;

  /**
   * Addresses of the global flux matrices:
   *
   *    0: \f$ F^{+,1} \f$
   *    1: \f$ F^{+,2} \f$
   *    2: \f$ F^{+,3} \f$
   **/
  typename seissol::tensor::fP<Cfg>::template Container<const real*> neighborFluxMatrices;

  /**
   * Addresses of the global stiffness matrices (multiplied by the inverse diagonal mass matrix):
   *
   *    0:  \f$ M^{-1} K^\xi \f$
   *    1:  \f$ M^{-1} K^\eta \f$
   *    2:  \f$ M^{-1} K^\zeta f$
   *
   *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks
   *(except for the additional flux matrix).
   **/
  typename seissol::tensor::kDivM<Cfg>::template Container<const real*> stiffnessMatrices;

  /**
   * Addresses of the transposed global stiffness matrices (multiplied by the inverse diagonal mass
   *matrix):
   *
   *    0:  \f$ M^{-1} ( K^\xi )^T \f$
   *    1:  \f$ M^{-1} ( K^\eta )^T \f$
   *    2:  \f$ M^{-1} ( K^\zeta )^T \f$
   *
   *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks
   *(except for the additional flux matrix).
   **/
  typename seissol::tensor::kDivMT<Cfg>::template Container<const real*>
      stiffnessMatricesTransposed;

  /**
   * Address of the (thread-local) local time stepping integration buffers used in the neighbor
   *integral computation
   **/
  real* integrationBufferLTS{nullptr};

  /**
   * Addresses of the global nodal flux matrices
   *
   *    0:  \f$ P^{+,1} \f$
   *    1:  \f$ P^{-,1,1} \f$
   *    2:  \f$ P^{-,1,2} \f$
   *    3 : \f$ P^{-,1,3} \f$
   *    4:  \f$ P^{+,2} \f$
   *    5:  \f$ P^{-,2,1} \f$
   *    6:  \f$ P^{-,2,2} \f$
   *    7 : \f$ P^{-,2,3} \f$
   *    [..]
   *    15: \f$ P^{-,4,3} \f$
   **/
  typename seissol::tensor::V3mTo2nTWDivM<Cfg>::template Container<const real*> nodalFluxMatrices;

  typename seissol::nodal::tensor::V3mTo2nFace<Cfg>::template Container<const real*> v3mTo2nFace;
  typename seissol::tensor::project2nFaceTo3m<Cfg>::template Container<const real*>
      project2nFaceTo3m;

  /**
   * Addresses of the global face to nodal matrices
   *
   *    0:  \f$ N^{+,1} \f$
   *    1:  \f$ N^{-,1,1} \f$
   *    2:  \f$ N^{-,1,2} \f$
   *    3 : \f$ N^{-,1,3} \f$
   *    4:  \f$ N^{+,2} \f$
   *    5:  \f$ N^{-,2,1} \f$
   *    6:  \f$ N^{-,2,2} \f$
   *    7 : \f$ N^{-,2,3} \f$
   *    [..]
   *    15: \f$ N^{-,4,3} \f$
   **/

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  typename seissol::tensor::plusFluxMatrices<Cfg>::template Container<real const*> plusFluxMatrices;
  typename seissol::tensor::minusFluxMatrices<Cfg>::template Container<const real*>
      minusFluxMatrices;
#endif // ACL_DEVICE

  typename seissol::tensor::V3mTo2n<Cfg>::template Container<real const*> faceToNodalMatrices;

  //! Modal basis to quadrature points
  real* evalAtQPMatrix{nullptr};

  //! Project function evaluated at quadrature points to modal basis
  real* projectQPMatrix{nullptr};

  //! Switch to nodal for plasticity
  real* vandermondeMatrix{nullptr};
  real* vandermondeMatrixInverse{nullptr};

  // A vector of ones. Note: It is only relevant for GPU computing.
  // It allows us to allocate this vector only once in the GPU memory
  real* replicateStresses{nullptr};

  real* selectAne{nullptr};
  real* selectEla{nullptr};

  real* resampleMatrix{nullptr};
  real* spaceWeights{nullptr};
  real* tpInverseFourierCoefficients{nullptr};
  real* tpGridPoints{nullptr};
  real* heatSource{nullptr};
};

// data for the cell local integration
template <typename Cfg>
struct LocalIntegrationData {
  // star matrices
  Real<Cfg> starMatrices[Cell::Dim][seissol::tensor::star<Cfg>::size(0)];

  // flux solver for element local contribution
  Real<Cfg> nApNm1[Cell::NumFaces][seissol::tensor::AplusT<Cfg>::size()];

  // equation-specific data
  typename seissol::model::MaterialTT<Cfg>::LocalSpecificData specific;
};

// data for the neighboring boundary integration
template <typename Cfg>
struct NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  Real<Cfg> nAmNm1[Cell::NumFaces][seissol::tensor::AminusT<Cfg>::size()];

  // equation-specific data
  typename seissol::model::MaterialTT<Cfg>::NeighborSpecificData specific;
};

// material constants per cell
struct CellMaterialData {
  seissol::model::Material* local;
  seissol::model::Material* neighbor[4];
};

struct DRFaceInformation {
  std::size_t meshFace;
  std::uint8_t plusSide;
  std::uint8_t minusSide;
  std::uint8_t faceRelation;
  bool plusSideOnThisRank;
};

template <typename Cfg>
struct DRGodunovData {
  Real<Cfg> dataTinvT[seissol::tensor::TinvT<Cfg>::size()];
  Real<Cfg> tractionPlusMatrix[seissol::tensor::tractionPlusMatrix<Cfg>::size()];
  Real<Cfg> tractionMinusMatrix[seissol::tensor::tractionMinusMatrix<Cfg>::size()];
  // When integrating quantities over the fault
  // we need to integrate over each physical element.
  // The integration is effectively done in the reference element, and the scaling factor of
  // the transformation, the surface Jacobian (e.g. |n^e(\chi)| in eq. (35) of Uphoff et al. (2023))
  // is incorporated. This explains the factor 2 (doubledSurfaceArea)
  //
  // Uphoff, C., May, D. A., & Gabriel, A. A. (2023). A discontinuous Galerkin method for
  // sequences of earthquakes and aseismic slip on multiple faults using unstructured curvilinear
  // grids. Geophysical Journal International, 233(1), 586-626.
  double doubledSurfaceArea;
};

template <typename Cfg>
struct DREnergyOutput {
  Real<Cfg> slip[seissol::tensor::slipInterpolated<Cfg>::size()];
  Real<Cfg> accumulatedSlip[seissol::dr::misc::NumPaddedPoints<Cfg>];
  Real<Cfg> frictionalEnergy[seissol::dr::misc::NumPaddedPoints<Cfg>];
  Real<Cfg> timeSinceSlipRateBelowThreshold[seissol::dr::misc::NumPaddedPoints<Cfg>];

  static std::vector<seissol::io::datatype::StructDatatype::MemberInfo> datatypeLayout() {
    return {
        seissol::io::datatype::StructDatatype::MemberInfo{
            "slip",
            offsetof(DREnergyOutput, slip),
            seissol::io::datatype::inferDatatype<decltype(slip)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "accumulatedSlip",
            offsetof(DREnergyOutput, accumulatedSlip),
            seissol::io::datatype::inferDatatype<decltype(accumulatedSlip)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "frictionalEnergy",
            offsetof(DREnergyOutput, frictionalEnergy),
            seissol::io::datatype::inferDatatype<decltype(frictionalEnergy)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "timeSinceSlipRateBelowThreshold",
            offsetof(DREnergyOutput, timeSinceSlipRateBelowThreshold),
            seissol::io::datatype::inferDatatype<decltype(timeSinceSlipRateBelowThreshold)>()},
    };
  }
};

template <typename Cfg>
struct CellDRMapping {
  unsigned side;
  unsigned faceRelation;
  Real<Cfg>* godunov;
  Real<Cfg>* fluxSolver;
};

template <typename Cfg>
struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  Real<Cfg>
      nodes[seissol::nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisFunctionDimension] * 3]{};
  Real<Cfg> dataT[seissol::tensor::T<Cfg>::size()]{};
  Real<Cfg> dataTinv[seissol::tensor::Tinv<Cfg>::size()]{};
  Real<Cfg> easiBoundaryConstant[seissol::tensor::easiBoundaryConstant<Cfg>::size()]{};
  Real<Cfg> easiBoundaryMap[seissol::tensor::easiBoundaryMap<Cfg>::size()]{};
};

template <typename Cfg>
struct CellBoundaryMapping {
  Real<Cfg>* nodes{nullptr};
  Real<Cfg>* dataT{nullptr};
  Real<Cfg>* dataTinv{nullptr};
  Real<Cfg>* easiBoundaryConstant{nullptr};
  Real<Cfg>* easiBoundaryMap{nullptr};

  CellBoundaryMapping() = default;
  CellBoundaryMapping(BoundaryFaceInformation<Cfg>& faceInfo)
      : nodes(faceInfo.nodes), dataT(faceInfo.dataT), dataTinv(faceInfo.dataTinv),
        easiBoundaryConstant(faceInfo.easiBoundaryConstant),
        easiBoundaryMap(faceInfo.easiBoundaryMap) {}
};

struct GravitationSetup {
  double acceleration = 9.81; // m/s
};

struct TravellingWaveParameters {
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  std::vector<int> varField;
  std::vector<std::complex<double>> ampField;
};

struct AcousticTravellingWaveParametersITM {
  double k;
  double itmStartingTime;
  double itmDuration;
  double itmVelocityScalingFactor;
};

struct PressureInjectionParameters {
  std::array<double, 3> origin;
  double magnitude;
  double width;
};

} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_TYPEDEFS_H_
