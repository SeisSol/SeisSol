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
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "generated_code/tensor.h"
#include <Eigen/Dense>
#include <complex>
#include <cstddef>
#include <vector>

namespace seissol {

namespace kernels {
constexpr std::size_t NumSpaceQuadraturePoints = (ConvergenceOrder + 1) * (ConvergenceOrder + 1);
} // namespace kernels

struct GlobalData {
  /**
   * Addresses of the global change of basis matrices (multiplied by the inverse diagonal mass
   *matrix):
   *
   *    0: \f$ M^{-1} R^1 \f$
   *    1: \f$ M^{-1} R^2 \f$
   *    2: \f$ M^{-1} R^3 \f$
   *    3: \f$ M^{-1} R^4 \f$
   **/
  seissol::tensor::rDivM::Container<const real*> changeOfBasisMatrices;

  /**
   * Addresses of the transposed global change of basis matrices left-multiplied with the local flux
   *matrix:
   *
   *    0: \f$ F^- ( R^1 )^T \f$
   *    1: \f$ F^- ( R^2 )^T \f$
   *    2: \f$ F^- ( R^3 )^T \f$
   *    3: \f$ F^- ( R^4 )^T \f$
   **/
  seissol::tensor::fMrT::Container<const real*> localChangeOfBasisMatricesTransposed;

  /**
   * Addresses of the transposed global change of basis matrices:
   *
   *    0: \f$ ( R^1 )^T \f$
   *    1: \f$ ( R^2 )^T \f$
   *    2: \f$ ( R^3 )^T \f$
   *    3: \f$ ( R^4 )^T \f$
   **/
  seissol::tensor::rT::Container<const real*> neighborChangeOfBasisMatricesTransposed;

  /**
   * Addresses of the global flux matrices:
   *
   *    0: \f$ F^{+,1} \f$
   *    1: \f$ F^{+,2} \f$
   *    2: \f$ F^{+,3} \f$
   **/
  seissol::tensor::fP::Container<const real*> neighborFluxMatrices;

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
  seissol::tensor::kDivM::Container<const real*> stiffnessMatrices;

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
  seissol::tensor::kDivMT::Container<const real*> stiffnessMatricesTransposed;

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
  seissol::tensor::V3mTo2nTWDivM::Container<const real*> nodalFluxMatrices;

  seissol::nodal::tensor::V3mTo2nFace::Container<const real*> v3mTo2nFace;
  seissol::tensor::project2nFaceTo3m::Container<const real*> project2nFaceTo3m;

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
  seissol::tensor::plusFluxMatrices::Container<real const*> plusFluxMatrices;
  seissol::tensor::minusFluxMatrices::Container<const real*> minusFluxMatrices;
#endif // ACL_DEVICE

  seissol::tensor::V3mTo2n::Container<real const*> faceToNodalMatrices;

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
};

struct CompoundGlobalData {
  GlobalData* onHost{nullptr};
  GlobalData* onDevice{nullptr};
};

// data for the cell local integration
struct LocalIntegrationData {
  // star matrices
  real starMatrices[3][seissol::tensor::star::size(0)];

  // flux solver for element local contribution
  real nApNm1[4][seissol::tensor::AplusT::size()];

  // equation-specific data
  seissol::model::MaterialT::LocalSpecificData specific;
};

// data for the neighboring boundary integration
struct NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  real nAmNm1[4][seissol::tensor::AminusT::size()];

  // equation-specific data
  seissol::model::MaterialT::NeighborSpecificData specific;
};

// material constants per cell
struct CellMaterialData {
  seissol::model::MaterialT local;
  seissol::model::MaterialT neighbor[4];
};

struct DRFaceInformation {
  std::size_t meshFace;
  std::uint8_t plusSide;
  std::uint8_t minusSide;
  std::uint8_t faceRelation;
  bool plusSideOnThisRank;
};

struct DRGodunovData {
  real dataTinvT[seissol::tensor::TinvT::size()];
  real tractionPlusMatrix[seissol::tensor::tractionPlusMatrix::size()];
  real tractionMinusMatrix[seissol::tensor::tractionMinusMatrix::size()];
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

struct DREnergyOutput {
  real slip[seissol::tensor::slipInterpolated::size()];
  real accumulatedSlip[seissol::dr::misc::NumPaddedPoints];
  real frictionalEnergy[seissol::dr::misc::NumPaddedPoints];
  real timeSinceSlipRateBelowThreshold[seissol::dr::misc::NumPaddedPoints];

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

struct CellDRMapping {
  unsigned side;
  unsigned faceRelation;
  real* godunov;
  real* fluxSolver;
};

struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  real nodes[seissol::nodal::tensor::nodes2D::Shape[0] * 3]{};
  real dataT[seissol::tensor::T::size()]{};
  real dataTinv[seissol::tensor::Tinv::size()]{};
  real easiBoundaryConstant[seissol::tensor::easiBoundaryConstant::size()]{};
  real easiBoundaryMap[seissol::tensor::easiBoundaryMap::size()]{};
};

struct CellBoundaryMapping {
  real* nodes{nullptr};
  real* dataT{nullptr};
  real* dataTinv{nullptr};
  real* easiBoundaryConstant{nullptr};
  real* easiBoundaryMap{nullptr};

  CellBoundaryMapping() = default;
  CellBoundaryMapping(BoundaryFaceInformation& faceInfo)
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
