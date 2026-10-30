// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_CURVEDSPLIT_H_
#define SEISSOL_SRC_KERNELS_CURVEDSPLIT_H_

#include "Common/Constants.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"

#include <array>
#include <cstddef>

namespace seissol::kernel {
/// The kernels of the corrector of a curved cell, which the generator writes where a cell may be
/// curved (Curvilinear); named here so that code for every build can name them.
struct splitVolume;
struct ownTraction;
} // namespace seissol::kernel

namespace seissol::kernels {

/**
 * The corrector of a curved cell (CellLocalInformation::curved).
 *
 * A curved cell is weight-adjusted with the weight J W, J its Jacobian determinant: it is tested
 * with the basis, its faces lift their fluxes without material and with the surface Jacobian alone
 * into its flux buffer, and the inverse weight it applies there carries 1 / J at every point (see
 * initializer::CurvedCell). Its volume term goes through the weight as well, in the split form:
 * the rows of the momentum in the weak form, the others in the strong form. In the energy, the two
 * volume terms then cancel at every point the rules have, and the faces leave the dissipation of
 * the flux; the strong form alone leaves a remainder there, the error of the rules on the
 * divergence theorem of a curved cell, with which it grows.
 *
 * So a curved cell leaves out the volume term of the predictor and forms its own from its time
 * integral into its flux buffer (kernel::splitVolume). The corrector of every face subtracts the
 * normal flux of the cell's own trace, as for every cell; the weak form of the momentum rows wants
 * the flux alone, so a curved cell adds that normal flux back in those rows, from the trace of
 * each face at the points of the face, with the scale the corrector has there
 * (kernel::ownTraction). A fault face leaves it out: the fault subtracts the own trace of a curved
 * side in the velocity alone (see FaultSubtractsOwnTrace). A constant state gives nothing, and on
 * a straight-sided cell the split form is the strong form.
 *
 * A class of its own, so that a build without the kernels never instantiates what names their
 * operands.
 */
template <bool Enabled, typename VolumeKernelT, typename OwnTractionKernelT>
class CurvedSplit {
  public:
  void setGlobalData(const CompoundGlobalData& global) { bindGlobals(*global.onHost); }

  void bindGlobals(const GlobalData& global) {
    volumePrototype_.bindGlobals(global);
    ownTractionPrototype_.bindGlobals(global);
  }

  /// Adds the volume term of the cell, from its time integral `timeIntegrated`, and the normal
  /// flux of its own trace on each face that is not a fault, to its flux buffer. The traces of the
  /// faces are the ones of this step (writeFaceTraces).
  void apply(const real* timeIntegrated, LTS::Ref& data) const {
    apply(timeIntegrated,
          data.get<LTS::LocalIntegration>(),
          data.get<LTS::CellInformation>(),
          data.get<LTS::FaceTraces>(),
          data.get<LTS::FluxBuffer>());
  }

  /// The same, with what the cell carries handed over one by one: `faceTraces` the traces of its
  /// faces, and `fluxBuffer` its flux buffer.
  void apply(const real* timeIntegrated,
             const LocalIntegrationData& localIntegration,
             const CellLocalInformation& cellInformation,
             const std::array<real*, Cell::NumFaces>& faceTraces,
             real* fluxBuffer) const {
    VolumeKernelT volume = volumePrototype_;
    volume.I = timeIntegrated;
    volume.fluxBuffer = fluxBuffer;
    for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
      volume.referenceGradients(dim) = localIntegration.referenceGradients[dim];
    }
    for (std::size_t a = 0; a < StarCoefficientCount; ++a) {
      volume.splitCoefficients(a) = localIntegration.splitCoefficients[a];
    }
    volume.execute();

    OwnTractionKernelT own = ownTractionPrototype_;
    own.fluxBuffer = fluxBuffer;
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (cellInformation.faceTypes[face] == FaceType::DynamicRupture) {
        continue;
      }
      own.faceTrace = faceTraces[face];
      bindFaceRotation(own, localIntegration.faceRotation[face]);
      for (std::size_t c = 0; c < OwnTractionCoefficientCount; ++c) {
        own.ownTractionCoefficients(OwnTractionCoefficientIndices[c]) =
            localIntegration.ownTractionCoefficients[face][c];
      }
      own.execute(face, cellInformation.traceSides[face]);
    }
  }

  /// What a curved cell runs on top of a straight-sided one in its local integral.
  [[nodiscard]] static PerformanceEstimate
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) {
    auto estimate = PerformanceEstimate::fromKernel<VolumeKernelT>();
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceTypes[face] != FaceType::DynamicRupture) {
        // both sides of a face cost the same
        estimate += PerformanceEstimate::fromKernel<OwnTractionKernelT>(face, 0);
      }
    }
    return estimate;
  }

  private:
  VolumeKernelT volumePrototype_;
  OwnTractionKernelT ownTractionPrototype_;
};

template <typename VolumeKernelT, typename OwnTractionKernelT>
class CurvedSplit<false, VolumeKernelT, OwnTractionKernelT> {
  public:
  void setGlobalData(const CompoundGlobalData& /*global*/) {}
  void bindGlobals(const GlobalData& /*global*/) {}
  void apply(const real* /*timeIntegrated*/, LTS::Ref& /*data*/) const {}
  void apply(const real* /*timeIntegrated*/,
             const LocalIntegrationData& /*localIntegration*/,
             const CellLocalInformation& /*cellInformation*/,
             const std::array<real*, Cell::NumFaces>& /*faceTraces*/,
             real* /*fluxBuffer*/) const {}
  [[nodiscard]] static PerformanceEstimate
      metrics(const std::array<FaceType, Cell::NumFaces>& /*faceTypes*/) {
    return {};
  }
};

/// The one of this build.
using CurvedSplitT = CurvedSplit<Curvilinear, kernel::splitVolume, kernel::ownTraction>;

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_CURVEDSPLIT_H_
