// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_GRAVITATIONALFREESURFACEBC_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_GRAVITATIONALFREESURFACEBC_H_

#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Numerical/ODEInt.h"
#include "Numerical/Quadrature.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#include <utility>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/LinearCK/DeviceAux/KernelsAux.h"

#include <Device/device.h>
#include <tuple>
#endif

namespace seissol {

class GravitationalFreeSurfaceBc {
  private:
  double gravitationalAcceleration;

  public:
  explicit GravitationalFreeSurfaceBc(double gravitationalAcceleration)
      : gravitationalAcceleration(gravitationalAcceleration) {};

  static std::pair<std::uint64_t, std::uint64_t>
      getFlopsDisplacementFace(unsigned face, [[maybe_unused]] FaceType faceType);

  template <typename MappingKrnl>
  void evaluate(unsigned faceIdx,
                MappingKrnl&& fsgKernelBase,
                const CellBoundaryMapping& boundaryMapping,
                real* displacementNodalData,
                real* integratedDisplacementNodalData,
                const real* derivatives,
                const real* power,
                double timeStepWidth,
                CellMaterialData& materialData) {
    // This function does two things:
    // 1: Compute eta (for all three dimensions) at the end of the timestep
    // 2: Compute the integral of eta in normal direction over the timestep
    // We do this by building up the Taylor series of eta.
    // Eta is defined by the ODE eta_t = u^R - 1/Z * (rho g eta - p^R)
    // Compute coefficients by differentiating ODE recursively, e.g.:
    // eta_tt = u^R_t - 1/Z * (rho eta_t g - p^R_t)
    // and substituting the previous coefficient eta_t
    // This implementation sums up the Taylor series directly without storing
    // all coefficients.
    if constexpr (multisim::MultisimEnabled) {
      // ... or maybe it even does work now?
      logError() << "The Free Surface Gravity BC kernel does not work with multiple simulations "
                    "yet. Or at least, nobody has yet tested the new implementation.";
    } else {
      kernel::fsgKernel kernel = std::forward<MappingKrnl>(fsgKernelBase);

      // Prepare kernel that projects volume data to face and rotates it to face-nodal basis.
      assert(boundaryMapping.nodes != nullptr);
      assert(boundaryMapping.dataTinv != nullptr);
      assert(boundaryMapping.dataT != nullptr);

      kernel.T = boundaryMapping.dataT;
      kernel.Tinv = boundaryMapping.dataTinv;
      kernel.faceDisplacement = displacementNodalData;
      kernel.Iint = integratedDisplacementNodalData;

      double coeffTmp = timeStepWidth;

      for (std::size_t i = 0; i < ConvergenceOrder; ++i) {
        kernel.dQ(i) = derivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
        if (i > 0) {
          kernel.coeff(i) = coeffTmp;
          coeffTmp *= timeStepWidth / static_cast<double>(i);
        }
        kernel.power(i) = power[i];
      }

      const double rho = materialData.local->getDensity();
      const double g = gravitationalAcceleration; // [m/s^2]
      const double z = std::sqrt(materialData.local->getLambdaBar() * rho);
      const real factor = (rho * g / z);

      kernel.invImpFactor = &factor;

      kernel.execute(faceIdx);
    }
  }

#ifdef ACL_DEVICE
  template <typename TimeKrnl, typename MappingKrnl>
  void evaluateOnDevice(unsigned faceIdx,
                        MappingKrnl&& fsgKernelBase,
                        recording::ConditionalPointersToRealsTable& dataTable,
                        recording::ConditionalMaterialTable& materialTable,
                        double timeStepWidth,
                        device::DeviceInstance& device,
                        seissol::parallel::runtime::StreamRuntime& runtime) {

    using namespace seissol::recording;
    auto* deviceStream = runtime.stream();
    ConditionalKey key(
        *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, faceIdx);
    if (dataTable.find(key) != dataTable.end()) {
      const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getSize()};
      // const double g = gravitationalAcceleration;

      auto* invImpedances =
          materialTable[key].get(inner_keys::Material::Id::InvImpedances)->getDeviceDataPtr();

      auto* rhos = materialTable[key].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
      auto* lambdas = materialTable[key].get(inner_keys::Material::Id::Lambda)->getDeviceDataPtr();
      kernels::time::aux::computeInvAcousticImpedance(
          invImpedances, rhos, lambdas, numElements, deviceStream);

      auto** TinvDataPtrs = dataTable[key].get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();
      auto** TDataPtrs = dataTable[key].get(inner_keys::Wp::Id::T)->getDeviceDataPtr();
      auto** derivativesPtrs =
          dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getDeviceDataPtr();

      auto** displacementsPtrs =
          dataTable[key].get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
      auto** integratedDisplacementNodalPtrs =
          dataTable[key].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();

      kernel::gpu_fsgKernel kernel = std::forward<MappingKrnl>(fsgKernelBase);
      for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        kernel.dQ(i) = const_cast<const real**>(derivativesPtrs);
        kernel.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
      }
      kernel.Tinv = const_cast<const real**>(TinvDataPtrs);
      kernel.T = const_cast<const real**>(TDataPtrs);
      kernel.faceDisplacement = displacementsPtrs;
      kernel.Iint = integratedDisplacementNodalPtrs;

      double coeffInt = timeStepWidth;
      double coeff = 1;

      kernel.power(0) = timeStepWidth;
      for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        coeffInt *= timeStepWidth / (i + 1);
        coeff *= timeStepWidth / i;

        kernel.coeff(i) = coeff;
        kernel.power(i) = coeffInt;
      }

      kernel.streamPtr = deviceStream;
      kernel.execute(faceIdx);
    }
  }
#endif // ACL_DEVICE
};

} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_LINEARCK_GRAVITATIONALFREESURFACEBC_H_
