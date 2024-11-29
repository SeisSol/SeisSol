#include <Kernels/Precision.h>
#include <init.h>
#include <omp.h>
#include <tensor.h>
#include <yateto.h>
#include <cstdio>
#include <cmath>

namespace seissol::kernels::local_flux::aux::details {

void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream) {

  const size_t workGroupSize = yateto::leadDim<seissol::nodal::init::nodes2D>();

  auto stream = reinterpret_cast<int*>(deviceStream);

  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t elementId = 0; elementId < numElements; ++elementId) {
    #pragma omp parallel for
    for (size_t tid = 0; tid < workGroupSize; ++tid) {
      const double rho = rhos[elementId];
      real* elementBoundaryDofs = dofsFaceBoundaryNodalPtrs[elementId];
      real* elementDisplacement = displacementDataPtrs[elementId];

      constexpr auto numNodes = seissol::nodal::tensor::nodes2D::Shape[0];
      if (tid < numNodes) {
        constexpr auto ldINodal = yateto::leadDim<seissol::init::INodal>();

        const auto pressureAtBnd = static_cast<real>(-1.0) * rho * g * elementDisplacement[tid];

        #pragma unroll
        for (int component{0}; component < 3; ++component) {
          elementBoundaryDofs[tid + component * ldINodal] =
            2.0 * pressureAtBnd - elementBoundaryDofs[tid + component * ldINodal];
        }
      }
    }
  }
}

void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  const size_t workGroupSize = yateto::leadDim<seissol::init::INodal>();

  constexpr auto ldINodalDim = yateto::leadDim<seissol::init::INodal>();
  constexpr auto INodalDim0 = seissol::tensor::INodal::Shape[0];
  constexpr auto INodalDim1 = seissol::tensor::INodal::Shape[1];

  constexpr auto ldConstantDim = yateto::leadDim<seissol::init::easiBoundaryConstant>();
  constexpr auto ConstantDim0 = seissol::tensor::easiBoundaryConstant::Shape[0];
  constexpr auto ConstantDim1 = seissol::tensor::easiBoundaryConstant::Shape[1];

  constexpr auto ldMapDim = yateto::leadDim<seissol::init::easiBoundaryMap>();
  constexpr auto MapDim0 = seissol::tensor::easiBoundaryMap::Shape[0];
  constexpr auto MapDim1 = seissol::tensor::easiBoundaryMap::Shape[1];
  constexpr auto MapDim2 = seissol::tensor::easiBoundaryMap::Shape[2];

  static_assert(INodalDim1 == ConstantDim0, "supposed to be equal");
  static_assert(INodalDim1 == MapDim0, "supposed to be equal");

  auto stream = reinterpret_cast<int*>(deviceStream);

  #pragma omp target teams nowait depend(inout: stream[0]) num_teams(numElements)
  {
    real resultTerm[INodalDim1][INodalDim0];
    real rightTerm[INodalDim1][ldConstantDim];
    real leftTerm[MapDim0][MapDim2];
    #pragma omp parallel num_threads(workGroupSize)
    {
      const int tid = omp_get_thread_num();
      const int elementId = omp_get_team_num();

      if (elementId < numElements) {
        real* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
        real* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
        auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

        for (int i = tid; i < (ldConstantDim * ConstantDim1); i += workGroupSize) {
          const auto b = i % ldConstantDim;
          const auto l = i / ldConstantDim;
          rightTerm[b][l] = easiBoundaryConstant[i];
        }
        #pragma omp barrier

        for (int i = 0; i < INodalDim1; ++i) {
          if (tid < INodalDim0) resultTerm[i][tid] = 0.0;
        }
        #pragma omp barrier

        for (int b = 0; b < MapDim1; ++b) {
          for (int l = 0; l < MapDim2; ++l) {
            if (tid < MapDim0) {
              leftTerm[tid][l] = easiBoundaryMap[tid + ldMapDim * (b + l * MapDim1)];
            }
          }
          #pragma omp barrier

          if (tid < MapDim2) {
            const real col = dofsFaceBoundaryNodal[tid + b * ldINodalDim];
            for (int a = 0; a < MapDim0; ++a) {
              resultTerm[a][tid] += leftTerm[a][tid] * col;
            }
          }
          #pragma omp barrier
        }

        if (tid < INodalDim0) {
          for (int a = 0; a < INodalDim1; ++a) {
            dofsFaceBoundaryNodal[tid + a * ldINodalDim] = resultTerm[a][tid] + rightTerm[a][tid];
          }
        }
      }
    }
  }
}
} // seissol::kernels::local_flux::aux::details


namespace seissol::kernels::time::aux {
void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream) {
  auto stream = reinterpret_cast<int*>(deviceStream);

  constexpr size_t workGroupSize = 9;
  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t elementId = 0; elementId < numElements; ++elementId) {
    #pragma omp parallel for
    for (int tid = 0; tid < workGroupSize; ++tid) {
      real* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
      real* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
      auto* T = TPtrs[elementId];
      auto* Tinv = TinvPtrs[elementId];

      constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv>();
      constexpr auto ldT = yateto::leadDim<seissol::init::T>();
      constexpr auto ldDisplacement = yateto::leadDim<seissol::init::displacementRotationMatrix>();

      const int i = tid % 3;
      const int j = tid / 3;

      displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
      displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
    }
  }
}

void initializeTaylorSeriesForGravitationalBoundary(
  real** prevCoefficientsPtrs,
  real** integratedDisplacementNodalPtrs,
  real** rotatedFaceDisplacementPtrs,
  double deltaTInt,
  size_t numElements,
  void* deviceStream) {

  const size_t workGroupSize = yateto::leadDim<seissol::nodal::init::nodes2D>();
  auto stream = reinterpret_cast<int*>(deviceStream);

  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t elementId = 0; elementId < numElements; ++elementId) {
    #pragma omp parallel for
    for (size_t tid = 0; tid < workGroupSize; ++tid) {
      auto* prevCoefficients = prevCoefficientsPtrs[elementId];
      auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
      const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

      assert(nodal::tensor::nodes2D::Shape[0] <= yateto::leadDim<seissol::init::rotatedFaceDisplacement>());

      constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];
      if (tid < num2dNodes) {
        prevCoefficients[tid] = rotatedFaceDisplacement[tid];
        integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
      }
    }
  }
}

void computeInvAcousticImpedance(double* invImpedances,
                                 double* rhos,
                                 double* lambdas,
                                 size_t numElements,
                                 void* deviceStream) {

  auto stream = reinterpret_cast<int*>(deviceStream);

  #pragma omp target teams distribute parallel for depend(inout: stream[0]) nowait
  for (size_t index = 0; index < numElements; ++index) {
      invImpedances[index] = 1.0 / std::sqrt(lambdas[index] * rhos[index]);
    }
}

void updateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                   real** prevCoefficientsPtrs,
                                   real** integratedDisplacementNodalPtrs,
                                   real** dofsFaceNodalPtrs,
                                   double* invImpedances,
                                   double* rhos,
                                   double g,
                                   double factorEvaluated,
                                   double factorInt,
                                   size_t numElements,
                                   void* deviceStream) {

  const size_t workGroupSize = yateto::leadDim<seissol::nodal::init::nodes2D>();

  auto stream = reinterpret_cast<int*>(deviceStream);

  #pragma omp target teams distribute depend(inout: stream[0]) nowait
  for (size_t elementId = 0; elementId < numElements; ++elementId) {
    #pragma omp parallel for
    for (size_t tid = 0; tid < workGroupSize; ++tid) {
      constexpr int pIdx = 0;
      constexpr int uIdx = 6;
      constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];

      if (tid < num2dNodes) {

        real* dofsFaceNodal = dofsFaceNodalPtrs[elementId];
        constexpr auto ldINodal = yateto::leadDim<seissol::init::INodal>();

        const auto uInside = dofsFaceNodal[tid + (uIdx + 0) * ldINodal];
        const auto vInside = dofsFaceNodal[tid + (uIdx + 1) * ldINodal];
        const auto wInside = dofsFaceNodal[tid + (uIdx + 2) * ldINodal];
        const auto pressureInside = dofsFaceNodal[tid + pIdx * ldINodal];

        real* prevCoefficients = prevCoefficientsPtrs[elementId];
#ifdef USE_ELASTIC
        const auto rho = rhos[elementId];
        const auto invImpedance = invImpedances[elementId];

        const double curCoeff = uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
#else
        const double curCoeff = uInside;
#endif
        prevCoefficients[tid] = curCoeff;

        constexpr auto ldFaceDisplacement = yateto::leadDim<seissol::init::faceDisplacement>();
        static_assert(num2dNodes <= ldFaceDisplacement, "");

        real* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];
        rotatedFaceDisplacement[tid + 0 * ldFaceDisplacement] += factorEvaluated * curCoeff;
        rotatedFaceDisplacement[tid + 1 * ldFaceDisplacement] += factorEvaluated * vInside;
        rotatedFaceDisplacement[tid + 2 * ldFaceDisplacement] += factorEvaluated * wInside;

        constexpr auto ldIntegratedFaceDisplacement = yateto::leadDim<seissol::init::averageNormalDisplacement>();
        static_assert(num2dNodes <= ldIntegratedFaceDisplacement, "");

        real* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
        integratedDisplacementNodal[tid] += factorInt * curCoeff;
      }
    }
  }
}
} // namespace seissol::kernels::time::aux
