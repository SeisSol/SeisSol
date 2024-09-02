#include "Kernels/precision.hpp"
#include <init.h>
#include <tensor.h>
#include <yateto.h>
#include <cuda.h>
#include <cstdio>


namespace seissol::kernels::local_flux::aux::details {

__global__ void kernelFreeSurfaceGravity(
  real** dofsFaceBoundaryNodalPtrs,
  real** displacementDataPtrs,
  double* rhos,
  double g,
  size_t numElements) {

  const int tid = threadIdx.x;
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
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

void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream) {
  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelFreeSurfaceGravity<<<grid, block, 0, stream>>>(dofsFaceBoundaryNodalPtrs,
                                                       displacementDataPtrs,
                                                       rhos,
                                                       g,
                                                       numElements);
}


__global__ void  kernelEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                                    real** easiBoundaryMapPtrs,
                                    real** easiBoundaryConstantPtrs,
                                    size_t numElements) {

  const int tid = threadIdx.x;
  const int elementId = blockIdx.x;

  constexpr auto ldINodalDim = yateto::leadDim<seissol::init::INodal>();
  constexpr auto iNodalDim0 = seissol::tensor::INodal::Shape[0];
  constexpr auto iNodalDim1 = seissol::tensor::INodal::Shape[1];
  __shared__ __align__(8) real resultTerm[iNodalDim1][iNodalDim0];

  constexpr auto ldConstantDim = yateto::leadDim<seissol::init::easiBoundaryConstant>();
  constexpr auto constantDim0 = seissol::tensor::easiBoundaryConstant::Shape[0];
  constexpr auto constantDim1 = seissol::tensor::easiBoundaryConstant::Shape[1];
  __shared__ __align__(8) real rightTerm[iNodalDim1][ldConstantDim];

  constexpr auto ldMapDim = yateto::leadDim<seissol::init::easiBoundaryMap>();
  constexpr auto mapDim0 = seissol::tensor::easiBoundaryMap::Shape[0];
  constexpr auto mapDim1 = seissol::tensor::easiBoundaryMap::Shape[1];
  constexpr auto mapDim2 = seissol::tensor::easiBoundaryMap::Shape[2];
  __shared__ __align__(8) real leftTerm[mapDim0][mapDim2];

  static_assert(iNodalDim1 == constantDim0, "supposed to be equal");
  static_assert(iNodalDim1 == mapDim0, "supposed to be equal");


  if (elementId < numElements) {
    real* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
    real* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
    auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

    for (int i = tid; i < (ldConstantDim * constantDim1); i += blockDim.x) {
      const auto b = i % ldConstantDim;
      const auto l = i / ldConstantDim;
      rightTerm[b][l] = easiBoundaryConstant[i];
    }
    __syncthreads();

    for (int i = 0; i < iNodalDim1; ++i) {
      if (tid < iNodalDim0) resultTerm[i][tid] = 0.0;
    }
    __syncthreads();

    for (int b = 0; b < mapDim1; ++b) {
      for (int l = 0; l < mapDim2; ++l) {
        if (tid < mapDim0) {
          leftTerm[tid][l] = easiBoundaryMap[tid + ldMapDim * (b + l * mapDim1)];
        }
      }
      __syncthreads();

      if (tid < mapDim2) {
        const real col = dofsFaceBoundaryNodal[tid + b * ldINodalDim];
        for (int a = 0; a < mapDim0; ++a) {
          resultTerm[a][tid] += leftTerm[a][tid] * col;
        }
      }
      __syncthreads();
    }

    if (tid < iNodalDim0) {
      for (int a = 0; a < iNodalDim1; ++a) {
        dofsFaceBoundaryNodal[tid + a * ldINodalDim] = resultTerm[a][tid] + rightTerm[a][tid];
      }
    }
  }
}


void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  dim3 block(yateto::leadDim<seissol::init::INodal>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelEasiBoundary<<<grid, block, 0, stream>>>(dofsFaceBoundaryNodalPtrs,
                                                 easiBoundaryMapPtrs,
                                                 easiBoundaryConstantPtrs,
                                                 numElements);
}
} // namespace seissol::kernels::local_flux::aux::details



namespace seissol::kernels::time::aux {
__global__ void kernelextractRotationMatrices(real** displacementToFaceNormalPtrs,
                                              real** displacementToGlobalDataPtrs,
                                              real** TPtrs,
                                              real** TinvPtrs,
                                              size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    real* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
    real* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
    auto* T = TPtrs[elementId];
    auto* Tinv = TinvPtrs[elementId];

    constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv>();
    constexpr auto ldT = yateto::leadDim<seissol::init::T>();
    constexpr auto ldDisplacement = yateto::leadDim<seissol::init::displacementRotationMatrix>();

    const int i = threadIdx.x;
    const int j = threadIdx.y;

    displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
    displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
  }
}

void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream) {
  dim3 block(3, 3, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelextractRotationMatrices<<<grid, block, 0, stream>>>(
    displacementToFaceNormalPtrs,
    displacementToGlobalDataPtrs,
    TPtrs,
    TinvPtrs,
    numElements);
}

__global__  void kernelInitializeTaylorSeriesForGravitationalBoundary(
  real** prevCoefficientsPtrs,
  real** integratedDisplacementNodalPtrs,
  real** rotatedFaceDisplacementPtrs,
  double deltaTInt,
  size_t numElements) {

  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    auto* prevCoefficients = prevCoefficientsPtrs[elementId];
    auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
    const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

    assert(nodal::tensor::nodes2D::Shape[0] <= yateto::leadDim<seissol::init::rotatedFaceDisplacement>());

    const int tid = threadIdx.x;
    constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];
    if (tid < num2dNodes) {
      prevCoefficients[tid] = rotatedFaceDisplacement[tid];
      integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
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

  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelInitializeTaylorSeriesForGravitationalBoundary<<<grid, block, 0, stream>>>(
    prevCoefficientsPtrs,
    integratedDisplacementNodalPtrs,
    rotatedFaceDisplacementPtrs,
    deltaTInt,
    numElements);
}

__global__ void kernelComputeInvAcousticImpedance(double* invImpedances,
                                                  double* rhos,
                                                  double* lambdas,
                                                  size_t numElements) {

  size_t index = threadIdx.x + blockIdx.x * blockDim.x;
  if (index < numElements) {
    invImpedances[index] = 1.0 / std::sqrt(lambdas[index] * rhos[index]);
  }
}

void computeInvAcousticImpedance(double* invImpedances,
                                 double* rhos,
                                 double* lambdas,
                                 size_t numElements,
                                 void* deviceStream) {
  constexpr size_t blockSize{256};
  dim3 block(blockSize, 1, 1);
  dim3 grid((numElements + blockSize - 1) / blockSize, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelComputeInvAcousticImpedance<<<grid, block, 0, stream>>>(invImpedances,
                                                                rhos,
                                                                lambdas,
                                                                numElements);
}

__global__ void kernelUpdateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                                    real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** dofsFaceNodalPtrs,
                                                    double* invImpedances,
                                                    double* rhos,
                                                    double g,
                                                    double factorEvaluated,
                                                    double factorInt,
                                                    size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    constexpr int pIdx = 0;
    constexpr int uIdx = 6;
    constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];

    const int tid = threadIdx.x;
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
  dim3 block(yateto::leadDim<seissol::nodal::init::nodes2D>(), 1, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<cudaStream_t>(deviceStream);
  kernelUpdateRotatedFaceDisplacement<<<grid, block, 0, stream>>>(
    rotatedFaceDisplacementPtrs,
    prevCoefficientsPtrs,
    integratedDisplacementNodalPtrs,
    dofsFaceNodalPtrs,
    invImpedances,
    rhos,
    g,
    factorEvaluated,
    factorInt,
    numElements);
}
} // namespace seissol::kernels::time::aux
