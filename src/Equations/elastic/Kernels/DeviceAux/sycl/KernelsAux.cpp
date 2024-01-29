#include <CL/sycl.hpp>
#include <Kernels/precision.hpp>
#include <generated_code/init.h>
#include <generated_code/tensor.h>
#include <yateto.h>
#include <cstdio>

namespace seissol::kernels::local_flux::aux::details {

void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream) {

  auto queue = reinterpret_cast<cl::sycl::queue*>(deviceStream);
  const size_t workGroupSize = yateto::leadDim<seissol::nodal::init::nodes2D>();
  cl::sycl::nd_range rng{{numElements * workGroupSize}, {workGroupSize}};

  queue->parallel_for(rng, [=](cl::sycl::nd_item<1> item) {
    const int tid = item.get_local_id(0);
    const int elementId = item.get_group().get_group_id(0);
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
  });
}

void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  auto queue = reinterpret_cast<cl::sycl::queue*>(deviceStream);
  const size_t workGroupSize = yateto::leadDim<seissol::init::INodal>();
  cl::sycl::nd_range rng{{numElements * workGroupSize}, {workGroupSize}};

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

  using LocalMemoryType = cl::sycl::accessor<real, 2, cl::sycl::access_mode::read_write, cl::sycl::access::target::local>;

  queue->submit([&](cl::sycl::handler& cgh) {
    LocalMemoryType resultTerm(cl::sycl::range<2>(INodalDim1, INodalDim0), cgh);
    LocalMemoryType rightTerm(cl::sycl::range<2>(INodalDim1, ldConstantDim), cgh);
    LocalMemoryType leftTerm(cl::sycl::range<2>(MapDim0, MapDim2), cgh);

    cgh.parallel_for(rng, [=](cl::sycl::nd_item<1> item) {

      const int tid = item.get_local_id(0);
      const int elementId = item.get_group().get_group_id(0);

      if (elementId < numElements) {
        real* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
        real* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
        auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

        for (int i = tid; i < (ldConstantDim * ConstantDim1); i += item.get_local_range(0)) {
          const auto b = i % ldConstantDim;
          const auto l = i / ldConstantDim;
          rightTerm[b][l] = easiBoundaryConstant[i];
        }
        item.barrier();

        for (int i = 0; i < INodalDim1; ++i) {
          if (tid < INodalDim0) resultTerm[i][tid] = 0.0;
        }
        item.barrier();

        for (int b = 0; b < MapDim1; ++b) {
          for (int l = 0; l < MapDim2; ++l) {
            if (tid < MapDim0) {
              leftTerm[tid][l] = easiBoundaryMap[tid + ldMapDim * (b + l * MapDim1)];
            }
          }
          item.barrier();

          if (tid < MapDim2) {
            const real col = dofsFaceBoundaryNodal[tid + b * ldINodalDim];
            for (int a = 0; a < MapDim0; ++a) {
              resultTerm[a][tid] += leftTerm[a][tid] * col;
            }
          }
          item.barrier();
        }

        if (tid < INodalDim0) {
          for (int a = 0; a < INodalDim1; ++a) {
            dofsFaceBoundaryNodal[tid + a * ldINodalDim] = resultTerm[a][tid] + rightTerm[a][tid];
          }
        }
      }
    });
  });
}
} // seissol::kernels::local_flux::aux::details


namespace seissol::kernels::time::aux {
void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream) {
  auto queue = reinterpret_cast<cl::sycl::queue*>(deviceStream);

  constexpr size_t workGroupSize = 9;
  cl::sycl::nd_range<1> rng{{numElements * workGroupSize}, {workGroupSize}};

  queue->parallel_for(rng, [=](cl::sycl::nd_item<1> item) {
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      real* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
      real* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
      auto* T = TPtrs[elementId];
      auto* Tinv = TinvPtrs[elementId];

      constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv>();
      constexpr auto ldT = yateto::leadDim<seissol::init::T>();
      constexpr auto ldDisplacement = yateto::leadDim<seissol::init::displacementRotationMatrix>();

      const int i = item.get_local_id(0) % 3;
      const int j = item.get_local_id(0) / 3;

      displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
      displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
    }
  });
}

void initializeTaylorSeriesForGravitationalBoundary(
  real** prevCoefficientsPtrs,
  real** integratedDisplacementNodalPtrs,
  real** rotatedFaceDisplacementPtrs,
  double deltaTInt,
  size_t numElements,
  void* deviceStream) {

  auto queue = reinterpret_cast<cl::sycl::queue*>(deviceStream);
  const size_t workGroupSize = yateto::leadDim<seissol::nodal::init::nodes2D>();
  cl::sycl::nd_range rng{{numElements * workGroupSize}, {workGroupSize}};

  queue->parallel_for(rng, [=](cl::sycl::nd_item<1> item) {
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      auto* prevCoefficients = prevCoefficientsPtrs[elementId];
      auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
      const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

      assert(nodal::tensor::nodes2D::Shape[0] <= yateto::leadDim<seissol::init::rotatedFaceDisplacement>());

      const int tid = item.get_local_id(0);
      constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];
      if (tid < num2dNodes) {
        prevCoefficients[tid] = rotatedFaceDisplacement[tid];
        integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
      }
    }
  });
}

void computeInvAcousticImpedance(double* invImpedances,
                                 double* rhos,
                                 double* lambdas,
                                 size_t numElements,
                                 void* deviceStream) {
  constexpr size_t blockSize{256};
  auto queue = reinterpret_cast<cl::sycl::queue*>(deviceStream);
  cl::sycl::nd_range rng{{numElements * blockSize}, {blockSize}};

  queue->parallel_for(rng, [=](cl::sycl::nd_item<1> item) {
    size_t index = item.get_global_id(0);
    if (index < numElements) {
      invImpedances[index] = 1.0 / cl::sycl::sqrt(lambdas[index] * rhos[index]);
    }
  });
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

  auto queue = reinterpret_cast<cl::sycl::queue*>(deviceStream);
  const size_t workGroupSize = yateto::leadDim<seissol::nodal::init::nodes2D>();
  cl::sycl::nd_range rng{{numElements * workGroupSize}, {workGroupSize}};

  queue->parallel_for(rng, [=](cl::sycl::nd_item<1> item) {
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      constexpr int pIdx = 0;
      constexpr int uIdx = 6;
      constexpr auto num2dNodes = seissol::nodal::tensor::nodes2D::Shape[0];

      const int tid = item.get_local_id(0);
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
  });
}
} // namespace seissol::kernels::time::aux
