#include "GravitationalFreeSurfaceBC.h"
#include "SeisSol.h"

namespace seissol {

double getGravitationalAcceleration() {
  return SeisSol::main.getGravitationSetup().acceleration;
}

std::pair<long long, long long>
GravitationalFreeSurfaceBc::getFlopsDisplacementFace(unsigned int face, FaceType faceType) {
  long long hardwareFlops = 0;
  long long nonZeroFlops = 0;

  constexpr auto numberOfNodes = nodal::tensor::nodes2D::Shape[0];

  // initialize integral of displacement
  hardwareFlops += 1 * numberOfNodes;
  nonZeroFlops += 1 * numberOfNodes;

  // Note: This neglects roughly 10 * CONVERGENCE_ORDER * numNodes2D flops
  // Before adjusting the range of the loop, check range of loop in computation!
  for (int order = 1; order < CONVERGENCE_ORDER + 1; ++order) {
#ifdef USE_ELASTIC
    constexpr auto flopsPerQuadpoint =
        4 + // Computing coefficient
        6 + // Updating displacement
        2; // Updating integral of displacement
#else
    constexpr auto flopsPerQuadpoint = 0;
#endif
    constexpr auto flopsUpdates = flopsPerQuadpoint * numberOfNodes;

    nonZeroFlops += kernel::projectDerivativeToNodalBoundaryRotated::nonZeroFlops(order - 1, face) + flopsUpdates;
    hardwareFlops += kernel::projectDerivativeToNodalBoundaryRotated::hardwareFlops(order - 1, face) + flopsUpdates;
  }

  // Two rotations: One to face-aligned, one to global
  hardwareFlops += 2 * kernel::rotateFaceDisplacement::HardwareFlops;
  nonZeroFlops += 2 * kernel::rotateFaceDisplacement::NonZeroFlops;

  return {nonZeroFlops, hardwareFlops};
}
} // namespace seissol
