#ifndef SEISSOL_BASICTYPEDEFS_HPP
#define SEISSOL_BASICTYPEDEFS_HPP

#include <Kernels/precision.hpp>

enum mpiTag { localIntegrationData = 0, neighboringIntegrationData = 1, timeData = 2 };

enum TimeClustering {
  // global time stepping
  single = 0,
  // offline clustering computed in pre-processing
  offline = 1,
  // online clustering resulting in a multi-rate scheme
  multiRate = 2,
  // online clustering aiming at LTS for slithers only
  slithers = 3
};

// face types
// Note: When introducting new types also change
// int seissol::initializers::time_stepping::LtsWeights::getBoundaryCondition
// and PUMLReader. Otherwise it might become a DR face...
enum class FaceType {
  // regular: inside the computational domain
  regular = 0,

  // free surface boundary
  freeSurface = 1,

  // free surface boundary with gravity
  freeSurfaceGravity = 2,

  // dynamic rupture boundary
  dynamicRupture = 3,

  // Dirichlet boundary
  dirichlet = 4,

  // absorbing/outflow boundary
  outflow = 5,

  // periodic boundary
  periodic = 6,

  // analytical boundary (from initial cond.)
  analytical = 7
};

enum SystemType { Host = 0, Device = 1 };

// plasticity information per cell
struct PlasticityData {
  // initial loading (stress tensor)
  real initialLoading[6];
  real cohesionTimesCosAngularFriction;
  real sinAngularFriction;
  real mufactor;
};

#endif // SEISSOL_BASICTYPEDEFS_HPP
