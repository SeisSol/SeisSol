#ifndef DR_TYPEDEFS
#define DR_TYPEDEFS

#include "Kernels/precision.hpp"

namespace seissol::dr {
constexpr unsigned int numberOfBoundaryGaussPoints =
    (CONVERGENCE_ORDER + 1) * (CONVERGENCE_ORDER + 1);

/**
 * Stores the different types of friction laws
 * The values resemble the identifiers used in the old fortran implementation
 */
enum class FrictionLawType : unsigned int {
  no_fault = 0,
  linear_slip_weakening = 2,
  linear_slip_weakening_bimaterial = 6,
  linear_slip_weakening_forced_rupture_time = 16,
  rate_and_state_aging_law = 3,
  rate_and_state_slip_law = 4,
  rate_and_state_fast_velocity_weakening = 103,
  imposed_slip_rates = 33,
  rate_and_state_velocity_weakening = 7,
  rate_and_state_aging_nucleation = 101,
};

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Uphoff equation (4.51)
 */
struct ImpedancesAndEta {
  real Zp, Zs, Zp_neig, Zs_neig, eta_p, eta_s, inv_eta_s, inv_Zp, inv_Zs, inv_Zp_neig, inv_Zs_neig;
};

constexpr unsigned int numberOfTPGridPoints = 60;
} // namespace seissol::dr

#endif
