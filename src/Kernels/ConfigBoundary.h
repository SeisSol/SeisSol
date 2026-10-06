// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_CONFIGBOUNDARY_H_
#define SEISSOL_SRC_KERNELS_CONFIGBOUNDARY_H_

#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/configboundary.h"
#include "Initializer/CellLocalInformation.h"
#include "Model/Common.h"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::kernels {

/**
 * @brief Whether the time integral of a face neighbor of the configuration `NeighborCfg` converts
 * into the configuration `Cfg` of the cell.
 *
 * The neighbor flux of the cell reads the time integral of a neighbor only through its trace on
 * the shared face, tested with the face basis of the cell. The conversion keeps these face
 * integrals exact. It goes through the canonical form of the family of both configurations (the
 * configurations whose materials pose the Riemann problem in the same material, and that fuse the
 * same number of simulations): toCanonical of the neighbor pads its time integral to the largest
 * order of the family in the build, selects the quantities of the Riemann problem and widens it to
 * double precision; fromCanonical of the cell, for the side of the neighbor on the shared face,
 * projects the trace on that face to the face basis of the cell and lifts it into the volume basis
 * of the cell, and narrows it to the precision of the cell. Quantities of the cell outside the
 * canonical form, such as memory variables, are zero. The code generator generates both kernels
 * for the configurations of a family that has more than one configuration in the build.
 */
template <typename Cfg, typename NeighborCfg>
constexpr bool Convertible =
    generated::ConfigBoundaryKernels<Cfg>::Host &&
    generated::ConfigBoundaryKernels<NeighborCfg>::Host &&
    model::CanNeighbor<model::MaterialOf<Cfg>, model::MaterialOf<NeighborCfg>> &&
    Cfg::NumSimulations == NeighborCfg::NumSimulations;

/**
 * @brief The faces of the cells of the configuration `Cfg` whose neighbor computes in another
 * configuration.
 *
 * The time integrals of such a neighbor are in its own configuration. They are computed there,
 * from the buffers or derivatives of the neighbor and the time basis of its solver, and converted
 * into the configuration of the cell.
 */
template <typename Cfg>
class ConfigBoundary {
  public:
  using real = Real<Cfg>;

  ConfigBoundary() = default;

  /// The boundaries towards the cells of all other configurations in `configs`.
  explicit ConfigBoundary(const std::vector<ConfigId>& configs);

  /// Whether a cell of the configuration `Cfg` can have a neighbor of another configuration.
  [[nodiscard]] bool empty() const { return empty_; }

  /**
   * Sets the time intervals of the next neighbor integration: the time step `timestep` of the
   * cluster, and for a neighbor with a larger time step, the part of its time step from
   * `subTimeStart` on, with `neighborTimestep` being its whole time step.
   */
  void setIntervals(double timestep, double subTimeStart, double neighborTimestep);

  /**
   * For the faces of a cell with `cellInformation` whose neighbor computes in another
   * configuration: computes the time integral of the neighbor from `timeDofs` into
   * `integrationBuffer`, and points `timeIntegrated` to it. The other faces are left as they are.
   */
  void computeIntegrals(const CellLocalInformation& cellInformation,
                        const std::array<void*, Cell::NumFaces>& timeDofs,
                        const std::array<real*, Cell::NumFaces>& integrationBuffer,
                        std::array<real*, Cell::NumFaces>& timeIntegrated) const;

  private:
  struct Neighbor {
    // the coefficients of the time basis of the neighbor for its step relation to the cell (GTS)
    // and for the sub-interval of a larger time step of the neighbor (LTS)
    std::vector<double> timeCoeffs;
    std::vector<double> subtimeCoeffs;
  };

  /// Computes the time integral of the neighbor of the configuration `NeighborCfg` on the face
  /// `face` from `timeDofs`, and converts it into `integrationBuffer`.
  template <typename NeighborCfg>
  static void computeIntegral(const CellLocalInformation& cellInformation,
                              std::size_t face,
                              const Neighbor& neighbor,
                              const void* timeDofs,
                              real* integrationBuffer);

  // indexed by the configuration of the neighbor
  std::vector<std::optional<Neighbor>> neighbors_;
  bool empty_{true};
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_CONFIGBOUNDARY_H_
