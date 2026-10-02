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
#include "Initializer/CellLocalInformation.h"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::kernels {

/**
 * @brief Converts the time integral of a face neighbor that computes in the configuration
 * `neighbor` into the configuration `cell` of the cell.
 *
 * The neighbor kernel of the cell reads the time integral of a neighbor only through its trace on
 * the shared face, tested with the basis functions of the cell. The conversion keeps these face
 * integrals exact. The modal bases are hierarchical, and orthogonal on the volume and on the
 * faces. So a neighbor of at most the order of the cell is copied, padded with zeros. For a
 * neighbor of a higher order, its trace on the face is projected to the face basis of the cell and
 * lifted back into the volume basis of the cell; the result agrees with the neighbor on that face
 * only.
 *
 * The quantities are matched by name. Quantities of the cell that the neighbor does not have, such
 * as the memory variables of an anelastic cell next to an elastic one, are zero; the neighbor flux
 * reads the elastic quantities only. The fused simulations are matched one by one.
 */
class NeighborConversion {
  public:
  NeighborConversion(ConfigId cell, ConfigId neighbor);

  /**
   * Writes `integral`, a time integral of a neighbor of the configuration `NeighborCfg`, into
   * `converted`, a time integral of the configuration `Cfg` of the cell. The neighbor touches the
   * cell with its side `neighborSide`.
   */
  template <typename Cfg, typename NeighborCfg>
  void apply(const Real<NeighborCfg>* integral,
             Real<Cfg>* converted,
             std::size_t neighborSide) const;

  private:
  // for every quantity of a time integral of the cell, the quantity of the neighbor with its name
  std::vector<std::optional<std::size_t>> quantities_;

  // for a neighbor of a higher order: per side of the neighbor, the map from its volume basis to
  // the volume basis of the cell (row-major); empty otherwise
  std::array<std::vector<double>, Cell::NumFaces> lift_;
};

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
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

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
    NeighborConversion conversion;

    // the coefficients of the time basis of the neighbor for its step relation to the cell (GTS)
    // and for the sub-interval of a larger time step of the neighbor (LTS)
    std::vector<double> timeCoeffs;
    std::vector<double> subtimeCoeffs;
  };

  // indexed by the configuration of the neighbor
  std::vector<std::optional<Neighbor>> neighbors_;
  bool empty_{true};
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_CONFIGBOUNDARY_H_
