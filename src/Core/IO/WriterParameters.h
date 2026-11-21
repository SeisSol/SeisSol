// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_IO_WRITERPARAMETERS_H_
#define SEISSOL_SRC_CORE_IO_WRITERPARAMETERS_H_

namespace seissol::initializer::parameters {

// an interval longer than any run
constexpr double VeryLongTime = 1.0e100;

enum class XdmfBackend : int { Posix, Hdf5 };

/**
 * @brief How the steps of a mesh output are distributed over files.
 *
 * `Snapshot` writes one file per step, each complete in itself. `Incremental` keeps the geometry
 * -- which does not move -- in a file of its own and lets every snapshot link to it, which is the
 * cheaper choice whenever the mesh is large against the fields written on it. `Monolith` puts
 * every step into a single file as a VTKHDF time series, so a run leaves one file behind instead
 * of one per output step.
 *
 * Only the VTKHDF format can do more than `Snapshot`; an Xdmf output (that is, one written
 * without a `vtkorder`) references its unchanging parts from every step anyway. An output whose
 * parameters choose none is written as `Monolith` in VTKHDF and as a snapshot per step in Xdmf.
 */
enum class TimeSeriesMode : int { Snapshot, Incremental, Monolith };

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_CORE_IO_WRITERPARAMETERS_H_
