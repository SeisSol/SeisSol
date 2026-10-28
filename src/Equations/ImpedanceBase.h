// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_IMPEDANCEBASE_H_
#define SEISSOL_SRC_EQUATIONS_IMPEDANCEBASE_H_

#include <Eigen/Dense>
#include <cstddef>

namespace seissol::model {

/// Admittance of one side of a face, mapping the interface traction to the interface velocity.
template <std::size_t Dim>
using AdmittanceMatrix = Eigen::Matrix<double, Dim, Dim>;

/// Maps the interface traction to the three stress components which do not take part in the
/// fault-normal Riemann problem (sigma_ss, sigma_dd, sigma_sd).
template <std::size_t Dim>
using LateralStressMatrix = Eigen::Matrix<double, 3, Dim>;

/**
 * Wave admittance of a half space of the given material, for a face whose normal is the local x
 * axis. Specialized in Equations/<equation>/Model/Impedance.h for every material that supports
 * dynamic rupture; a specialization provides
 *
 *   TractionIndices, VelocityIndices  the quantities of the face-aligned frame that make up the
 *                                     interface traction and velocity
 *   Dim                               their number: 3, or 4 where the fluid couples in
 *   Matrix, LateralMatrix             AdmittanceMatrix<Dim>, LateralStressMatrix<Dim>
 *   signature()                       sign of each stored traction component relative to its
 *                                     energy conjugate one; the identity where they agree
 *   admittance(materialLocal, lateralStress = nullptr)
 *
 * `materialLocal` has to be rotated into the fault-local frame already (Bond matrix). If
 * `lateralStress` is given, it receives the reconstruction of the stress components outside the
 * Riemann problem from the traction difference across the face.
 *
 * The specializations only depend on the material parameters -- neither on the MaterialT of the
 * build nor on generated tensors -- so all of them compile in every build. The anisotropic and
 * poroelastic closed forms are checked against the eigendecomposition of the normal Jacobian in
 * tests/Model (the poroelastic one in a poroelastic build); the elastic one, which the viscoelastic
 * material reuses, against the isotropic limit of the anisotropic one.
 *
 * NOTE ON NAMING: what the `Zplus`/`Zminus` tensors and `dr::ImpedanceMatrices::impedance{,Neig}`
 * store is *not* an impedance but an admittance Y, mapping traction to velocity. The quantity
 * carrying the dimension of an impedance is `eta` (see
 * Initializer/Model/DynamicRuptureImpedance.h).
 */
template <typename MaterialT>
struct ImpedanceCompute;

/// Admittance of one side; `materialLocal` in the fault-local frame. See ImpedanceCompute.
template <typename MaterialT>
typename ImpedanceCompute<MaterialT>::Matrix computeAdmittance(
    const MaterialT& materialLocal,
    typename ImpedanceCompute<MaterialT>::LateralMatrix* lateralStress = nullptr) {
  return ImpedanceCompute<MaterialT>::admittance(materialLocal, lateralStress);
}

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_IMPEDANCEBASE_H_
