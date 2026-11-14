// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/SderivLibrary.h"

#include <map>
#include <string>

namespace seissol::expr {

namespace {

// The stress of easi's special components: their formulas, with the operations of easi's C++ in
// its order, so that a program computes easi's bits (as the easi walker's translation of the
// components does, Reader/Scripting/EasiSpecial.cpp).
constexpr const char* Stress = R"sderiv(
# The stress of a fault optimally oriented under friction, normalized as easi's special components
# give it (easi's Special.cpp, OptimalStress.cpp and AndersonianStress.cpp), to easi's bits.

type normalized = {b_xx, b_yy, b_zz, b_xy, b_yz, b_xz}

# Aochi and Madariaga (2003, 2004), extended to dip slip faults: easi's !STRESS_STR_DIP_SLIP_AM.
def str_dip_slip_am(p : {mu_d, mu_s, strike, dip, s_zz, cohesion, R, DipSlipFaulting, s2ratio})
    : normalized =
  let phi = pi / 4 - 0.5 * atan(p.mu_s) in
  let strike = p.strike * pi / 180 in
  let dip = p.dip * pi / 180 in
  let phi_xyz = phi - dip in
  let s2 = sin(2 * phi) in
  let c2 = cos(2 * phi) in
  let c2xyz = cos(2 * phi_xyz) in
  let c2bis = c2 - c2xyz in
  let szz_inv = 1 / p.s_zz in
  select(le(abs(p.DipSlipFaulting), 0x1p-52),
         _am_strike_slip(p, phi, strike, dip, s2, c2, szz_inv),
         _am_dip_slip(p, phi_xyz, strike, s2, c2bis, c2xyz, szz_inv))

def _am_strike_slip(p, phi, strike, dip, s2, c2, szz_inv) =
  let ds = (p.mu_d * p.s_zz + p.R * (p.cohesion + (p.mu_s - p.mu_d) * p.s_zz)) /
           (s2 + p.mu_d * c2 + p.R * (p.mu_s - p.mu_d) * c2) in
  let sm = p.s_zz in
  let s11 = sm + ds in
  let s22 = sm - ds + 2 * ds * p.s2ratio in
  let s33 = sm - ds in
  let cd = cos(phi) in
  let sd = sin(phi) in
  let ci = cos(dip) in
  let si = sin(dip) in
  let cs = cos(strike) in
  let ss = sin(strike) in
  {b_xx = szz_inv * (cd * cd * cs * cs * s33 * si * si + cs * cs * s11 * sd * sd * si * si +
                     2 * cd * cs * s11 * sd * si * ss - 2 * cd * cs * s33 * sd * si * ss +
                     cd * cd * s11 * ss * ss + ci * ci * cs * cs * s22 + s33 * sd * sd * ss * ss),
   b_xy = szz_inv * (-cd * cd * cs * s33 * si * si * ss - cs * s11 * sd * sd * si * si * ss +
                     cd * cs * cs * s11 * sd * si - cd * cs * cs * s33 * sd * si -
                     cd * s11 * sd * si * ss * ss + cd * s33 * sd * si * ss * ss +
                     cd * cd * cs * s11 * ss - ci * ci * cs * s22 * ss + cs * s33 * sd * sd * ss),
   b_xz = szz_inv * (cd * cd * ci * cs * s33 * si + ci * cs * s11 * sd * sd * si +
                     cd * ci * s11 * sd * ss - cd * ci * s33 * sd * ss - cs * ci * s22 * si),
   b_yy = szz_inv * (cd * cd * s33 * si * si * ss * ss + s11 * sd * sd * si * si * ss * ss -
                     2 * cd * cs * s11 * sd * si * ss + 2 * cd * cs * s33 * sd * si * ss +
                     cd * cd * cs * cs * s11 + ci * ci * s22 * ss * ss + cs * cs * s33 * sd * sd),
   b_yz = szz_inv * (-cd * cd * ci * s33 * si * ss - ci * s11 * sd * sd * si * ss +
                     cd * ci * cs * s11 * sd - cd * ci * cs * s33 * sd + ss * ci * s22 * si),
   b_zz = szz_inv * (cd * cd * ci * ci * s33 + ci * ci * s11 * sd * sd + si * si * s22)}

# ds from R (A5, Aochi and Madariaga 2003), sigma 1 and 3 in the yz plane, s_zz and the mean stress
# related by a rotation by phi + dip (A3)
def _am_dip_slip(p, phi_xyz, strike, s2, c2bis, c2xyz, szz_inv) =
  let ds = (p.mu_d * p.s_zz + p.R * (p.cohesion + (p.mu_s - p.mu_d) * p.s_zz)) /
           (s2 + p.mu_d * c2bis + p.R * (p.mu_s - p.mu_d) * c2bis) in
  let sm = p.s_zz + ds * c2xyz in
  let s11 = sm + ds in
  let s22 = sm - ds + 2 * ds * p.s2ratio in
  let s33 = sm - ds in
  let cd = cos(phi_xyz) in
  let sd = sin(phi_xyz) in
  let cs = cos(strike) in
  let ss = sin(strike) in
  {b_xx = szz_inv * (cd * cd * cs * cs * s11 + cs * cs * s33 * sd * sd + ss * ss * s22),
   b_yy = szz_inv * (cd * cd * s11 * ss * ss + s33 * sd * sd * ss * ss + cs * cs * s22),
   b_zz = szz_inv * (cd * cd * s33 + sd * sd * s11),
   b_xy = szz_inv * (-cd * cd * cs * s11 * ss - cs * s33 * sd * sd * ss + ss * s22 * cs),
   b_yz = szz_inv * (-ss * cd * s11 * sd + ss * sd * s33 * cd),
   b_xz = szz_inv * (cs * cd * s11 * sd - cs * sd * s33 * cd)}

# Ulrich et al. (2018): easi's !OptimalStress.
def optimal(p : {mu_d, mu_s, strike, dip, rake, effectiveConfiningStress, cohesion, R, s2ratio})
    : normalized =
  let phi = pi / 4 - 0.5 * atan(p.mu_s) in
  let strike = p.strike * pi / 180 in
  let dip = p.dip * pi / 180 in
  let rake = p.rake * pi / 180 in
  let s2 = sin(2 * phi) in
  let c2 = cos(2 * phi) in
  let alpha = (2 * p.s2ratio - 1) / 3 in
  let confining = abs(p.effectiveConfiningStress) in
  let ds = (p.mu_d * confining + p.R * (p.cohesion + (p.mu_s - p.mu_d) * confining)) /
           (s2 + p.mu_d * (alpha + c2) + p.R * (p.mu_s - p.mu_d) * (alpha + c2)) in
  let sm = confining - alpha * ds in
  let s11 = sm + ds in
  let s22 = sm - ds + 2 * ds * p.s2ratio in
  let s33 = sm - ds in
  let cd = cos(-phi) in
  let sd = sin(-phi) in
  let cr = cos(rake) in
  let sr = sin(rake) in
  let ci = cos(dip) in
  let si = sin(dip) in
  let cs = cos(strike) in
  let ss = sin(strike) in
  {b_xx = -((ci * ci * s11 * sr * sr + s33 * si * si) * cd * cd +
            2 * ci * sd * si * sr * (s33 - s11) * cd +
            (s33 * sd * sd * sr * sr + cr * cr * s22) * ci * ci + sd * sd * si * si * s11) *
              cs * cs +
          2 * cr * cs *
              (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd -
               ci * sr * (-s33 * sd * sd + s22)) *
              ss -
          ss * ss * (cd * cd * cr * cr * s11 + cr * cr * s33 * sd * sd + s22 * sr * sr),
   b_xy = cr *
              (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd -
               ci * sr * (-s33 * sd * sd + s22)) *
              cs * cs +
          ss *
              ((-cd * cd * s11 + ci * ci * s22 - s33 * sd * sd) * cr * cr +
               (ci * ci * s11 * sr * sr + s33 * si * si) * cd * cd +
               2 * ci * sd * si * sr * (s33 - s11) * cd + s33 * ci * ci * sd * sd * sr * sr +
               sd * sd * si * si * s11 - s22 * sr * sr) *
              cs -
          ss * ss * cr *
              (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd -
               ci * sr * (-s33 * sd * sd + s22)),
   b_xz = cd * cs * sd * sr * (s33 - s11) * si * si -
          ((cs * (-s11 * sr * sr + s33) * ci + cr * sr * ss * s11) * cd * cd -
           cs * ((s33 * sr * sr - s11) * sd * sd + cr * cr * s22) * ci -
           cr * sr * ss * (-s33 * sd * sd + s22)) *
              si -
          cd * ci * sd * (ci * cs * sr - cr * ss) * (s33 - s11),
   b_yy = -((ci * ci * s11 * sr * sr + s33 * si * si) * cd * cd +
            2 * ci * sd * si * sr * (s33 - s11) * cd +
            (s33 * sd * sd * sr * sr + cr * cr * s22) * ci * ci + sd * sd * si * si * s11) *
              ss * ss -
          2 * cr * cs *
              (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd -
               ci * sr * (-s33 * sd * sd + s22)) *
              ss -
          cs * cs * (cd * cd * cr * cr * s11 + cr * cr * s33 * sd * sd + s22 * sr * sr),
   b_yz = -cd * sd * sr * ss * (s33 - s11) * si * si -
          ((-ss * (-s11 * sr * sr + s33) * ci + cr * cs * sr * s11) * cd * cd +
           ss * ((s33 * sr * sr - s11) * sd * sd + cr * cr * s22) * ci -
           cr * cs * sr * (-s33 * sd * sd + s22)) *
              si +
          cd * ci * sd * (ci * sr * ss + cr * cs) * (s33 - s11),
   b_zz = -(cd * cd * s11 * sr * sr + s33 * sd * sd * sr * sr + cr * cr * s22) * si * si -
          2 * cd * ci * sd * sr * (s11 - s33) * si -
          ci * ci * (cd * cd * s33 + s11 * sd * sd)}

# Ulrich et al. (2018), the vertical principal stress S_v the first (1), second (2) or third:
# easi's !AndersonianStress. Compressive stress is negative.
def andersonian(p : {mu_d, mu_s, SH_max, sig_zz, cohesion, S, s2ratio, S_v}) : normalized =
  let phi = pi / 4 - 0.5 * atan(p.mu_s) in
  let sh_max = p.SH_max * pi / 180 in
  let s2 = sin(2 * phi) in
  let c2 = cos(2 * phi) in
  let alpha = (2 * p.s2ratio - 1) / 3 in
  let confining = abs(p.sig_zz) in
  let r = 1 / (1 + p.S) in
  let ds = (p.mu_d * confining + r * (p.cohesion + (p.mu_s - p.mu_d) * confining)) /
           (s2 + p.mu_d * (alpha + c2) + r * (p.mu_s - p.mu_d) * (alpha + c2)) in
  let sm = confining - alpha * ds in
  let s11 = sm + ds in
  let s22 = sm - ds + 2 * ds * p.s2ratio in
  let s33 = sm - ds in
  let vertical1 = eq(p.S_v, 1) in
  let vertical2 = eq(p.S_v, 2) in
  let sh = select(vertical1, s22, s11) in
  let shmin = select(vertical1, s33, select(vertical2, s33, s22)) in
  let sv = select(vertical1, s11, select(vertical2, s22, s33)) in
  let cs = cos(sh_max) in
  let ss = sin(sh_max) in
  {b_xx = -(cs * cs * shmin + ss * ss * sh) * confining / sv,
   b_xy = -cs * ss * (sh - shmin) * confining / sv,
   b_xz = 0,
   b_yy = -(ss * ss * shmin + cs * cs * sh) * confining / sv,
   b_yz = 0,
   b_zz = -confining}
)sderiv";

constexpr const char* Rock = R"sderiv(
# Rocks.

# The elastic material of a poroelastic one, the fluid locked in its pores: Gassmann's undrained
# moduli (Biot), lambda + alpha^2 M with alpha = 1 - K / K_s and 1 / M = phi / K_f + (alpha - phi)
# / K_s, K the bulk modulus of the drained matrix, and the density of solid and fluid together.
def gassmann(p : {bulk_solid, bulk_fluid, porosity, rho, rho_fluid, mu, lambda}) : elastic =
  let bulk_matrix = p.lambda + 2 / 3 * p.mu in
  let a = (1 - bulk_matrix / p.bulk_solid) * (1 - bulk_matrix / p.bulk_solid) in
  let b = p.porosity / p.bulk_fluid in
  let c = (1 - p.porosity) / p.bulk_solid in
  let d = bulk_matrix / (p.bulk_solid * p.bulk_solid) in
  {rho = p.porosity * p.rho_fluid + (1 - p.porosity) * p.rho,
   mu = p.mu,
   lambda = p.lambda + a / (b + c - d)}
)sderiv";

constexpr const char* Rotation = R"sderiv(
# Rotations by strike, dip and rake, in degrees.

# The rotation of SeisSol's point sources (SourceTerm/PointSource.cpp), r[i][j] its row i and
# column j: a moment tensor m of the fault is transpose(r) @ m @ r in x, y and z.
def matrix(strike, dip, rake) =
  let s = strike * pi / 180 in
  let d = dip * pi / 180 in
  let r = rake * pi / 180 in
  let cs = cos(s) in
  let ss = sin(s) in
  let cd = cos(d) in
  let sd = sin(d) in
  let cr = cos(r) in
  let sr = sin(r) in
  [[cr * cs + cd * sr * ss, cd * cs * sr - cr * ss, -sd * sr],
   [cd * cr * ss - cs * sr, sr * ss + cd * cr * cs, -cr * sd],
   [sd * ss, cs * sd, cd]]

# The tensor m of the fault in x, y and z, by the rotation r of matrix().
def rotate(m, r) = transpose(r) @ m @ r
)sderiv";

constexpr const char* Motion = R"sderiv(
# Ground motion.

# n angles from 0, 180 / n degrees apart, below 180.
def angles(n) = [k * 180 / n for k in range(n)]

# The horizontal motion (a1, a2) along each of the angles, in degrees from the x axis:
# |a1 cos(t) + a2 sin(t)|. Its largest value over time per angle, kept in a state, and the
# largest of these give RotD100.
def rotated(a1, a2, angles) =
  [abs(a1 * cos(t * pi / 180) + a2 * sin(t * pi / 180)) for t in angles]

# Of the largest motions along the angles (`peaks`, of angles() as above): RotD50, their median;
# and GMRotD50 (Boore, Watson-Lamprey and Abrahamson 2006), the median over the angles below 90
# degrees of the geometric mean of the largest motions along an angle and along the angle 90
# degrees on, of an even number of angles.
def rotd50(peaks) = median(peaks)
def gmrotd50(peaks) =
  median([sqrt(peaks[k] * peaks[k + len(peaks) / 2]) for k in range(len(peaks) / 2)])
)sderiv";

} // namespace

const std::map<std::string, std::string>& libraryModules() {
  static const std::map<std::string, std::string> Modules = {
      {"stress", Stress}, {"rock", Rock}, {"rotation", Rotation}, {"motion", Motion}};
  return Modules;
}

} // namespace seissol::expr
