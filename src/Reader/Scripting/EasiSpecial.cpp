// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "EasiSpecial.h"

#include "Expr/Ir.h"

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

namespace {

using expr::Fn;
using expr::NodeId;

/// A value of a formula: a node, with the arithmetic of C++ on it, in the order C++ evaluates an
/// expression -- computed where the operands are numbers, with the same operation.
class E {
  public:
  E(expr::Arena& arena, NodeId node) : arena_(&arena), node_(node) {}

  [[nodiscard]] NodeId node() const { return node_; }
  [[nodiscard]] std::optional<double> number() const {
    const auto& node = (*arena_)[node_];
    if (node.kind == expr::Kind::Const) {
      return node.value;
    }
    return std::nullopt;
  }
  [[nodiscard]] E of(double value) const { return E(*arena_, arena_->konst(value)); }

  friend E operator+(const E& a, const E& b) { return a.binary(Fn::Add, b, std::plus<>()); }
  friend E operator-(const E& a, const E& b) { return a.binary(Fn::Sub, b, std::minus<>()); }
  friend E operator*(const E& a, const E& b) { return a.binary(Fn::Mul, b, std::multiplies<>()); }
  friend E operator/(const E& a, const E& b) { return a.binary(Fn::Div, b, std::divides<>()); }
  friend E operator+(const E& a, double b) { return a + a.of(b); }
  friend E operator-(const E& a, double b) { return a - a.of(b); }
  friend E operator*(const E& a, double b) { return a * a.of(b); }
  friend E operator/(const E& a, double b) { return a / a.of(b); }
  friend E operator+(double a, const E& b) { return b.of(a) + b; }
  friend E operator-(double a, const E& b) { return b.of(a) - b; }
  friend E operator*(double a, const E& b) { return b.of(a) * b; }
  friend E operator/(double a, const E& b) { return b.of(a) / b; }
  friend E operator-(const E& a) { return a.unary(Fn::Neg, std::negate<>()); }
  friend E sin(const E& a) {
    return a.unary(Fn::Sin, [](double x) { return std::sin(x); });
  }
  friend E cos(const E& a) {
    return a.unary(Fn::Cos, [](double x) { return std::cos(x); });
  }
  friend E atan(const E& a) {
    return a.unary(Fn::Atan, [](double x) { return std::atan(x); });
  }
  friend E fabs(const E& a) {
    return a.unary(Fn::Abs, [](double x) { return std::fabs(x); });
  }
  /// 1 where `a <= b`, else 0
  friend E operator<=(const E& a, double b) {
    return a.binary(Fn::Le, a.of(b), [](double x, double y) { return x <= y ? 1.0 : 0.0; });
  }
  /// 1 where `a == b`, else 0
  friend E operator==(const E& a, double b) {
    return a.binary(Fn::Eq, a.of(b), [](double x, double y) { return x == y ? 1.0 : 0.0; });
  }
  /// `condition ? a : b`
  friend E select(const E& condition, const E& a, const E& b) {
    if (const auto known = condition.number(); known.has_value()) {
      return *known != 0.0 ? a : b;
    }
    return E(*a.arena_, a.arena_->pw(Fn::Select, condition.node_, a.node_, b.node_));
  }

  private:
  template <typename F>
  [[nodiscard]] E unary(Fn fn, F f) const {
    if (const auto x = number(); x.has_value()) {
      return of(f(*x));
    }
    return E(*arena_, arena_->pw(fn, node_));
  }
  template <typename F>
  [[nodiscard]] E binary(Fn fn, const E& other, F f) const {
    const auto x = number();
    const auto y = other.number();
    if (x.has_value() && y.has_value()) {
      return of(f(*x, *y));
    }
    return E(*arena_, arena_->pw(fn, node_, other.node_));
  }

  expr::Arena* arena_;
  NodeId node_;
};

/// The outputs of every special component: a stress, normalized as easi has it.
struct Stress {
  E bXx;
  E bYy;
  E bZz;
  E bXy;
  E bYz;
  E bXz;
};

/// `condition ? then() : otherwise()`, of every component: the one branch where the condition is
/// a number.
template <typename Then, typename Otherwise>
Stress choose(const E& condition, Then then, Otherwise otherwise) {
  if (const auto known = condition.number(); known.has_value()) {
    return *known != 0.0 ? then() : otherwise();
  }
  const Stress a = then();
  const Stress b = otherwise();
  return Stress{select(condition, a.bXx, b.bXx),
                select(condition, a.bYy, b.bYy),
                select(condition, a.bZz, b.bZz),
                select(condition, a.bXy, b.bXy),
                select(condition, a.bYz, b.bYz),
                select(condition, a.bXz, b.bXz)};
}

constexpr double Pi = 3.141592653589793;

// ------------------------------------------------------ STRESS_STR_DIP_SLIP_AM ---

struct StressStrDipSlipAmIn {
  E muD;
  E muS;
  E strike;
  E dip;
  E sZz;
  E cohesion;
  E r;
  E dipSlipFaulting;
  E s2ratio;
};

// easi's STRESS_STR_DIP_SLIP_AM::evaluate(): the normalized stress following the method of Aochi
// and Madariaga (2004), extended to dip slip faults
Stress stressStrDipSlipAm(const StressStrDipSlipAmIn& i) {
  // most favorable direction (A4, AM2003)
  const E phi = Pi / 4.0 - 0.5 * atan(i.muS);
  const E strikeRad = i.strike * Pi / 180.0;
  const E dipRad = i.dip * Pi / 180.0;
  const E phiXyz = phi - dipRad;
  const E s2 = sin(2.0 * phi);
  const E c2 = cos(2.0 * phi);
  const E c2xyz = cos(2.0 * phiXyz);
  const E c2bis = c2 - c2xyz;
  const E szzInv = 1.0 / i.sZz;

  const auto strikeSlip = [&]() {
    const E ds = (i.muD * i.sZz + i.r * (i.cohesion + (i.muS - i.muD) * i.sZz)) /
                 (s2 + i.muD * c2 + i.r * (i.muS - i.muD) * c2);
    const E sm = i.sZz;
    const E s11 = sm + ds;
    const E s22 = sm - ds + 2.0 * ds * i.s2ratio;
    const E s33 = sm - ds;

    const E cd = cos(phi);
    const E sd = sin(phi);

    const E ci = cos(dipRad);
    const E si = sin(dipRad);

    const E cs = cos(strikeRad);
    const E ss = sin(strikeRad);

    const E bXx =
        szzInv * (cd * cd * cs * cs * s33 * si * si + cs * cs * s11 * sd * sd * si * si +
                  2 * cd * cs * s11 * sd * si * ss - 2 * cd * cs * s33 * sd * si * ss +
                  cd * cd * s11 * ss * ss + ci * ci * cs * cs * s22 + s33 * sd * sd * ss * ss);
    const E bXy =
        szzInv * (-cd * cd * cs * s33 * si * si * ss - cs * s11 * sd * sd * si * si * ss +
                  cd * cs * cs * s11 * sd * si - cd * cs * cs * s33 * sd * si -
                  cd * s11 * sd * si * ss * ss + cd * s33 * sd * si * ss * ss +
                  cd * cd * cs * s11 * ss - ci * ci * cs * s22 * ss + cs * s33 * sd * sd * ss);
    const E bXz = szzInv * (cd * cd * ci * cs * s33 * si + ci * cs * s11 * sd * sd * si +
                            cd * ci * s11 * sd * ss - cd * ci * s33 * sd * ss - cs * ci * s22 * si);
    const E bYy =
        szzInv * (cd * cd * s33 * si * si * ss * ss + s11 * sd * sd * si * si * ss * ss -
                  2 * cd * cs * s11 * sd * si * ss + 2 * cd * cs * s33 * sd * si * ss +
                  cd * cd * cs * cs * s11 + ci * ci * s22 * ss * ss + cs * cs * s33 * sd * sd);
    const E bYz = szzInv * (-cd * cd * ci * s33 * si * ss - ci * s11 * sd * sd * si * ss +
                            cd * ci * cs * s11 * sd - cd * ci * cs * s33 * sd + ss * ci * s22 * si);
    const E bZz = szzInv * (cd * cd * ci * ci * s33 + ci * ci * s11 * sd * sd + si * si * s22);
    return Stress{bXx, bYy, bZz, bXy, bYz, bXz};
  };

  const auto dipSlip = [&]() {
    // ds (delta_sigma) is deduced from R (A5, Aochi and Madariaga 2003), assuming that sig1 and
    // sig3 are in the yz plane; sigzz and sigma_ini are then related by a phi+dip rotation (using
    // A3, AM03): s_zz = sm - ds * cos(2.*(Phi+dip_rad)). Now we have to assume that P = sm (not
    // any more equal to s_zz), and we can obtain the new expression of ds:
    const E ds = (i.muD * i.sZz + i.r * (i.cohesion + (i.muS - i.muD) * i.sZz)) /
                 (s2 + i.muD * c2bis + i.r * (i.muS - i.muD) * c2bis);
    const E sm = i.sZz + ds * c2xyz;

    const E s11 = sm + ds;
    const E s22 = sm - ds + 2.0 * ds * i.s2ratio;
    const E s33 = sm - ds;
    const E cd = cos(phiXyz);
    const E sd = sin(phiXyz);
    const E cs = cos(strikeRad);
    const E ss = sin(strikeRad);

    const E bXx = szzInv * (cd * cd * cs * cs * s11 + cs * cs * s33 * sd * sd + ss * ss * s22);
    const E bYy = szzInv * (cd * cd * s11 * ss * ss + s33 * sd * sd * ss * ss + cs * cs * s22);
    const E bZz = szzInv * (cd * cd * s33 + sd * sd * s11);
    const E bXy = szzInv * (-cd * cd * cs * s11 * ss - cs * s33 * sd * sd * ss + ss * s22 * cs);
    const E bYz = szzInv * (-ss * cd * s11 * sd + ss * sd * s33 * cd);
    const E bXz = szzInv * (cs * cd * s11 * sd - cs * sd * s33 * cd);
    return Stress{bXx, bYy, bZz, bXy, bYz, bXz};
  };

  return choose(
      fabs(i.dipSlipFaulting) <= std::numeric_limits<double>::epsilon(), strikeSlip, dipSlip);
}

// --------------------------------------------------------------- OptimalStress ---

struct OptimalStressIn {
  E muD;
  E muS;
  E strike;
  E dip;
  E rake;
  E effectiveConfiningStress;
  E cohesion;
  E r;
  E s2ratio;
};

// easi's OptimalStress::evaluate(): the normalized stress following the method of Ulrich et al.
// (2018)
Stress optimalStress(const OptimalStressIn& i) {
  // most favorable direction (A4, AM2003)
  const E phi = Pi / 4.0 - 0.5 * atan(i.muS);
  const E strikeRad = i.strike * Pi / 180.0;
  const E dipRad = i.dip * Pi / 180.0;
  const E rakeRad = i.rake * Pi / 180.0;
  const E s2 = sin(2.0 * phi);
  const E c2 = cos(2.0 * phi);
  const E alpha = (2.0 * i.s2ratio - 1.0) / 3.0;
  const E effectiveConfiningStress = fabs(i.effectiveConfiningStress);

  const E ds = (i.muD * effectiveConfiningStress +
                i.r * (i.cohesion + (i.muS - i.muD) * effectiveConfiningStress)) /
               (s2 + i.muD * (alpha + c2) + i.r * (i.muS - i.muD) * (alpha + c2));
  const E sm = effectiveConfiningStress - alpha * ds;
  // sii are all positive
  const E s11 = sm + ds;
  const E s22 = sm - ds + 2.0 * ds * i.s2ratio;
  const E s33 = sm - ds;

  const E cd = cos(-phi);
  const E sd = sin(-phi);

  const E cr = cos(rakeRad);
  const E sr = sin(rakeRad);

  const E ci = cos(dipRad);
  const E si = sin(dipRad);

  const E cs = cos(strikeRad);
  const E ss = sin(strikeRad);

  const E bXx = -((ci * ci * s11 * sr * sr + s33 * si * si) * cd * cd +
                  2 * ci * sd * si * sr * (s33 - s11) * cd +
                  (s33 * sd * sd * sr * sr + cr * cr * s22) * ci * ci + sd * sd * si * si * s11) *
                    cs * cs +
                2 * cr * cs *
                    (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd -
                     ci * sr * (-s33 * sd * sd + s22)) *
                    ss -
                ss * ss * (cd * cd * cr * cr * s11 + cr * cr * s33 * sd * sd + s22 * sr * sr);
  const E bXy =
      cr *
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
          (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd - ci * sr * (-s33 * sd * sd + s22));
  const E bXz = cd * cs * sd * sr * (s33 - s11) * si * si -
                ((cs * (-s11 * sr * sr + s33) * ci + cr * sr * ss * s11) * cd * cd -
                 cs * ((s33 * sr * sr - s11) * sd * sd + cr * cr * s22) * ci -
                 cr * sr * ss * (-s33 * sd * sd + s22)) *
                    si -
                cd * ci * sd * (ci * cs * sr - cr * ss) * (s33 - s11);
  const E bYy = -((ci * ci * s11 * sr * sr + s33 * si * si) * cd * cd +
                  2 * ci * sd * si * sr * (s33 - s11) * cd +
                  (s33 * sd * sd * sr * sr + cr * cr * s22) * ci * ci + sd * sd * si * si * s11) *
                    ss * ss -
                2 * cr * cs *
                    (cd * cd * ci * sr * s11 + sd * si * (s33 - s11) * cd -
                     ci * sr * (-s33 * sd * sd + s22)) *
                    ss -
                cs * cs * (cd * cd * cr * cr * s11 + cr * cr * s33 * sd * sd + s22 * sr * sr);
  const E bYz = -cd * sd * sr * ss * (s33 - s11) * si * si -
                ((-ss * (-s11 * sr * sr + s33) * ci + cr * cs * sr * s11) * cd * cd +
                 ss * ((s33 * sr * sr - s11) * sd * sd + cr * cr * s22) * ci -
                 cr * cs * sr * (-s33 * sd * sd + s22)) *
                    si +
                cd * ci * sd * (ci * sr * ss + cr * cs) * (s33 - s11);
  const E bZz = -(cd * cd * s11 * sr * sr + s33 * sd * sd * sr * sr + cr * cr * s22) * si * si -
                2 * cd * ci * sd * sr * (s11 - s33) * si -
                ci * ci * (cd * cd * s33 + s11 * sd * sd);
  return Stress{bXx, bYy, bZz, bXy, bYz, bXz};
}

// ----------------------------------------------------------- AndersonianStress ---

struct AndersonianStressIn {
  E muD;
  E muS;
  E shMax;
  E sigZz;
  E cohesion;
  E s;
  E s2ratio;
  E sV;
};

// easi's AndersonianStress::evaluate(): the normalized stress following the method of Ulrich et
// al. (2018)
Stress andersonianStress(const AndersonianStressIn& i) {
  // most favorable direction (A4, AM2003)
  const E phi = Pi / 4.0 - 0.5 * atan(i.muS);
  const E shMaxRad = i.shMax * Pi / 180.0;
  const E s2 = sin(2.0 * phi);
  const E c2 = cos(2.0 * phi);
  const E alpha = (2.0 * i.s2ratio - 1.0) / 3.0;
  // we assume that szz is the effectiveConfiningStress  = sii/3 which is not necessarily true,
  // but then we scale the stress to get b_zz = szz
  const E effectiveConfiningStress = fabs(i.sigZz);
  const E r = 1.0 / (1.0 + i.s);
  const E ds = (i.muD * effectiveConfiningStress +
                r * (i.cohesion + (i.muS - i.muD) * effectiveConfiningStress)) /
               (s2 + i.muD * (alpha + c2) + r * (i.muS - i.muD) * (alpha + c2));
  const E sm = effectiveConfiningStress - alpha * ds;
  // sii are all positive
  const E s11 = sm + ds;
  const E s22 = sm - ds + 2.0 * ds * i.s2ratio;
  const E s33 = sm - ds;
  // S_v == 1: s22, s33, s11; S_v == 2: s11, s33, s22; else s11, s22, s33
  const E vertical1 = i.sV == 1;
  const E vertical2 = i.sV == 2;
  const E sHMax = select(vertical1, s22, s11);
  const E sHMin = select(vertical1, s33, select(vertical2, s33, s22));
  const E sVert = select(vertical1, s11, select(vertical2, s22, s33));

  const E cs = cos(shMaxRad);
  const E ss = sin(shMaxRad);

  // minus signs added to fall back to SeisSol convention (compressive stress<0)
  const E bXx = -(cs * cs * sHMin + ss * ss * sHMax) * effectiveConfiningStress / sVert;
  const E bXy = -cs * ss * (sHMax - sHMin) * effectiveConfiningStress / sVert;
  const E bXz = s11.of(0);
  const E bYy = -(ss * ss * sHMin + cs * cs * sHMax) * effectiveConfiningStress / sVert;
  const E bYz = s11.of(0);
  const E bZz = -effectiveConfiningStress;
  return Stress{bXx, bYy, bZz, bXy, bYz, bXz};
}

const std::vector<std::string>& stressNames() {
  static const std::vector<std::string> Names = {"b_xx", "b_yy", "b_zz", "b_xy", "b_yz", "b_xz"};
  return Names;
}

} // namespace

const EasiSpecial* easiSpecialOf(const std::string& tag) {
  static const std::vector<EasiSpecial> Specials = {
      {"!STRESS_STR_DIP_SLIP_AM",
       {"mu_d", "mu_s", "strike", "dip", "s_zz", "cohesion", "R", "DipSlipFaulting", "s2ratio"},
       stressNames()},
      {"!OptimalStress",
       {"mu_d",
        "mu_s",
        "strike",
        "dip",
        "rake",
        "effectiveConfiningStress",
        "cohesion",
        "R",
        "s2ratio"},
       stressNames()},
      {"!AndersonianStress",
       {"mu_d", "mu_s", "SH_max", "sig_zz", "cohesion", "S", "s2ratio", "S_v"},
       stressNames()}};
  for (const auto& special : Specials) {
    if (special.tag == tag) {
      return &special;
    }
  }
  return nullptr;
}

std::map<std::string, NodeId> easiSpecialValues(expr::Arena& arena,
                                                const EasiSpecial& special,
                                                const std::map<std::string, NodeId>& in) {
  const auto at = [&](const char* name) { return E(arena, in.at(name)); };
  const Stress stress = [&]() {
    if (special.tag == "!STRESS_STR_DIP_SLIP_AM") {
      return stressStrDipSlipAm(StressStrDipSlipAmIn{at("mu_d"),
                                                     at("mu_s"),
                                                     at("strike"),
                                                     at("dip"),
                                                     at("s_zz"),
                                                     at("cohesion"),
                                                     at("R"),
                                                     at("DipSlipFaulting"),
                                                     at("s2ratio")});
    }
    if (special.tag == "!OptimalStress") {
      return optimalStress(OptimalStressIn{at("mu_d"),
                                           at("mu_s"),
                                           at("strike"),
                                           at("dip"),
                                           at("rake"),
                                           at("effectiveConfiningStress"),
                                           at("cohesion"),
                                           at("R"),
                                           at("s2ratio")});
    }
    return andersonianStress(AndersonianStressIn{at("mu_d"),
                                                 at("mu_s"),
                                                 at("SH_max"),
                                                 at("sig_zz"),
                                                 at("cohesion"),
                                                 at("S"),
                                                 at("s2ratio"),
                                                 at("S_v")});
  }();
  return {{"b_xx", stress.bXx.node()},
          {"b_yy", stress.bYy.node()},
          {"b_zz", stress.bZz.node()},
          {"b_xy", stress.bXy.node()},
          {"b_yz", stress.bYz.node()},
          {"b_xz", stress.bXz.node()}};
}

} // namespace seissol::reader::scripting
