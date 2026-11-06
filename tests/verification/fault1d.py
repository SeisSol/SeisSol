#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The exact solution of a rupture on the fault1d meshes.

fault1d is a column that is periodic along the fault plane x = 0 and open
across it. A fault that starts out uniform stays uniform there: nothing varies
along it, so nothing is transferred along it, and each side radiates plane
waves straight away from it. The traction then depends on the slip rate V of
the same instant only,

    tau(t) = tau_load(t) - eta V(t),    eta = 1 / (1 / Z_left + 1 / Z_right),

where Z = rho c_s is the shear impedance of a side, tau_load the traction the
fault would carry if it were welded (the initial stress and the nucleations),
and eta V the radiation damping of the two half-spaces. With the friction law,
every point of the fault solves the same ordinary differential equation, and the
wave it sends out follows from it: at the distance |x| on a side of impedance Z
and shear wave speed c_s, the particle velocity is eta / Z V(t - |x| / c_s), and
the shear stress eta V(t - |x| / c_s).

In an anisotropic medium, the impedances are matrices, the radiation damping
depends on the direction of the slip, which then leaves the direction of
tau_load, and the slip changes the normal stress; each side sends out three
waves, polarized and as fast as the Christoffel matrix of the direction x
makes them (see Medium and Problem). The receivers are compared with all of
this in the fault frame that the documentation defines.

The equation is integrated here with steps far below the time step of the case,
so that what remains of the difference to SeisSol is the error of SeisSol. The
friction laws follow their documentation (docs/dynamic-rupture.rst and
docs/slip-rate-on-dr.rst) and the publications it points to, not the code they
are meant to check:

* the nucleation ramps with the smooth step G(t) of the SCEC benchmarks TPV101
  to TPV104;
* a forced rupture time weakens the friction linearly over t_0, as the time
  weakening of the SCEC benchmarks with forced rupture does;
* the regularized Yoffe function is the Yoffe function convolved with a
  triangle (Tinti et al., 2005), here by quadrature rather than in closed form;
* the Gaussian source time function is the rate of G(t), stretched to the rise
  time;
* thermal pressurization is the convolution of the frictional heating with the
  Green's function of the diffusion across a Gaussian shear zone, rather than a
  sum over wave numbers.

The delta pulse (FL=35), whose output is a single impulse, is not covered.
"""

import math
import re
from pathlib import Path

import numpy as np

# the slip rate above which a point counts as ruptured (docs/fault-output.rst)
RUPTURE_SLIP_RATE = 1e-3


class Unsupported(ValueError):
    """The inputs pose a problem that this module does not solve."""


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------

_TAG = re.compile(r"!(\w+)")
_ENTRY = re.compile(r"^(\s*)(?:-\s+)?([A-Za-z_]\w*)\s*:\s*(.*)$")


def read_constant_maps(path):
    """Read the constant maps of an easi file, per group.

    Only what the cases use is understood: a !ConstantMap, or an !Any of
    !GroupFilter components that each hold one. The result maps a group number,
    or None for every group, to the values of the map. Anything else is
    refused, rather than read as something it is not.
    """
    maps = {}
    groups = None
    map_indent = None
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines()):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        for tag in _TAG.findall(line):
            if tag not in ("ConstantMap", "Any", "GroupFilter"):
                raise Unsupported(f"{path}:{number + 1}: !{tag} is not supported")
            if tag == "GroupFilter":
                groups = None
        if map_indent is not None and indent > map_indent:
            entry = _ENTRY.match(line)
            if entry is None:
                raise Unsupported(f"{path}:{number + 1}: expected name: value")
            for group in groups or [None]:
                maps.setdefault(group, {})[entry.group(2)] = float(entry.group(3))
            continue
        map_indent = None
        entry = _ENTRY.match(line)
        if entry is None:
            continue
        key, value = entry.group(2), entry.group(3).strip()
        if key == "map" and not value:
            map_indent = indent
        elif key == "groups":
            groups = [int(item) for item in re.findall(r"-?\d+", value)]
    return maps


# the stress components in Voigt notation, in the order of the stiffness c_ij
VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")
# the Voigt indices of the components (x, x), (y, x) and (z, x): the traction on
# a plane x = const, and the strains of a plane wave that travels along x
_ALONG_X = (0, 5, 4)


class Medium:
    """A homogeneous elastic medium, by its density and its stiffness in Voigt
    notation, ordered as VOIGT.

    Along x, it carries three plane waves: their polarizations q are the
    eigenvectors of the Christoffel matrix Gamma_ik = C_i1k1, and their speeds
    c follow from its eigenvalues rho c^2. A wave leaving a plane x = const
    with the velocity q (q . v) carries the traction -rho c q (q . v) on it, so
    that the impedance Z = sum rho c q q^T = (rho Gamma)^(1/2) maps the
    velocity of the waves leaving a plane to the traction they carry, and the
    admittance Y = Z^-1 maps back.
    """

    def __init__(self, rho, stiffness):
        self.rho = float(rho)
        self.stiffness = np.array(stiffness, dtype=float)
        christoffel = self.stiffness[np.ix_(_ALONG_X, _ALONG_X)]
        values, vectors = np.linalg.eigh(christoffel)
        if values[0] <= 0.0:
            raise Unsupported("the medium carries no wave along x")
        self.speeds = np.sqrt(values / self.rho)
        self.polarizations = vectors
        self.impedance = vectors @ np.diag(self.rho * self.speeds) @ vectors.T
        self.admittance = vectors @ np.diag(1.0 / (self.rho * self.speeds)) @ vectors.T

    @classmethod
    def isotropic(cls, rho, mu, lam):
        stiffness = np.zeros((6, 6))
        stiffness[:3, :3] = lam
        for index in range(3):
            stiffness[index, index] = lam + 2.0 * mu
            stiffness[index + 3, index + 3] = mu
        return cls(rho, stiffness)


def medium_of(maps, group):
    """The medium of a group: rho, and either mu and lambda or the 21 c_ij of an
    anisotropic medium."""
    values = dict(maps.get(None, {}))
    values.update(maps.get(group, {}))
    try:
        if "c11" in values:
            stiffness = np.zeros((6, 6))
            for i in range(6):
                for j in range(i, 6):
                    stiffness[i, j] = stiffness[j, i] = values[f"c{i + 1}{j + 1}"]
            return Medium(values["rho"], stiffness)
        return Medium.isotropic(values["rho"], values["mu"], values["lambda"])
    except KeyError as error:
        raise Unsupported(f"the material of group {group} lacks {error}") from error


class Inputs:
    """The parameters of a case: the fault map first, then the parameter file.

    A value the fault map gives holds for the whole fault, which is all a
    uniform fault can use. Names of the parameter file are case-insensitive,
    those of easi are not.
    """

    def __init__(self, dynamic_rupture, fault):
        self.parameters = {key.lower(): value for key, value in dynamic_rupture.items()}
        self.fault = dict(fault)

    def has(self, *names):
        return any(
            name in self.fault or name.lower() in self.parameters for name in names
        )

    def get(self, *names, default=None):
        for name in names:
            if name in self.fault:
                return float(self.fault[name])
        for name in names:
            if name.lower() in self.parameters:
                return float(self.parameters[name.lower()])
        if default is None:
            raise Unsupported(f"the case gives none of {', '.join(names)}")
        return float(default)


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------


def smooth_step(t, rise):
    """The ramp G of the SCEC benchmarks TPV101 to TPV104, from 0 at t <= 0 to 1."""
    if rise <= 0.0:
        return 1.0 if t >= 0.0 else 0.0
    if t <= 0.0:
        return 0.0
    if t >= rise:
        return 1.0
    return math.exp((t - rise) ** 2 / (t * (t - 2.0 * rise)))


def smooth_step_rate(t, rise):
    """The time derivative of smooth_step."""
    if rise <= 0.0 or t <= 0.0 or t >= rise:
        return 0.0
    value = smooth_step(t, rise)
    if value == 0.0:
        return 0.0
    product = t * (t - 2.0 * rise)
    return value * (-2.0 * rise * rise * (t - rise) / (product * product))


class Loading:
    """The traction on the welded fault: the initial stress and the nucleations.

    The fault is the plane x = 0, so the traction of a stress tensor is its
    first column: the normal stress s_xx, negative in compression, and the
    shear traction (s_xy, s_xz) along y and z.
    """

    def __init__(self, inputs):
        if inputs.has("T_n", "T_s", "T_d", "Tnuc_n", "Tnuc_s", "Tnuc_d"):
            raise Unsupported("tractions in fault coordinates are not supported")
        self.initial = self._traction(inputs, "s_")
        self.sources = []
        count = int(inputs.get("nucleationcount", default=1))
        for index in range(count):
            suffix = "" if index == 0 else str(index + 1)
            prefix = f"nuc{suffix}_"
            if not any(f"{prefix}{c}" in inputs.fault for c in ("xx", "xy", "xz")):
                continue
            rise = inputs.get(
                f"nuc{suffix}_rise_time", f"t{suffix or ''}_0", default=0.0
            )
            onset = inputs.get(f"nuc{suffix}_onset", f"s{suffix or ''}_0", default=0.0)
            self.sources.append((self._traction(inputs, prefix), onset, rise))

    @staticmethod
    def _traction(inputs, prefix):
        return tuple(
            inputs.get(f"{prefix}{component}", default=0.0)
            for component in ("xx", "xy", "xz")
        )

    def __call__(self, t):
        normal, along_y, along_z = self.initial
        for (dn, dy, dz), onset, rise in self.sources:
            fraction = smooth_step(t - onset, rise)
            normal += fraction * dn
            along_y += fraction * dy
            along_z += fraction * dz
        return normal, along_y, along_z

    def direction_is_fixed(self):
        """Whether the shear traction keeps its direction, so that the slip does."""
        _, y0, z0 = self.initial
        for (_, dy, dz), _, _ in self.sources:
            if abs(y0 * dz - z0 * dy) > 1e-12 * math.hypot(y0, z0) * math.hypot(dy, dz):
                return False
        return True


# --------------------------------------------------------------------------
# friction laws
# --------------------------------------------------------------------------


def _log_sinh(x):
    """log(sinh(x)) for x > 0, without overflow."""
    return x + math.log(-math.expm1(-2.0 * x) / 2.0)


def _arsinh_exp(log_z):
    """asinh(exp(log_z)), also where exp(log_z) overflows."""
    if log_z > 20.0:
        return log_z + math.log(2.0) + math.exp(-2.0 * log_z) / 4.0
    return math.asinh(math.exp(log_z))


def _arsinh_exp_slope(log_z):
    """The derivative of _arsinh_exp."""
    return 1.0 / math.sqrt(1.0 + math.exp(-2.0 * min(log_z, 350.0)))


class SlipWeakening:
    """Linear slip weakening (FL=16, FL=6) and its TP proxy (FL=1058).

    The strength is -C - min(0, sigma_n) mu, with mu a function of the
    accumulated slip S, and of the time where a forced rupture time is given.
    The friction coefficient does not depend on the slip rate, so the slip rate
    is explicit: with the normal stress sigma_n - k V of a slip rate V,
    V = max(0, tau - strength(sigma_n)) / (eta + mu k), for the traction tau
    and the radiation damping eta along the slip. With FL=6, the strength
    itself relaxes towards that of the slip-weakening law at the rate
    (|V| + v*) / L, the form of the Prakash-Clifton regularization that Cochard
    and Rice (2000) use for faults between different media.
    """

    def __init__(self, law, inputs):
        self.law = law
        self.mu_s = inputs.get("mu_s")
        self.mu_d = inputs.get("mu_d")
        self.d_c = inputs.get("d_c")
        self.cohesion = inputs.get("cohesion", default=0.0)
        if inputs.get("lsw_healingThreshold", default=-1.0) >= 0.0:
            raise Unsupported("healing is not supported")
        self.forced = None
        if inputs.has("forced_rupture_time"):
            if law == 1058:
                raise Unsupported("a forced rupture time with FL=1058 is not supported")
            self.forced = (
                inputs.get("forced_rupture_time"),
                inputs.get("t_0", default=0.0),
            )
        self.exponent = inputs.get("TpProxyExponent", default=1.0 / 3.0)
        if law == 6:
            self.v_star = inputs.get("pc_vstar", "v_star")
            self.length = inputs.get("pc_prakashLength", "L")

    def initial_state(self, loading):
        # slip along y and z, accumulated slip, and with FL=6 the strength
        state = [0.0, 0.0, 0.0]
        if self.law == 6:
            state.append(self.strength(0.0, state, loading(0.0)[0]))
        return state

    def strength(self, t, state, normal):
        return -self.cohesion - min(0.0, normal) * self.friction(t, state)

    def friction(self, t, state):
        slip = state[2]
        if self.law == 1058:
            return self.mu_d + (self.mu_s - self.mu_d) / (1.0 + slip / self.d_c) ** (
                self.exponent
            )
        weakening = min(slip / self.d_c, 1.0)
        if self.forced is not None:
            time, rise = self.forced
            if rise > 0.0:
                weakening = max(weakening, min(max((t - time) / rise, 0.0), 1.0))
            elif t >= time:
                weakening = 1.0
        return self.mu_s - (self.mu_s - self.mu_d) * weakening

    def evaluate(self, t, state, normal, shear, eta, pressure, coupling=0.0):
        friction = self.friction(t, state)
        if self.law == 6:
            rate = max(0.0, shear - state[3]) / eta
        else:
            # the strength grows by mu k V as the slip rate clamps the fault
            divisor = eta + (friction * coupling if normal < 0.0 else 0.0)
            if divisor <= 0.0:
                raise Unsupported("the fault loses its strength faster than it slips")
            rate = max(0.0, shear - self.strength(t, state, normal)) / divisor
        normal_stress = normal - coupling * rate
        if (normal_stress < 0.0) != (normal < 0.0):
            raise Unsupported("the slip rate opens or closes the fault")
        result = {
            "slip-rate": rate,
            "friction": friction,
            "traction": shear - eta * rate,
            "normal-stress": normal_stress,
        }
        extra = []
        if self.law == 6:
            strength = self.strength(t, state, normal_stress)
            extra.append(-(rate + self.v_star) / self.length * (state[3] - strength))
        return result, extra

    def weakened(self, series):
        """When the fault weakened fully: the accumulated slip reaches d_c."""
        return _first_crossing(series["time"], series["accumulated-slip"], self.d_c)


class RateAndState:
    """The rate-and-state laws FL=3 (aging), FL=4 (slip), FL=7 (severe
    velocity weakening) and FL=103 (fast velocity weakening), optionally with
    thermal pressurization.

    The friction coefficient of FL=3, 4 and 103 is
    a asinh(V / (2 V_0) exp(Psi / a)), with Psi = f_0 + b ln(V_0 psi / L) for
    aging and slip law, and Psi = psi for FL=103; that of FL=7 is
    f_0 + a V / (V + V_0) - b psi / (psi + L). The slip rate solves
    tau - eta V = (sigma_eff + k V) f(V, psi), for the traction tau and the
    radiation damping eta along the slip and the effective normal stress
    sigma_eff + k V that the slip rate leaves, positive in compression here. It
    has at most one root where f and the normal stress grow with V, and, for
    k < 0, where eta outweighs |k| f. The asinh always has one; FL=7 has none
    where the traction does not exceed sigma_eff f(0, psi), and the fault is
    locked there.
    """

    def __init__(self, law, inputs, loading):
        self.law = law
        self.a = inputs.get("rs_a")
        self.b = inputs.get("rs_b")
        self.length = inputs.get("rs_sl0")
        self.v0 = inputs.get("rs_sr0")
        self.f0 = inputs.get("rs_f0")
        if law == 103:
            self.mu_w = inputs.get("rs_muW", "mu_w")
            self.v_w = inputs.get("rs_srW")
        initial_rate = math.hypot(
            inputs.get("rs_iniSlipRate1"), inputs.get("rs_iniSlipRate2")
        )
        # the state that holds the traction at the start, the initial stress
        # with every nucleation in effect then, at the initial slip rate
        normal, along_y, along_z = loading(0.0)
        friction = math.hypot(along_y, along_z) / -normal
        if law == 7:
            direct = self.a * initial_rate / (initial_rate + self.v0)
            weakening = (self.f0 + direct - friction) / self.b
            if weakening < 1.0:
                self.psi0 = self.length * weakening / (1.0 - weakening)
            else:
                # no state slips at the initial slip rate: locked, and healed
                self.psi0 = self.length * initial_rate / self.v0
            return
        psi = self.a * (
            math.log(2.0)
            + _log_sinh(friction / self.a)
            - math.log(initial_rate / self.v0)
        )
        if law == 103:
            self.psi0 = psi
        else:
            self.psi0 = self.length / self.v0 * math.exp((psi - self.f0) / self.b)

    def initial_state(self, loading):
        return [0.0, 0.0, 0.0, self.psi0]

    def big_psi(self, psi):
        if self.law == 103:
            return psi
        return self.f0 + self.b * math.log(self.v0 * psi / self.length)

    def severe_friction(self, rate, psi):
        """The friction coefficient of FL=7."""
        return (
            self.f0
            + self.a * rate / (rate + self.v0)
            - self.b * psi / (psi + self.length)
        )

    def slip_rate(self, shear, sigma, psi, eta, coupling=0.0):
        """Solve shear - eta V = (sigma + coupling V) f(V, psi) for V, by Newton
        in log V; 0 where no slip rate solves it."""
        if shear <= 0.0:
            return 0.0
        if sigma + coupling * shear / eta <= 0.0:
            raise Unsupported("the slip rate may open the fault")
        if self.law == 7:

            def friction(u):
                rate = math.exp(u)
                slope = self.a * rate * self.v0 / (rate + self.v0) ** 2
                return self.severe_friction(rate, psi), slope

        else:
            offset = self.big_psi(psi) / self.a - math.log(2.0 * self.v0)

            def friction(u):
                log_z = offset + u
                return self.a * _arsinh_exp(log_z), self.a * _arsinh_exp_slope(log_z)

        def residual(u):
            value, slope = friction(u)
            rate = math.exp(u)
            normal = sigma + coupling * rate
            return (
                normal * value + eta * rate - shear,
                normal * slope + (coupling * value + eta) * rate,
            )

        upper = math.log(shear / eta)
        lower = upper - 1.0
        while residual(lower)[0] > 0.0:
            lower -= 2.0 * (upper - lower)
            if lower < -800.0:
                return 0.0
        u = 0.5 * (lower + upper)
        for _ in range(200):
            value, slope = residual(u)
            if value > 0.0:
                upper = u
            else:
                lower = u
            step = value / slope
            candidate = u - step
            if not lower < candidate < upper:
                candidate = 0.5 * (lower + upper)
            if abs(candidate - u) < 1e-15 * max(1.0, abs(u)):
                return math.exp(candidate)
            u = candidate
        return math.exp(u)

    def state_rate(self, rate, psi):
        if self.law == 3:
            return 1.0 - rate * psi / self.length
        if self.law == 7:
            return rate - self.v0 / self.length * psi
        if self.law == 4:
            product = rate * psi / self.length
            return -product * math.log(product) if product > 0.0 else 0.0
        # FL=103
        ratio = rate / self.v0
        f_lv = self.f0 - (self.b - self.a) * math.log(ratio)
        mu_ss = self.mu_w + (f_lv - self.mu_w) / (1.0 + (rate / self.v_w) ** 8) ** 0.125
        psi_ss = self.a * (math.log(2.0) - math.log(ratio) + _log_sinh(mu_ss / self.a))
        return -rate / self.length * (psi - psi_ss)

    def evaluate(self, t, state, normal, shear, eta, pressure, coupling=0.0):
        sigma = -(normal + pressure)
        if sigma <= 0.0:
            raise Unsupported("the effective normal stress is not compressive")
        psi = state[3]
        rate = self.slip_rate(shear, sigma, psi, eta, coupling)
        traction = shear - eta * rate
        if self.law == 7:
            # where the fault is locked, the friction exceeds the traction
            friction = self.severe_friction(rate, psi)
        else:
            friction = traction / (sigma + coupling * rate)
        result = {
            "slip-rate": rate,
            "friction": friction,
            "traction": traction,
            "normal-stress": normal - coupling * rate,
            "state": psi,
        }
        return result, [self.state_rate(rate, psi)]

    def weakened(self, series):
        """When the friction coefficient falls to mu_w + 0.05 (f_0 - mu_w),
        not before the rupture time; only FL=103 defines it."""
        if self.law != 103:
            return None
        start = series["rupture-time"]
        if start is None:
            return None
        threshold = self.mu_w + 0.05 * (self.f0 - self.mu_w)
        time = series["time"]
        after = time >= start
        crossing = _first_crossing(time[after], -series["friction"][after], -threshold)
        return crossing


class ThermalPressurization:
    """Temperature and pore pressure at the centre of a Gaussian shear zone.

    The heat tau V is released over a shear zone of the profile
    exp(-y^2 / (2 w^2)) / (sqrt(2 pi) w) and diffuses across it, as the pore
    pressure does, which the heating drives with Lambda dT/dt. At the centre,

        T = T_0 + 1 / (rho c) int tau V(t') G_th(t - t') dt',
        p = p_0 + Lambda / (rho c) int tau V(t')
                  (a_hy G_hy - a_th G_th)(t - t') / (a_hy - a_th) dt',

    with G_a(s) = 1 / sqrt(2 pi (w^2 + 2 a s)). The heating is integrated as a
    product rule over the steps of the integration, with the kernels integrated
    exactly over each step; on steps of equal width, the kernels depend on the
    distance in steps only and are tabulated once.
    """

    def __init__(self, inputs):
        self.a_th = inputs.get("tp_thermalDiffusivity", "alpha_th")
        self.rho_c = inputs.get("tp_heatCapacity", "rho_c")
        self.lam = inputs.get("tp_undrainedTPResponse", "tp_lambda")
        self.t0 = inputs.get("tp_iniTemp", "iniTemp")
        if inputs.get("tp_iniPressure", "iniPressure") != 0.0:
            # how it enters the normal stress of the fault file is a convention
            # of the input, which this module does not reproduce
            raise Unsupported("an initial pore pressure is not supported")
        self.a_hy = inputs.get("tp_hydraulicDiffusivity")
        self.width = inputs.get("tp_halfWidthShearZone")
        if abs(self.a_hy - self.a_th) <= 1e-9 * max(self.a_hy, self.a_th):
            raise Unsupported("equal hydraulic and thermal diffusivity")

    def _root(self, diffusivity, lag):
        return np.sqrt(self.width * self.width + 2.0 * diffusivity * lag)

    def _kernels(self, low, high):
        """The integrals of G_th and of the pressure kernel over the lags [low, high]."""
        root_th_high, root_th_low = self._root(self.a_th, high), self._root(
            self.a_th, low
        )
        root_hy_high, root_hy_low = self._root(self.a_hy, high), self._root(
            self.a_hy, low
        )
        scale = 1.0 / math.sqrt(2.0 * math.pi)
        thermal = 2.0 * scale * (high - low) / (root_th_high + root_th_low)
        pressure = (
            scale
            * ((root_hy_high - root_hy_low) - (root_th_high - root_th_low))
            / (self.a_hy - self.a_th)
        )
        return thermal, pressure

    def prepare(self, step, count):
        """Tabulate the kernels for `count` steps of width `step`.

        A step j < n, seen from t_n + c step, lies at the lags
        [(n - j - 1 + c) step, (n - j + c) step]; the step n itself contributes
        from the lag 0 to c step.
        """
        distance = np.arange(1, count + 2, dtype=np.float64)
        self.history = {
            fraction: self._kernels(
                (distance - 1.0 + fraction) * step, (distance + fraction) * step
            )
            for fraction in (0.0, 0.5, 1.0)
        }
        self.local = {
            fraction: self._kernels(np.float64(0.0), np.float64(fraction * step))
            for fraction in (0.5, 1.0)
        }
        self.heating = np.zeros(count + 1)
        self.recorded = 0

    def changes(self, fraction=0.0, current=0.0):
        """Temperature and pressure change at the fraction of the step that
        follows the recorded ones, with the heating `current` within it."""
        count = self.recorded
        past = self.heating[:count][::-1]
        thermal, pressure = self.history[fraction]
        temperature = float(np.dot(past, thermal[:count]))
        change = float(np.dot(past, pressure[:count]))
        if fraction > 0.0:
            thermal, pressure = self.local[fraction]
            temperature += current * float(thermal)
            change += current * float(pressure)
        return temperature / self.rho_c, self.lam * change / self.rho_c

    def record(self, heating):
        """Close a step with its average heating."""
        self.heating[self.recorded] = heating
        self.recorded += 1


# --------------------------------------------------------------------------
# imposed slip rates
# --------------------------------------------------------------------------

_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(24)


def _yoffe_convolution(t, tau_s, tau_r, cumulative):
    """The Yoffe function of duration tau_r, or its integral, convolved with
    the triangle of half-width tau_s, both of unit area.

    With s = tau_r sin^2(phi), the Yoffe function becomes (4 / pi) cos^2(phi)
    dphi, which removes its singularities at both ends; the triangle (or its
    integral) is piecewise polynomial in s, so the integral is split at its
    kinks and taken by Gauss-Legendre on each piece, which is exact to rounding.
    """

    def triangle(u):
        if u <= 0.0:
            return 0.0
        if u <= tau_s:
            return u * u / (2.0 * tau_s * tau_s) if cumulative else u / (tau_s * tau_s)
        if u <= 2.0 * tau_s:
            rest = 2.0 * tau_s - u
            return (
                1.0 - rest * rest / (2.0 * tau_s * tau_s)
                if cumulative
                else rest / (tau_s * tau_s)
            )
        return 1.0 if cumulative else 0.0

    def angle(s):
        return math.asin(math.sqrt(min(max(s / tau_r, 0.0), 1.0)))

    breaks = sorted({0.0, math.pi / 2.0} | {angle(t - k * tau_s) for k in (0, 1, 2)})
    total = 0.0
    for low, high in zip(breaks[:-1], breaks[1:]):
        if high <= low:
            continue
        half, mid = 0.5 * (high - low), 0.5 * (high + low)
        for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS):
            phi = mid + half * node
            s = tau_r * math.sin(phi) ** 2
            total += weight * half * math.cos(phi) ** 2 * triangle(t - s)
    return 4.0 / math.pi * total


class ImposedSlipRate:
    """The imposed slip rates, FL=33 (regularized Yoffe) and FL=34 (Gaussian).

    The fault slips by (strike_slip, dip_slip) after rupture_onset; the traction
    is whatever the slip rate leaves of the initial one, which these cases do
    not give, so it is -eta V.
    """

    def __init__(self, law, inputs, loading):
        if any(loading.initial) or loading.sources:
            raise Unsupported("imposed slip rates with an initial stress")
        self.law = law
        self.slip = math.hypot(inputs.get("strike_slip"), inputs.get("dip_slip"))
        self.direction = (inputs.get("strike_slip"), inputs.get("dip_slip"))
        self.onset = inputs.get("rupture_onset")
        if law == 33:
            self.tau_s = inputs.get("tau_S")
            self.tau_r = max(inputs.get("tau_R"), self.tau_s)
        else:
            self.rise = inputs.get("rupture_rise_time")

    def __call__(self, t):
        local = t - self.onset
        if self.law == 33:
            rate = _yoffe_convolution(local, self.tau_s, self.tau_r, cumulative=False)
            slip = _yoffe_convolution(local, self.tau_s, self.tau_r, cumulative=True)
        else:
            rate = smooth_step_rate(local, self.rise)
            slip = smooth_step(local, self.rise)
        return self.slip * rate, self.slip * slip


# --------------------------------------------------------------------------
# the problem and its solution
# --------------------------------------------------------------------------


def _first_crossing(time, values, threshold):
    """The first time `values` reaches `threshold` from below, interpolated."""
    above = np.nonzero(values >= threshold)[0]
    if above.size == 0:
        return None
    index = int(above[0])
    if index == 0:
        return float(time[0])
    before, after = values[index - 1], values[index]
    fraction = (threshold - before) / (after - before) if after != before else 1.0
    return float(time[index - 1] + fraction * (time[index] - time[index - 1]))


class Problem:
    """The 0D problem of a case on fault1d.

    `dynamic_rupture` maps the names of the DynamicRupture section of the
    parameter file to their values, `fault` is the constant map of the fault,
    and `left` and `right` are the media at x < 0 and at x > 0.

    A slip rate V = v(x > 0) - v(x < 0), which lies in the fault plane, leaves
    the traction t = -H V on the fault, H = (Y_left + Y_right)^-1, since the
    waves that leave the two sides carry t = -Z_right v(x > 0) and
    t = Z_left v(x < 0). Its tangential block H_tt is the radiation damping,
    which an anisotropic medium makes depend on the direction of the slip, and
    its normal row H_nt changes the normal stress by -H_nt V. Friction acts
    against the slip, so that the traction the slip leaves is the strength S
    along the direction n of the slip: tau_load - V H_tt n = S n. For an
    isotropic medium, H_tt = eta I, and n is the direction of tau_load.
    """

    def __init__(self, dynamic_rupture, fault, left, right):
        self.inputs = Inputs(dynamic_rupture, fault)
        self.law = int(self.inputs.get("FL"))
        self.left, self.right = left, right
        self.impedance = np.linalg.inv(left.admittance + right.admittance)
        self.impedance = 0.5 * (self.impedance + self.impedance.T)
        self.eta_shear = self.impedance[1:, 1:]
        self.eta_normal = self.impedance[0, 1:]
        scale = float(np.max(np.abs(self.impedance)))
        self.isotropic = (
            abs(self.eta_shear[0, 0] - self.eta_shear[1, 1])
            + abs(self.eta_shear[0, 1])
            + float(np.max(np.abs(self.eta_normal)))
            <= 1e-14 * scale
        )
        self.eta = float(self.eta_shear[0, 0])
        self._direction = None
        self.loading = Loading(self.inputs)
        self.thermal = None
        if self.law in (2, 16, 6, 1058):
            self.friction = SlipWeakening(
                16 if self.law == 2 else self.law, self.inputs
            )
        elif self.law in (3, 4, 7, 103):
            self.friction = RateAndState(self.law, self.inputs, self.loading)
            if int(self.inputs.get("thermalPress", default=0)):
                self.thermal = ThermalPressurization(self.inputs)
        elif self.law in (33, 34):
            self.friction = ImposedSlipRate(self.law, self.inputs, self.loading)
        else:
            raise Unsupported(f"friction law {self.law} has no exact solution here")

    @property
    def frame(self):
        """The map from the y and z components of a vector in the fault plane
        to the strike and dip components the fault receivers report.

        As docs/dynamic-rupture.rst and docs/parameter-file.rst state it: the
        normal n points from the side of the reference point (or, with
        refPointMethod=1, against the reference direction) to the other side;
        the strike is (n_y, -n_x, 0) / |.|, the dip
        (-s_y n_z, s_x n_z, s_y n_x - n_y s_x) / |.|; the traction is sigma n,
        and the slip rate the velocity on the side n points to minus that on
        the other, which is the side of x > 0 for n = e_x.
        """
        reference = self.inputs.get("XRef", default=0.0)
        if reference == 0.0:
            raise Unsupported("the reference point does not tell the sides apart")
        method = int(self.inputs.get("refPointMethod", default=0))
        side = 1.0 if reference > 0.0 else -1.0
        n_x = -side if method == 0 else side
        normal = np.array((n_x, 0.0, 0.0))
        strike = np.array((normal[1], -normal[0], 0.0)) / math.hypot(*normal[:2])
        dip = np.array(
            (
                -strike[1] * normal[2],
                strike[0] * normal[2],
                strike[1] * normal[0] - normal[1] * strike[0],
            )
        )
        dip /= np.linalg.norm(dip)
        return n_x * np.array((strike[1:], dip[1:]))

    @property
    def imposed(self):
        return isinstance(self.friction, ImposedSlipRate)

    # the solution ------------------------------------------------------

    def solve(self, end, step):
        """Integrate up to `end` with steps of `step`; returns the series."""
        count = int(round(end / step))
        self._direction = None
        if self.imposed:
            return self._imposed(count, step)
        return self._integrate(count, step)

    def _imposed(self, count, step):
        time = np.arange(count + 1) * step
        rate, slip = np.array([self.friction(t) for t in time]).T
        direction = self.frame.T @ np.asarray(self.friction.direction)
        direction /= np.linalg.norm(direction)
        load = np.array([self.loading(t) for t in time])
        vector = np.outer(rate, direction)
        traction = load[:, 1:] - vector @ self.eta_shear.T
        normal = load[:, 0] - vector @ self.eta_normal
        return self._finish(
            {
                "time": time,
                "slip-rate": rate,
                "slip-rate-y": vector[:, 0],
                "slip-rate-z": vector[:, 1],
                "traction": np.hypot(traction[:, 0], traction[:, 1]),
                # the traction the fault would carry if it were welded
                "locked-traction": np.hypot(load[:, 1], load[:, 2]),
                "normal-stress": normal,
                "effective-normal-stress": normal,
                "slip": slip,
                "accumulated-slip": slip,
            }
        )

    def _rates(self, t, state, pressure):
        normal, along_y, along_z = self.loading(t)
        load = np.array((along_y, along_z))
        shear = math.hypot(along_y, along_z)
        direction = load / shear if shear > 0.0 else np.array((1.0, 0.0))
        if not self.isotropic and self._direction is not None:
            direction = self._direction
        for _ in range(200):
            eta = float(direction @ self.eta_shear @ direction)
            coupling = float(self.eta_normal @ direction)
            result, extra = self.friction.evaluate(
                t, state, normal, float(direction @ load), eta, pressure, coupling
            )
            rate = result["slip-rate"]
            if self.isotropic or rate == 0.0:
                if rate == 0.0 and shear > 0.0:
                    direction = load / shear
                break
            # the direction in which the traction the slip leaves is the strength
            strength = result["traction"]
            update = np.linalg.solve(strength * np.eye(2) + rate * self.eta_shear, load)
            update /= np.linalg.norm(update)
            change = float(np.max(np.abs(update - direction)))
            direction = update
            if change <= 1e-15:
                break
        else:
            raise Unsupported("the direction of the slip does not settle")
        self._direction = direction
        rates = [rate * direction[0], rate * direction[1], rate] + extra
        return rates, result, direction

    def _integrate(self, count, step):
        state = self.friction.initial_state(self.loading)
        thermal = self.thermal
        if thermal is not None:
            thermal.prepare(step, count)
        rows = []

        def observe(t, state):
            pressure = temperature = 0.0
            if thermal is not None:
                temperature, pressure = thermal.changes()
            _, result, direction = self._rates(t, state, pressure)
            row = dict(result)
            row["time"] = t
            row["effective-normal-stress"] = result["normal-stress"] + pressure
            row["slip-rate-y"] = result["slip-rate"] * direction[0]
            row["slip-rate-z"] = result["slip-rate"] * direction[1]
            row["slip"] = math.hypot(state[0], state[1])
            row["accumulated-slip"] = state[2]
            if thermal is not None:
                row["pressure-change"] = pressure
                row["temperature"] = thermal.t0 + temperature
            return row

        rows.append(observe(0.0, state))
        for index in range(count):
            t = index * step
            heating = []
            slopes = []
            current = rows[-1]["traction"] * rows[-1]["slip-rate"]
            for fraction, previous in ((0.0, None), (0.5, 0), (0.5, 1), (1.0, 2)):
                if previous is None:
                    probe = state
                else:
                    scale = step * (1.0 if fraction == 1.0 else 0.5)
                    probe = [s + scale * k for s, k in zip(state, slopes[previous])]
                pressure = 0.0
                if thermal is not None:
                    _, pressure = thermal.changes(fraction, current)
                rates, result, _ = self._rates(t + fraction * step, probe, pressure)
                current = result["traction"] * result["slip-rate"]
                slopes.append(rates)
                heating.append(current)
            state = [
                s + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
                for s, k1, k2, k3, k4 in zip(state, *slopes)
            ]
            if thermal is not None:
                weights = (1.0, 2.0, 2.0, 1.0)
                average = sum(w * h for w, h in zip(weights, heating)) / 6.0
                thermal.record(average)
            rows.append(observe((index + 1) * step, state))
        series = {key: np.array([row[key] for row in rows]) for key in rows[0]}
        return self._finish(series)

    def _finish(self, series):
        time = series["time"]
        series["peak-slip-rate"] = np.maximum.accumulate(series["slip-rate"])
        series["rupture-time"] = _first_crossing(
            time, series["slip-rate"], RUPTURE_SLIP_RATE
        )
        if not self.imposed:
            series["dynamic-stress-time"] = self.friction.weakened(series)
        return series

    # what the receivers see --------------------------------------------

    def radiated(self, series, x):
        """The particle velocity (x, y, z) and the stress (ordered as VOIGT)
        of the waves at the distance x, as the receivers record them: without
        the initial stress, and zero before the waves arrive.

        On the fault, the slip rate V leaves the velocity Y H V on the side
        x > 0 and -Y H V on the side x < 0, Y being the admittance of the side.
        Each of the three waves of a side carries its part q (q . v) away at
        its own speed c, with the strain du/dx = -/+ q (q . v) / c on the side
        x > 0 / x < 0, and the stress follows from the strain.
        """
        medium = self.left if x < 0.0 else self.right
        side = -1.0 if x < 0.0 else 1.0
        time = series["time"]
        slip_rate = np.zeros((time.size, 3))
        slip_rate[:, 1] = series["slip-rate-y"]
        slip_rate[:, 2] = series["slip-rate-z"]
        start = side * slip_rate @ (medium.admittance @ self.impedance).T
        velocity = np.zeros((time.size, 3))
        gradient = np.zeros((time.size, 3))
        for speed, polarization in zip(medium.speeds, medium.polarizations.T):
            delayed = np.interp(
                time - abs(x) / speed, time, start @ polarization, left=0.0
            )
            velocity += np.outer(delayed, polarization)
            gradient += np.outer(-side * delayed / speed, polarization)
        stress = gradient @ medium.stiffness[:, _ALONG_X].T
        return velocity, stress


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------


def self_test():
    """Hold the integration against what is known in closed form.

    Returns the failures, as text.
    """
    failures = []
    medium = Medium.isotropic(2700.0, 3.24e10, 3.24e10)
    eta = medium.impedance[1, 1] / 2.0

    def check(label, deviation, tolerance):
        status = "ok" if deviation <= tolerance else "FAILED"
        print(
            f"[{status}] {label}: deviation {deviation:.3g}, tolerance {tolerance:.3g}"
        )
        if deviation > tolerance:
            failures.append(label)

    # linear slip weakening under a constant overload Delta: until d_c, the
    # slip grows as Delta / W (exp(W t / eta) - 1), with the weakening rate
    # W = (mu_s - mu_d) sigma / d_c; then at the rate (tau - tau_d) / eta
    sigma, tau, mu_s, mu_d, d_c = 120.0e6, 85.0e6, 0.677, 0.525, 0.08
    problem = Problem(
        {"FL": 16},
        {"s_xx": -sigma, "s_xy": tau, "mu_s": mu_s, "mu_d": mu_d, "d_c": d_c},
        medium,
        medium,
    )
    series = problem.solve(0.06, 1e-5)
    overload, weakening = tau - mu_s * sigma, (mu_s - mu_d) * sigma / d_c
    weakened = eta / weakening * math.log(1.0 + weakening * d_c / overload)
    exact = np.where(
        series["time"] < weakened,
        overload / weakening * np.expm1(weakening * series["time"] / eta),
        d_c + (tau - mu_d * sigma) / eta * (series["time"] - weakened),
    )
    check(
        "slip weakening under a constant overload",
        float(np.max(np.abs(series["slip"] - exact)) / np.max(exact)),
        1e-8,
    )

    # the imposed slip is the integral of the imposed slip rate, and all of it
    for law, extra in (
        (33, {"tau_S": 0.01, "tau_R": 0.05}),
        (34, {"rupture_rise_time": 0.05}),
    ):
        fault = {"strike_slip": 1.0, "dip_slip": 0.0, "rupture_onset": 0.01, **extra}
        problem = Problem({"FL": law, "XRef": -1.0e5}, fault, medium, medium)
        series = problem.solve(0.1, 1e-5)
        integral = np.concatenate(
            (
                [0.0],
                np.cumsum(
                    np.diff(series["time"])
                    * 0.5
                    * (series["slip-rate"][1:] + series["slip-rate"][:-1])
                ),
            )
        )
        check(
            f"slip and slip rate of FL={law} agree",
            float(np.max(np.abs(integral - series["slip"]))),
            1e-5,
        )
        check(f"FL={law} imposes all of the slip", abs(series["slip"][-1] - 1.0), 1e-12)

    # the kernels of thermal pressurization, tabulated per step, against a
    # direct quadrature of the convolution for a prescribed heating
    fault = {"tp_hydraulicDiffusivity": 4.0e-4, "tp_halfWidthShearZone": 0.005}
    parameters = {
        "tp_thermalDiffusivity": 1.0e-6,
        "tp_heatCapacity": 2.7e6,
        "tp_undrainedTPResponse": 0.1e6,
        "tp_iniTemp": 483.15,
        "tp_iniPressure": 0.0,
    }
    thermal = ThermalPressurization(Inputs(parameters, fault))
    step, count = 1e-4, 1000
    thermal.prepare(step, count)

    def heating(t):
        return 1.0e8 * math.sin(math.pi * t / 0.1) ** 2

    for index in range(count):
        # the average of the heating over the step, by Gauss-Legendre
        start = index * step
        nodes = start + 0.5 * step * (_GAUSS_NODES + 1.0)
        thermal.record(0.5 * float(np.dot(_GAUSS_WEIGHTS, [heating(t) for t in nodes])))
    temperature, pressure = thermal.changes()
    lags = np.linspace(0.0, count * step, 400001)
    weights = np.full(lags.size, lags[1] - lags[0])
    weights[[0, -1]] *= 0.5
    history = np.array([heating(count * step - lag) for lag in lags])

    def green(diffusivity):
        return 1.0 / np.sqrt(
            2.0 * math.pi * (thermal.width**2 + 2.0 * diffusivity * lags)
        )

    exact_temperature = (
        float(np.dot(weights, history * green(thermal.a_th))) / thermal.rho_c
    )
    kernel = (
        thermal.a_hy * green(thermal.a_hy) - thermal.a_th * green(thermal.a_th)
    ) / (thermal.a_hy - thermal.a_th)
    exact_pressure = (
        thermal.lam * float(np.dot(weights, history * kernel)) / thermal.rho_c
    )
    check(
        "temperature of thermal pressurization",
        abs(temperature - exact_temperature) / exact_temperature,
        1e-6,
    )
    check(
        "pore pressure of thermal pressurization",
        abs(pressure - exact_pressure) / exact_pressure,
        1e-6,
    )

    # the integration of a rate-and-state law has converged
    fault = {
        "s_xx": -120.0e6,
        "s_xy": 75.0e6,
        "nuc_xy": 40.0e6,
        "rs_a": 0.008,
        "rs_sl0": 0.02,
    }
    parameters = {
        "FL": 3,
        "t_0": 0.0625,
        "rs_f0": 0.6,
        "rs_sr0": 1e-6,
        "rs_b": 0.012,
        "rs_iniSlipRate1": 1e-16,
        "rs_iniSlipRate2": 0.0,
    }
    coarse = Problem(parameters, fault, medium, medium).solve(0.1, 2e-5)
    fine = Problem(parameters, fault, medium, medium).solve(0.1, 1e-5)
    check(
        "rate and state with half the step",
        float(
            np.max(np.abs(coarse["slip-rate"] - fine["slip-rate"][::2]))
            / np.max(fine["slip-rate"])
        ),
        1e-8,
    )

    # severe velocity weakening: the initial state holds the initial traction at
    # the initial slip rate; the healing then locks the fault, and while it is
    # locked, the state relaxes as psi_0 exp(-V_0 t / L)
    fault = {"s_xx": -120.0e6, "s_xy": 70.8e6, "rs_a": 0.008, "rs_sl0": 1.0e-7}
    parameters = {
        "FL": 7,
        "rs_f0": 0.6,
        "rs_sr0": 1e-6,
        "rs_b": 0.012,
        "rs_iniSlipRate1": 1e-16,
        "rs_iniSlipRate2": 0.0,
    }
    problem = Problem(parameters, fault, medium, medium)
    law = problem.friction
    check(
        "the initial state of FL=7",
        abs(law.severe_friction(1e-16, law.psi0) - 70.8e6 / 120.0e6),
        1e-14,
    )
    series = problem.solve(0.1, 1e-4)
    healed = law.psi0 * np.exp(-law.v0 / law.length * series["time"])
    check(
        "healing of FL=7 on the locked fault",
        float(np.max(np.abs(series["state"] - healed)) / law.psi0),
        1e-10,
    )
    check("FL=7 stays locked", float(np.max(series["slip-rate"][1:])), 0.0)

    # the tilted media of the anisotropic cases: the waves they radiate solve
    # rho dv/dt = dt/dx, d(sigma)/dt = C dv/dx with the traction t continuous
    # across the fault and the velocity jumping by the slip rate; what remains is
    # the error of the differences taken, and of the interpolation in time
    maps = read_constant_maps(
        Path(__file__).parent / "cases" / "material-fault1d-anisotropic.yaml"
    )
    left, right = medium_of(maps, 1), medium_of(maps, 2)
    # the source time function ends with a kink in the slope of the slip rate,
    # which the window leaves out, and starts smoothly
    fault = {"strike_slip": 1.0, "dip_slip": 0.5, "rupture_onset": 0.002}
    fault["rupture_rise_time"] = 0.04
    problem = Problem({"FL": 34, "XRef": -1.0e5}, fault, left, right)
    step, offset = 2e-7, 0.002
    series = problem.solve(0.03, step)
    slip_rate = np.zeros((series["time"].size, 3))
    slip_rate[:, 1] = series["slip-rate-y"]
    slip_rate[:, 2] = series["slip-rate-z"]
    below, below_stress = problem.radiated(series, -1e-9)
    above, above_stress = problem.radiated(series, 1e-9)
    check(
        "the waves of a tilted pair of media continue the slip rate",
        float(np.max(np.abs(above - below - slip_rate)) / np.max(np.abs(slip_rate))),
        1e-9,
    )
    traction = above_stress[:, _ALONG_X]
    check(
        "the waves of a tilted pair of media continue the traction, -H V",
        float(
            max(
                np.max(np.abs(traction - below_stress[:, _ALONG_X])),
                np.max(np.abs(traction + slip_rate @ problem.impedance)),
            )
            / np.max(np.abs(traction))
        ),
        1e-9,
    )
    # the one-sided differences at the ends of the series are left out
    inner = slice(2, -2)
    for x in (-17.3, 13.1):
        medium = left if x < 0.0 else right
        velocity, stress = problem.radiated(series, x)
        _, ahead = problem.radiated(series, x + offset)
        velocity_ahead, _ = problem.radiated(series, x + offset)
        velocity_behind, behind = problem.radiated(series, x - offset)
        velocity_gradient = (velocity_ahead - velocity_behind) / (2.0 * offset)
        traction_gradient = (ahead - behind)[:, _ALONG_X] / (2.0 * offset)
        rate = np.gradient(stress, step, axis=0)[inner]
        momentum = medium.rho * np.gradient(velocity, step, axis=0)[inner]
        momentum -= traction_gradient[inner]
        hooke = rate - velocity_gradient[inner] @ medium.stiffness[:, _ALONG_X].T
        check(
            f"the waves at x = {x} solve the wave equation",
            float(
                max(
                    np.max(np.abs(momentum)) / np.max(np.abs(traction_gradient)),
                    np.max(np.abs(hooke)) / np.max(np.abs(rate)),
                )
            ),
            1e-5,
        )

    # slip weakening in the tilted medium: the traction the slip leaves runs
    # along the slip and is the strength under the normal stress it leaves
    fault = {
        "s_xx": -120.0e6,
        "s_xy": 70.0e6,
        "s_xz": 30.0e6,
        "nuc_xy": 20.0e6,
        "mu_s": 0.677,
        "mu_d": 0.525,
        "d_c": 0.08,
    }
    problem = Problem({"FL": 16, "t_0": 0.0625, "XRef": -1.0e5}, fault, left, left)
    series = problem.solve(0.09, 1e-5)
    slipping = series["slip-rate"] > 1e-3 * np.max(series["slip-rate"])
    slip_rate = np.array([series["slip-rate-y"], series["slip-rate-z"]]).T[slipping]
    load = np.array([problem.loading(t) for t in series["time"][slipping]])
    left_traction = load[:, 1:] - slip_rate @ problem.eta_shear
    normal = load[:, 0] - slip_rate @ problem.eta_normal
    strength = -normal * series["friction"][slipping]
    cross = (
        left_traction[:, 0] * slip_rate[:, 1] - left_traction[:, 1] * slip_rate[:, 0]
    )
    check(
        "the slip of the tilted medium runs along the traction it leaves",
        float(
            np.max(
                np.abs(cross)
                / (
                    np.linalg.norm(left_traction, axis=1)
                    * np.linalg.norm(slip_rate, axis=1)
                )
            )
        ),
        1e-12,
    )
    check(
        "the traction of the tilted medium is the strength",
        float(
            max(
                np.max(np.abs(np.linalg.norm(left_traction, axis=1) - strength)),
                np.max(np.abs(normal - series["normal-stress"][slipping])),
            )
            / np.max(strength)
        ),
        1e-12,
    )
    return failures


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("self-test", help="check the integration against closed forms")
    args = parser.parse_args(argv)
    if args.command == "self-test":
        failures = self_test()
        if failures:
            print(f"{len(failures)} of the checks failed", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
