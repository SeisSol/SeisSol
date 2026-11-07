# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Antiplane ruptures by a spectral boundary integral method.

This is the exact solution of the ruptures on fault2d-periodic. The fault is
the plane x = 0 between two halves of the same isotropic solid, of the shear
modulus mu and the shear wave speed c, and slips along z only, uniformly
along z. The solid is periodic along y over the length L, so that the slip
delta(y, t) is a Fourier series in y, and each of its modes, of the wave
number k, sends out the antiplane wave of that mode on both sides. In the
Laplace domain, the traction it leaves on the fault is

    tau_k = -(mu / 2) gamma delta_k,    gamma = sqrt(s^2 / c^2 + k^2),

which splits into the radiation damping -eta V, eta = mu / (2 c), and the
stress that the slip of the past transfers along the fault,

    f_k(t) = -(mu |k| / 2) [delta_k(t) - int_0^t C(|k| c t') V_k(t - t') dt'],

C(T) = int_T^inf J1(u) / u du (Perrin, Rice and Zheng, 1995; Lapusta et al.,
2000). Where the solid is periodic along x as well, over the width W, the
fault repeats at x = n W, and gamma becomes gamma coth(gamma W / 2): every
image adds, from the time tau_n = n W / c on,

    -(mu / c) [V_k(t - tau_n) + int_tau_n^t G_n'(t') V_k(t - t') dt'
               + a^2 int_tau_n^t G_n(t') delta_k(t - t') dt'],

with a = c |k| and G_n(t') = J0(a sqrt(t'^2 - tau_n^2)), whose Laplace
transform is exp(-tau_n c gamma) / (c gamma): the image is the time domain
form of -(mu / c) (s^2 + a^2) times it, applied to delta_k.

Every point of the fault then solves the problem of fault1d, with the
traction tau_load + f and the images in place of tau_load, here with linear
slip weakening. The slip rate is linear between the steps of the
integration, and the slip its integral; the images take the slip as linear
between the steps as well, which differs from that integral by the square of
the step. A step solves the friction law together with the transfer, which
depends on the slip rate of the step through delta_k(t) and the first
weight of the convolution, by fixed-point iteration.

Where a rupture front arrives, the slip rate rises as the square root of the
time since, and its largest deviation converges as the square root of the
spacing of the points. For fault2d-antiplane, with the two steps per time
step of the case that the comparison takes, the slip rate with 1024 points
along y deviates from the one with 2048 by up to 2.5e-2 of its largest
value and the slip by up to 7e-5, while half the step changes them by 5e-4
and 2e-6. The self-test holds the method against a uniform slip across
periodic faults, the finite differences of a mode of a given slip, the
quadrature of the images of a constant slip rate, the images of a mode that
decays within a fraction of the period, and the forced rupture front of
fault1d.Front.

The integration needs scipy, for the Bessel functions.
"""

import math
from pathlib import Path

import fault1d
import numpy as np

# Gauss-Legendre nodes and weights on [0, 1], for the weights of the kernels
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(6)
_NODES = 0.5 * (_NODES + 1.0)
_WEIGHTS = 0.5 * _WEIGHTS


def _transfer_kernel(argument):
    """C(T) = int_T^inf J1(u) / u du = 1 - int_0^T J0(u) du + J1(T)."""
    from scipy import special

    integral, _ = special.itj0y0(argument)
    return 1.0 - integral + special.j1(argument)


def _image_weights(a, delay, step, count):
    """The weights of the image at the delay tau, per wave number (a = c |k|)
    and node j step, j = 0 .. count - 1: those of the slip rate, from
    G'(t') = -a J1(a R) t' / R, and those of the slip, from
    a^2 G(t') = a^2 J0(a R), R = sqrt(t'^2 - tau^2), each the integral of the
    kernel against the hat function of the node over t' >= tau, and each also
    against the rising half of the hat alone, which the node at t' = t takes.

    The integration runs over R, in which the kernels oscillate evenly, as
    G'(t') dt' = -a J1(a R) dR and G(t') dt' = J0(a R) R / t' dR, on pieces
    of at most a radian of the fastest mode: right after tau, a step of t'
    spans many of them.
    """
    from scipy import special

    rate_weights = np.zeros((a.size, count))
    slip_weights = np.zeros((a.size, count))
    rate_rising = np.zeros((a.size, count))
    slip_rising = np.zeros((a.size, count))
    nodes = np.arange(int(math.floor(delay / step)), count - 1)
    if nodes.size == 0:
        return rate_weights, slip_weights, rate_rising, slip_rising
    lower = nodes * step
    upper = lower + step
    radius_low = np.sqrt(np.maximum(np.maximum(lower, delay) ** 2 - delay**2, 0.0))
    radius_high = np.sqrt(upper**2 - delay**2)
    pieces = np.ceil(float(np.max(a)) * (radius_high - radius_low)).astype(int)
    pieces = np.maximum(pieces, 1)
    owner = np.repeat(np.arange(nodes.size), pieces)
    within = np.arange(owner.size) - np.repeat(np.cumsum(pieces) - pieces, pieces)
    width = ((radius_high - radius_low) / pieces)[owner]
    radius = (radius_low[owner] + within * width)[:, None] + _NODES * width[:, None]
    weight = _WEIGHTS * width[:, None]
    time = np.sqrt(radius * radius + delay * delay)
    rising = (time - lower[owner, None]) / step
    falling = (upper[owner, None] - time) / step
    # in chunks of pieces, which bounds the memory of the modes times the points
    for chunk in range(0, owner.size, 256):
        part = slice(chunk, chunk + 256)
        x = a[:, None, None] * radius[None, part]
        rate = -a[:, None, None] * special.j1(x) * weight[None, part]
        slip = (a * a)[:, None, None] * special.j0(x)
        slip = slip * (weight * radius / time)[None, part]
        columns = nodes[owner[part]]
        for weights, halves, kernel in (
            (rate_weights, rate_rising, rate),
            (slip_weights, slip_rising, slip),
        ):
            up = np.einsum("kpg,pg->kp", kernel, rising[part])
            np.add.at(weights, (slice(None), columns + 1), up)
            np.add.at(halves, (slice(None), columns + 1), up)
            np.add.at(
                weights,
                (slice(None), columns),
                np.einsum("kpg,pg->kp", kernel, falling[part]),
            )
    return rate_weights, slip_weights, rate_rising, slip_rising


def _hat_halves(kernel, step, count):
    """The integrals of kernel(t') against the hat functions of the nodes
    j step, j = 0 .. count - 1: their rising halves, below the node, and
    their falling halves, above it, each (wave numbers, nodes)."""
    lower = np.arange(count - 1) * step
    upper = lower + step
    rising, falling = 0.0, 0.0
    for node, weight in zip(_NODES, _WEIGHTS):
        t = lower + node * step
        values = kernel(t) * (weight * step)
        rising = rising + values * (t - lower) / step
        falling = falling + values * (upper - t) / step
    shape = np.shape(rising)[:-1] + (count,)
    up, down = np.zeros(shape), np.zeros(shape)
    up[..., 1:] = rising
    down[..., :-1] = falling
    return up, down


def _smooth_step(t, rise):
    """fault1d.smooth_step for arrays."""
    t, rise = np.broadcast_arrays(
        np.asarray(t, dtype=float), np.asarray(rise, dtype=float)
    )
    inside = (t > 0.0) & (t < rise) & (rise > 0.0)
    safe_t = np.where(inside, t, 0.5)
    safe_rise = np.where(inside, rise, 1.0)
    value = np.exp((safe_t - safe_rise) ** 2 / (safe_t * (safe_t - 2.0 * safe_rise)))
    return np.where(inside, value, np.where(t >= rise, 1.0, 0.0))


def fault_along(model, dynamic_rupture, y):
    """The parameters of the fault at the points (0, y, 0) as the solver
    takes them: per easi parameter an array along y, and the parameters of
    the parameter file, with lowercase names, under the key None."""
    import easimodel

    rows = [easimodel.evaluate_at(model, (0.0, float(value), 0.0)) for value in y]
    if any(row is None for row in rows):
        raise fault1d.Unsupported("the fault is not given everywhere along y")
    fault = {name: np.array([row[name] for row in rows]) for name in rows[0]}
    fault[None] = {key.lower(): value for key, value in dynamic_rupture.items()}
    return fault


class Antiplane:
    """The spectral boundary integral method of an antiplane rupture on the
    plane x = 0.

    `medium` is the isotropic fault1d.Medium of both sides, `length` the
    period L along y, `width` the period W along x, or None for a solid that
    extends to infinity along x, and `count` the number of points along y.
    """

    def __init__(self, medium, length, width=None, count=512):
        if not medium.is_isotropic():
            raise fault1d.Unsupported("the spectral method takes an isotropic solid")
        self.mu = float(medium.stiffness[3, 3])
        self.speed = math.sqrt(self.mu / medium.rho)
        self.eta = 0.5 * self.mu / self.speed
        self.length, self.width, self.count = length, width, count
        self.y = -0.5 * length + np.arange(count) * (length / count)
        self.k = 2.0 * math.pi / length * np.arange(count // 2 + 1)

    def _weights(self, step, steps):
        a = self.speed * self.k[:, None]
        rising, falling = _hat_halves(
            lambda t: _transfer_kernel(a * t[None, :]), step, steps + 1
        )
        # the hat of the node at t' = t counts with its rising half only
        self.direct, self.direct_end = rising + falling, rising
        self.images = []
        if self.width is None:
            return
        if not self.width / self.speed > step:
            raise fault1d.Unsupported(
                "the step is not shorter than the time a wave takes across the period"
            )
        n = 1
        while n * self.width / self.speed < steps * step:
            delay = n * self.width / self.speed
            self.images.append(
                (delay, *_image_weights(a[:, 0], delay, step, steps + 1))
            )
            n += 1

    def _memory(self, rate_k, slip_k, n, step):
        """What the slip up to t_n transfers to t_{n+1}, per mode: the
        convolution with the transfer kernel without the weight of t' = 0,
        which the slip rate of t_{n+1} takes, and the images."""
        t = (n + 1) * step
        past = np.einsum("km,km->k", self.direct[:, 1 : n + 1], rate_k[:, n:0:-1])
        past = past + self.direct_end[:, n + 1] * rate_k[:, 0]
        history = np.zeros(self.k.size, dtype=complex)
        for delay, rate_weights, slip_weights, rate_end, slip_end in self.images:
            if delay >= t:
                break
            # V at t - delay, linear between the steps
            position = (t - delay) / step
            lower = min(int(math.floor(position)), n)
            fraction = position - lower
            upper = min(lower + 1, n)
            delayed_rate = (1.0 - fraction) * rate_k[:, lower]
            delayed_rate = delayed_rate + fraction * rate_k[:, upper]
            # the hat of the node at t' = t counts with its rising half only
            convolution = (
                np.einsum("km,km->k", rate_weights[:, 1 : n + 1], rate_k[:, n:0:-1])
                + rate_end[:, n + 1] * rate_k[:, 0]
                + np.einsum("km,km->k", slip_weights[:, 1 : n + 1], slip_k[:, n:0:-1])
                + slip_end[:, n + 1] * slip_k[:, 0]
            )
            history = history - (self.mu / self.speed) * (delayed_rate + convolution)
        return past, history

    def _transfer(self, slip_k, rate_k, past, history):
        """The stress the slip transfers at t_{n+1}, per mode, from the slip
        and the slip rate of t_{n+1} and what _memory gives."""
        damping = 0.5 * self.mu * self.k
        return -damping * (slip_k - (self.direct[:, 0] * rate_k + past)) + history

    def transfer(self, rate_k, slip_k, step):
        """The stress that a given slip transfers along the fault, per mode
        and step, from its slip rate and its slip per mode and step; it starts
        from rest."""
        steps = rate_k.shape[1] - 1
        self._weights(step, steps)
        result = np.zeros_like(rate_k)
        for n in range(steps):
            past, history = self._memory(rate_k, slip_k, n, step)
            result[:, n + 1] = self._transfer(
                slip_k[:, n + 1], rate_k[:, n + 1], past, history
            )
        return result

    def solve(self, fault, end, step, iterations=50):
        """Integrate up to `end` with steps of `step`.

        `fault` gives the parameters of the fault at the points along y, as
        fault_along does. Returns the times and, per time and point, the slip
        rate along z, the slip, the shear traction along z and the normal
        stress, which an antiplane slip leaves as the loading gives it.
        """
        steps = int(round(end / step))
        self._weights(step, steps)
        friction = _SlipWeakening(fault)
        loading = _Loading(fault)

        modes = self.k.size
        rate_k = np.zeros((modes, steps + 1), dtype=complex)
        slip_k = np.zeros((modes, steps + 1), dtype=complex)
        rate = np.zeros((steps + 1, self.count))
        slip = np.zeros((steps + 1, self.count))
        traction = np.zeros((steps + 1, self.count))
        normal_stress = np.zeros((steps + 1, self.count))
        accumulated = np.zeros(self.count)
        time = np.arange(steps + 1) * step

        shear, normal = loading(0.0)
        rate[0] = friction.rate(shear, normal, accumulated, 0.0, self.eta)
        rate_k[:, 0] = np.fft.rfft(rate[0])
        traction[0] = shear - self.eta * rate[0]
        normal_stress[0] = normal

        for n in range(steps):
            t = time[n + 1]
            past, history = self._memory(rate_k, slip_k, n, step)
            shear, normal = loading(t)
            guess = rate[n]
            for _ in range(iterations):
                new_slip = slip[n] + 0.5 * step * (rate[n] + guess)
                new_accumulated = accumulated + 0.5 * step * (
                    np.abs(rate[n]) + np.abs(guess)
                )
                transfer_k = self._transfer(
                    np.fft.rfft(new_slip), np.fft.rfft(guess), past, history
                )
                transfer = np.fft.irfft(transfer_k, self.count)
                update = friction.rate(
                    shear + transfer, normal, new_accumulated, t, self.eta
                )
                change = np.max(np.abs(update - guess))
                guess = update
                if change <= 1e-12 * max(1.0, np.max(np.abs(guess))):
                    break
            else:
                raise fault1d.Unsupported(f"the step at t = {t:.6g} does not settle")
            rate[n + 1] = guess
            slip[n + 1] = slip[n] + 0.5 * step * (rate[n] + guess)
            accumulated = accumulated + 0.5 * step * (np.abs(rate[n]) + np.abs(guess))
            rate_k[:, n + 1] = np.fft.rfft(rate[n + 1])
            slip_k[:, n + 1] = np.fft.rfft(slip[n + 1])
            traction[n + 1] = shear + transfer - self.eta * rate[n + 1]
            normal_stress[n + 1] = normal
        return {
            "time": time,
            "slip-rate": rate,
            "slip": slip,
            "traction": traction,
            "normal-stress": normal_stress,
        }


class _Loading:
    """The traction along z and the normal stress on the welded fault, per
    point: the initial stress and the nucleations, as fault1d.Loading."""

    def __init__(self, fault):
        parameters = fault[None]
        size = fault["s_xx"].shape

        def field(name):
            values = fault.get(name)
            return np.zeros(size) if values is None else values

        tractions = ("T_n", "T_s", "T_d", "Tnuc_n", "Tnuc_s", "Tnuc_d")
        if any(name in fault for name in tractions):
            raise fault1d.Unsupported(
                "tractions in fault coordinates are not supported"
            )
        if np.any(field("s_xy") != 0.0):
            raise fault1d.Unsupported("an antiplane rupture takes no traction along y")
        self.normal = field("s_xx")
        self.shear = field("s_xz")
        self.sources = []
        count = int(parameters.get("nucleationcount", 1))
        for index in range(count):
            suffix = "" if index == 0 else str(index + 1)
            prefix = f"nuc{suffix}_"
            if not any(f"{prefix}{c}" in fault for c in ("xx", "xy", "xz")):
                continue
            if np.any(field(f"{prefix}xy") != 0.0):
                raise fault1d.Unsupported(
                    "an antiplane rupture takes no traction along y"
                )
            rise = fault.get(f"nuc{suffix}_rise_time")
            if rise is None:
                rise = float(parameters.get(f"t{suffix}_0", 0.0))
            onset = fault.get(f"nuc{suffix}_onset")
            if onset is None:
                onset = float(parameters.get(f"s{suffix}_0", 0.0))
            self.sources.append(
                (field(f"{prefix}xx"), field(f"{prefix}xz"), onset, rise)
            )

    def __call__(self, t):
        normal, shear = self.normal, self.shear
        for dn, dz, onset, rise in self.sources:
            fraction = _smooth_step(t - onset, rise)
            normal = normal + fraction * dn
            shear = shear + fraction * dz
        return shear, normal


class _SlipWeakening:
    """Linear slip weakening (FL=16, or FL=2), with or without a forced
    rupture time, as fault1d.SlipWeakening, per point."""

    def __init__(self, fault):
        parameters = fault[None]
        law = int(parameters.get("fl", 16))
        if law not in (2, 16):
            raise fault1d.Unsupported(
                f"the spectral method takes linear slip weakening, not FL={law}"
            )
        if float(parameters.get("lsw_healingthreshold", -1.0)) >= 0.0:
            raise fault1d.Unsupported("healing is not supported")
        self.mu_s, self.mu_d, self.d_c = fault["mu_s"], fault["mu_d"], fault["d_c"]
        self.cohesion = fault.get("cohesion", 0.0)
        self.forced = None
        if "forced_rupture_time" in fault:
            self.forced = (
                fault["forced_rupture_time"],
                float(parameters.get("t_0", 0.0)),
            )

    def rate(self, shear, normal, accumulated, t, eta):
        weakening = np.minimum(accumulated / self.d_c, 1.0)
        if self.forced is not None:
            onset, rise = self.forced
            if rise > 0.0:
                weakening = np.maximum(weakening, np.clip((t - onset) / rise, 0.0, 1.0))
            else:
                weakening = np.where(t >= onset, 1.0, weakening)
        friction = self.mu_s - (self.mu_s - self.mu_d) * weakening
        strength = -self.cohesion - np.minimum(normal, 0.0) * friction
        return np.sign(shear) * np.maximum(np.abs(shear) - strength, 0.0) / eta


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------


def _klein_gordon_traction(medium, k, half_width, slip, end, step):
    """The traction mu dU/dx (0, t) of the mode U(x, t) cos(k y) of an
    antiplane slip slip(t) cos(k y) on the faults x = n W of a solid periodic
    across them, from U_tt = c^2 (U_xx - k^2 U) on 0 < x < W / 2, with
    U = slip / 2 on the fault and U = 0 halfway to the next, by finite
    differences of the second order at `step`, and at the Courant number 1,
    or just below, at which they carry the waves along x without error."""
    mu = float(medium.stiffness[3, 3])
    speed = math.sqrt(mu / medium.rho)
    points = int(math.floor(half_width / (speed * step)))
    spacing = half_width / points
    courant = (speed * step / spacing) ** 2
    previous = np.zeros(points + 1)
    current = np.zeros(points + 1)
    traction = [0.0]
    for n in range(1, int(round(end / step)) + 1):
        following = np.empty_like(current)
        following[1:-1] = (
            2.0 * current[1:-1]
            - previous[1:-1]
            + courant * (current[2:] - 2.0 * current[1:-1] + current[:-2])
            - (speed * step * k) ** 2 * current[1:-1]
        )
        following[0] = 0.5 * slip(n * step)
        following[-1] = 0.0
        previous, current = current, following
        traction.append(
            mu * (-3.0 * current[0] + 4.0 * current[1] - current[2]) / (2.0 * spacing)
        )
    return np.array(traction)


def self_test():
    """Hold the method against what is known in closed form or computed
    otherwise.

    Returns the failures, as text.
    """
    import easimodel

    failures = []
    medium = fault1d.Medium.isotropic(2700.0, 3.24e10, 3.24e10)

    def check(label, deviation, tolerance):
        status = "ok" if deviation <= tolerance else "FAILED"
        print(
            f"[{status}] {label}: deviation {deviation:.3g}, tolerance {tolerance:.3g}"
        )
        if deviation > tolerance:
            failures.append(label)

    # a uniform slip under a constant friction across the faults x = n W: the
    # mode k = 0 alone, which transfers nothing along the fault; the wave each
    # side sends out arrives at the next fault after T = W / c, from both of
    # its sides, so that eta V(t) = tau(t) - strength - 2 eta sum V(t - n T)
    width, end, step = 400.0 / 3.0, 0.12, 2e-5
    solver = Antiplane(medium, 800.0, width, 8)
    loading = {"s_xx": -120.0e6, "s_xz": 70.0e6, "nuc_xz": 10.0e6}
    friction = {"mu_s": 0.6, "mu_d": 0.6, "d_c": 1.0}
    fault = {name: np.full(8, value) for name, value in {**loading, **friction}.items()}
    fault[None] = {"t_0": 0.02}
    series = solver.solve(fault, end, step)
    delay = width / solver.speed

    def rate(t):
        if t <= 0.0:
            return 0.0
        shear = 70.0e6 + 10.0e6 * fault1d.smooth_step(t, 0.02)
        echoes = sum(rate(t - n * delay) for n in range(1, int(t / delay) + 1))
        return max(0.0, (shear - 0.6 * 120.0e6) / solver.eta - 2.0 * echoes)

    # the slip rate is linear between the steps, which it is not where the
    # fault starts or stops slipping, so that the echoes of these times are
    # off by a fraction of a step there
    exact = np.array([rate(t) for t in series["time"]])
    check(
        "a uniform slip across periodic faults",
        float(np.max(np.abs(series["slip-rate"][:, 0] - exact)) / np.max(exact)),
        3e-3,
    )

    # a mode of a prescribed slip, against the finite differences of its wave
    # equation across the faults: the transfer kernel and the images
    duration = 0.01

    def slip(t):
        return (t / duration) ** 3 * math.exp(-t / duration)

    def slip_rate(t):
        return (
            (3.0 - t / duration)
            * (t / duration) ** 2
            * math.exp(-t / duration)
            / duration
        )

    step = 6.103515625e-05 / 2.0
    time = np.arange(int(round(end / step)) + 1) * step
    solver = Antiplane(medium, 800.0, width, 16)
    mode = 5
    rate_k = np.zeros((solver.k.size, time.size), dtype=complex)
    slip_k = np.zeros_like(rate_k)
    rate_k[mode] = [slip_rate(t) for t in time]
    slip_k[mode] = [slip(t) for t in time]
    computed = solver.transfer(rate_k, slip_k, step)[mode].real
    computed = computed - solver.eta * rate_k[mode].real
    exact = _klein_gordon_traction(medium, solver.k[mode], 0.5 * width, slip, end, step)
    check(
        f"the traction of the mode k = {solver.k[mode]:.4g} / m across periodic faults",
        float(np.max(np.abs(computed - exact)) / np.max(np.abs(exact))),
        1e-4,
    )

    # a mode that decays within a fraction of the width does not reach the
    # next fault, but for the part of the slip faster than c |k|, which is
    # tiny here: its images transfer next to nothing, however fast their
    # kernels oscillate right after they arrive
    solver = Antiplane(medium, 12.5, width, 16)
    rate_k = np.zeros((solver.k.size, time.size), dtype=complex)
    slip_k = np.zeros_like(rate_k)
    mode = 4
    rate_k[mode] = [slip_rate(t) for t in time]
    slip_k[mode] = [slip(t) for t in time]
    periodic = solver.transfer(rate_k, slip_k, step)[mode]
    solver.width = None
    alone = solver.transfer(rate_k, slip_k, step)[mode]
    check(
        f"the images of the mode k = {solver.k[mode]:.4g} / m",
        float(np.max(np.abs(periodic - alone)) / np.max(np.abs(alone))),
        1e-5,
    )

    # the image of a slip rate of 1 m/s from t = 0 on, of a mode whose kernels
    # oscillate some thirty times in the first step after it arrives, against
    # the quadrature of -(mu / c) [J0(a R(t)) + a^2 int G(t') (t - t') dt']
    from scipy import integrate, special

    width = 800.0 / 3.0
    solver = Antiplane(medium, 4.0 * math.pi, width, 16)
    mode = 8
    rate_k = np.zeros((solver.k.size, time.size), dtype=complex)
    slip_k = np.zeros_like(rate_k)
    rate_k[mode] = 1.0
    slip_k[mode] = time
    periodic = solver.transfer(rate_k, slip_k, step)[mode].real
    solver.width = None
    alone = solver.transfer(rate_k, slip_k, step)[mode].real
    a, delay = solver.speed * solver.k[mode], width / solver.speed
    deviation = 0.0
    for index in range(int(delay / step) + 1, time.size, 61):
        t = time[index]
        integral, _ = integrate.quad(
            lambda u, t=t: special.j0(a * math.sqrt(max(u * u - delay * delay, 0.0)))
            * (t - u),
            delay,
            t,
            limit=2000,
            epsabs=1e-10 / (a * a),
            epsrel=1e-10,
        )
        image = special.j0(a * math.sqrt(t * t - delay * delay)) + a * a * integral
        deviation = max(
            deviation, abs((periodic - alone)[index] * solver.speed / solver.mu + image)
        )
    check(
        f"the image of a constant slip rate of the mode k = {solver.k[mode]:.4g} / m",
        deviation,
        1e-8,
    )

    # the forced rupture front of fault2d-supershear.yaml, which fault1d.Front
    # solves, up to the time to which its solution holds
    rupture = {"FL": 16, "t_0": 0.02, "XRef": -1.0e5}
    model = easimodel.load(Path(__file__).parent / "cases" / "fault2d-supershear.yaml")
    solver = Antiplane(medium, 1200.0, None, 256)
    series = solver.solve(fault_along(model, rupture, solver.y), end, step)
    front = fault1d.Front(rupture, model, medium, medium, 1200.0, math.inf)
    for y in (300.0, -356.25):
        point = (0.0, y, 0.0)
        exact = front.problem_at(point).solve(end, step)["slip-rate"]
        held = series["time"] <= front.window(point)
        computed = series["slip-rate"][:, int(round((y + 600.0) / (1200.0 / 256)))]
        check(
            f"the forced rupture front at y = {y}",
            float(np.max(np.abs(computed[held] - exact[held])) / np.max(exact[held])),
            3e-3,
        )
    return failures


def main(argv=None):
    import argparse
    import importlib.util
    import sys

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "self-test", help="check the method against independent solutions"
    )
    args = parser.parse_args(argv)
    if args.command == "self-test":
        if importlib.util.find_spec("scipy") is None:
            print("skipped: the method needs scipy, which is not installed")
            return 77
        failures = self_test()
        if failures:
            print(f"{len(failures)} of the checks failed", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
