# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Antiplane and in-plane ruptures by a spectral boundary integral method.

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
slip weakening (FL=16, or FL=2), its proxy of thermal pressurization
(FL=1058) and its Prakash-Clifton regularization (FL=6), each with the
instantaneous healing of SeisSol, one of the rate-and-state laws FL=3, 4, 7
and 103, with thermal pressurization or without, or an imposed slip rate,
FL=33 and 34. The slip rate is linear between the steps of the integration,
and the slip its integral; the images take the slip as linear between the
steps as well, which differs from that integral by the square of the step. A
step advances the state of the law, the accumulated slip of slip weakening,
with FL=6 also the regularized strength, or psi of rate and state, in closed
form at the mean slip rate of the step, as SeisSol does, and solves the
friction law together with the transfer, which depends on the slip rate of
the step through delta_k(t) and the first weight of the convolution, by
fixed-point iteration; with thermal pressurization, together with the pore
pressure, from the frictional heating of fault1d.ThermalPressurization, mean
over the step. Besides the slip rate, the slip and the traction, the solution
holds what the fault receivers take from their Gauss point: the friction
coefficient, the accumulated slip, the largest slip rate so far, the state,
the pore pressure and the temperature, and the times of the rupture and of
the dynamic stress.

Off the fault, at 0 < x < W, the displacement of a mode is

    u_k = (delta_k / 2) sinh(gamma (W / 2 - x)) / sinh(gamma W / 2)
        = (delta_k / 2) sum_n [exp(-gamma (x + n W)) - exp(-gamma ((n + 1) W - x))],

what the fault and its images send out over the distances d = x + n W and
(n + 1) W - x, n >= 0; where the solid extends to infinity along x, it is
sign(x) (delta_k / 2) exp(-gamma |x|). exp(-gamma d) is the propagator
P_d f = f(t - tau) - int_tau^t a tau J1(a R) / R f(t - t') dt' of the delay
tau = d / c, R = sqrt(t'^2 - tau^2), and gamma exp(-gamma d) is 1 / c times
the kernel of an image of that delay: the velocity is the sum of the
+-P_d V_k / 2, the stress mu du/dy is i k mu u_k, and mu du/dx is
-(mu / 2 c) times the sum of the images of all distances.

An in-plane rupture (Inplane) slips along y, and sends out P waves, of the
speed c_p = r c, as well as S waves. With the potentials of the two, a slip
of the mode leaves at the distance d of the fault (Geubelle and Rice, 1995)

    sigma_xy = -(mu / 2) (c^2 / s^2) [(k^2 + beta^2)^2 / beta e^(-beta d)
                                      - 4 k^2 alpha e^(-alpha d)] delta_k,

beta = sqrt(s^2 / c^2 + k^2) and alpha = sqrt(s^2 / c_p^2 + k^2), whereas in
a solid of one medium, the slip leaves the normal stress on the fault as it
is. On the fault, d = 0, this is the radiation damping and

    f_k(t) = -mu |k| kappa [delta_k(t) - int_0^t C_II(|k| c t') V_k(t - t') dt'],

kappa = 1 - 1 / r^2, with the kernel C_II of _inplane_transfer_kernel, whose
static part is the stiffness mu |k| / (2 (1 - nu)) of mode II. Taken apart
into the waves of each speed, the factor c^2 / s^2 leaves parts that grow
with t^2 and cancel only together; as a sum over the speeds c' between c and
c_p instead, the difference of the two waves is a near field, which acts
between their arrivals only and is bounded:

    sigma_xy = -(mu / 2) [(s / c^2) e^(-beta d) / beta V_k
                          + 4 k^2 c^2 L(Q_d) delta_k],

L(Q_d) = (beta e^(-beta d) - alpha e^(-alpha d)) / s^2, the Laplace transform
of Q_d(t) = t [F(w_P) - F(w_S)] - t / d while only the P wave has arrived,
with F' = [w J0(|k| w) + |k| d^2 J1(|k| w)] / (w^2 + d^2)^(3/2), w_P =
sqrt(c_p^2 t^2 - d^2), and w_S = sqrt(c^2 t^2 - d^2) once the S wave has
arrived, 0 before. An image adds the first term as an antiplane image does,
of the S waves, and 4 k^2 c^3 Q_d in place of a^2 G_n, from the arrival of
its P wave on. Off the fault, the velocity along y is the sum of the
+-P_d V_k / 2 of the S waves and of k^2 c^2 E_d V_k, with E_d the near field
of (e^(-beta d) - e^(-alpha d)) / s^2, -d t [F(w_P) - F(w_S)] with
F' = w J0(|k| w) / (w^2 + d^2)^(3/2), and the displacement accordingly; the
velocity along x is i k [c^2 Q_d V_k - (c / 2) G_d V_k], with G_d the G_n of
the distance d, the same from both sides; the divergence is
i k (c^2 / c_p^2) P'_d delta_k, with the propagator P' of the P waves;
sigma_xy is the sum of the images of all distances, without their factor 2;
and sigma_zz = lambda div u, sigma_yy = sigma_zz + 2 mu i k u_y and
sigma_xx = 2 (lambda + mu) div u - sigma_yy. What is odd in x, u_y and the
normal stresses, takes the sign of the side.

Where a rupture front arrives, the slip rate rises as the square root of the
time since, and its largest deviation converges as the square root of the
spacing of the points. For fault2d-antiplane, with the two steps per time
step of the case that the comparison takes, the slip rate with 1024 points
along y deviates from the one with 2048 by up to 2.5e-2 of its largest
value and the slip by up to 7e-5, while half the step changes them by 5e-4
and 2e-6; with the other laws, 2048 points change the slip rate by 4.5e-3 to
2.2e-2 and the slip by up to 2.8e-5, and half the step, with the aging law,
by 1.6e-4 and 7e-7. The self-test holds the method against a uniform slip
across periodic faults, with healing and without, the finite differences of
a mode of a given slip, on the fault and off it, antiplane and in plane, the
quadrature of the images of a constant slip rate, the images of a mode that
decays within a fraction of the period, the forced rupture front of
fault1d.Front, and the uniform slip of every other law against
fault1d.Problem, with what the Gauss points report.

The integration needs scipy, for the Bessel functions.
"""

import math
import re
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


def _delayed_weights(a, delay, step, count, kernels):
    """The weights of kernels that start at the delay tau, per wave number
    (a = c |k|) and node j step, j = 0 .. count - 1: the integral of each
    kernel against the hat function of the node over t' >= tau, and against
    the rising half of the hat alone, which the node at t' = t takes.

    A kernel is a pair (f, measure): in R = sqrt(t'^2 - tau^2), in which the
    kernels oscillate evenly, its integral is that of f(a, a R) over R, each
    point R of the quadrature weighted by measure(weight, R, t'), or by its
    weight where measure is None. The integration runs on pieces of at most
    a radian of the fastest mode: right after tau, a step of t' spans many of
    them.
    """
    weights = [(np.zeros((a.size, count)), np.zeros((a.size, count))) for _ in kernels]
    nodes = np.arange(int(math.floor(delay / step)), count - 1)
    if nodes.size == 0:
        return weights
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
        columns = nodes[owner[part]]
        for (full, halves), (function, measure) in zip(weights, kernels):
            if measure is None:
                kernel = function(a, x) * weight[None, part]
            else:
                kernel = function(a, x)
                kernel = kernel * measure(weight, radius, time)[None, part]
            up = np.einsum("kpg,pg->kp", kernel, rising[part])
            np.add.at(full, (slice(None), columns + 1), up)
            np.add.at(halves, (slice(None), columns + 1), up)
            np.add.at(
                full,
                (slice(None), columns),
                np.einsum("kpg,pg->kp", kernel, falling[part]),
            )
    return weights


def _image_weights(a, delay, step, count):
    """The weights of the image at the delay tau (see _delayed_weights):
    those of the slip rate, from G'(t') = -a J1(a R) t' / R, and those of the
    slip, from a^2 G(t') = a^2 J0(a R), with G'(t') dt' = -a J1(a R) dR and
    G(t') dt' = J0(a R) R / t' dR; each as the weights of the whole hats and
    of their rising halves."""
    from scipy import special

    (rate_weights, rate_rising), (slip_weights, slip_rising) = _delayed_weights(
        a,
        delay,
        step,
        count,
        (
            (lambda a, x: -a[:, None, None] * special.j1(x), None),
            (
                lambda a, x: (a * a)[:, None, None] * special.j0(x),
                lambda weight, radius, time: weight * radius / time,
            ),
        ),
    )
    return rate_weights, slip_weights, rate_rising, slip_rising


def _propagator_weights(a, delay, step, count):
    """The weights of a tau J1(a R) / R, the kernel of the propagator of the
    delay tau (see _delayed_weights), with dt' = R / t' dR: those of the whole
    hats and of their rising halves."""
    from scipy import special

    ((weights, rising),) = _delayed_weights(
        a,
        delay,
        step,
        count,
        (
            (
                lambda a, x: (a * delay)[:, None, None] * special.j1(x),
                lambda weight, radius, time: weight / time,
            ),
        ),
    )
    return weights, rising


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


def _convolve(weights, rising, values):
    """sum_j w_j f(t - j step) per mode and node of time, with the node at
    t' = t weighted by the rising half of its hat; f is zero before t = 0."""
    from scipy import signal

    count = values.shape[1]
    result = np.empty(values.shape, dtype=complex)
    # in chunks of modes, which bounds the memory of the transforms
    for chunk in range(0, values.shape[0], 64):
        part = slice(chunk, chunk + 64)
        result[part] = signal.fftconvolve(weights[part], values[part], axes=1)[
            :, :count
        ]
    return result + (rising - weights) * values[:, :1]


def _delayed(values, delay, step):
    """f(t - tau) per mode and node of time, linear between the nodes and
    zero before t = 0."""
    count = values.shape[1]
    position = np.arange(count) - delay / step
    lower = np.floor(position).astype(int)
    fraction = position - lower
    result = np.zeros(values.shape, dtype=complex)
    inside = position >= 0.0
    low = lower[inside]
    high = np.minimum(low + 1, count - 1)
    result[:, inside] = (1.0 - fraction[inside]) * values[:, low]
    result[:, inside] += fraction[inside] * values[:, high]
    return result


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


def _smooth_step_rate(t, rise):
    """fault1d.smooth_step_rate for arrays."""
    t, rise = np.broadcast_arrays(
        np.asarray(t, dtype=float), np.asarray(rise, dtype=float)
    )
    inside = (t > 0.0) & (t < rise) & (rise > 0.0)
    safe_t = np.where(inside, t, 0.5)
    safe_rise = np.where(inside, rise, 1.0)
    product = safe_t * (safe_t - 2.0 * safe_rise)
    value = np.exp((safe_t - safe_rise) ** 2 / product)
    rate = value * (
        -2.0 * safe_rise * safe_rise * (safe_t - safe_rise) / (product * product)
    )
    return np.where(inside, rate, 0.0)


def _yoffe_rate(t, tau_s, tau_r):
    """fault1d._yoffe_convolution of the rate, for arrays: the pieces between
    the kinks of the triangle, at s = t - 2 tau_s, t - tau_s and t, in
    ascending order of the angle phi of s = tau_r sin^2(phi), each by the
    same Gauss-Legendre rule."""
    t, tau_s, tau_r = np.broadcast_arrays(
        np.asarray(t, dtype=float),
        np.asarray(tau_s, dtype=float),
        np.asarray(tau_r, dtype=float),
    )
    nodes, weights = fault1d._GAUSS_NODES, fault1d._GAUSS_WEIGHTS

    def angle(s):
        return np.arcsin(np.sqrt(np.clip(s / tau_r, 0.0, 1.0)))

    breaks = (
        np.zeros(t.shape),
        angle(t - 2.0 * tau_s),
        angle(t - tau_s),
        angle(t),
        np.full(t.shape, 0.5 * math.pi),
    )
    total = np.zeros(t.shape)
    for low, high in zip(breaks[:-1], breaks[1:]):
        half = 0.5 * (high - low)[..., None]
        phi = 0.5 * (high + low)[..., None] + half * nodes
        u = t[..., None] - tau_r[..., None] * np.sin(phi) ** 2
        width = tau_s[..., None]
        triangle = np.where(
            u <= 0.0,
            0.0,
            np.where(u <= width, u, np.maximum(2.0 * width - u, 0.0)) / (width * width),
        )
        total += np.sum(weights * half * np.cos(phi) ** 2 * triangle, axis=-1)
    return 4.0 / math.pi * total


def _first_crossings(time, values, threshold):
    """fault1d._first_crossing per column of `values`, against a threshold
    per column or one for all, NaN where the values stay below it."""
    threshold = np.broadcast_to(threshold, values.shape[1:])
    above = values >= threshold
    reached = above.any(axis=0)
    index = np.argmax(above, axis=0)
    columns = np.arange(values.shape[1])
    before = values[np.maximum(index - 1, 0), columns]
    after = values[index, columns]
    with np.errstate(divide="ignore", invalid="ignore"):
        fraction = np.where(
            after != before, (threshold - before) / (after - before), 1.0
        )
    previous = time[np.maximum(index - 1, 0)]
    crossing = np.where(
        index == 0, time[0], previous + fraction * (time[index] - previous)
    )
    return np.where(reached, crossing, np.nan)


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
    extends to infinity along x, and `count` the number of points along y,
    which start at y = -L / 2 + `offset`.
    """

    # the direction of the slip, in its y and z components
    axis = "z"
    direction = np.array((0.0, 1.0))

    def __init__(self, medium, length, width=None, count=512, offset=0.0):
        if not medium.is_isotropic():
            raise fault1d.Unsupported("the spectral method takes an isotropic solid")
        self.mu = float(medium.stiffness[3, 3])
        self.speed = math.sqrt(self.mu / medium.rho)
        # the fastest wave, which reaches farthest
        self.fastest = self.speed
        self.eta = 0.5 * self.mu / self.speed
        self.length, self.width, self.count = length, width, count
        self.y = -0.5 * length + offset + np.arange(count) * (length / count)
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
        rate along the axis of the slip, z here, the slip, the shear traction
        along it and the normal stress, which the slip leaves as the loading
        gives it, the
        friction coefficient, the accumulated slip and the largest slip rate
        so far, and with rate and state the state psi; with thermal
        pressurization, the normal stress is the effective one, and the pore
        pressure, as SeisSol states it (-p), and the temperature come along;
        and per point the times the fault ruptures, where the slip rate first
        reaches fault1d.RUPTURE_SLIP_RATE, and weakens (see the weakened
        method of the law), NaN where it does not.
        """
        steps = int(round(end / step))
        self._weights(step, steps)
        friction = _friction_law(fault, self.axis)
        loading = _Loading(fault, self.axis)
        thermal = None
        if int(fault[None].get("thermalpress", 0)):
            if not isinstance(friction, _RateAndState):
                raise fault1d.Unsupported("thermal pressurization takes rate and state")
            thermal = _ThermalPressurization(fault)
            thermal.prepare(step, steps, self.count)

        modes = self.k.size
        rate_k = np.zeros((modes, steps + 1), dtype=complex)
        slip_k = np.zeros((modes, steps + 1), dtype=complex)
        rate = np.zeros((steps + 1, self.count))
        slip = np.zeros((steps + 1, self.count))
        traction = np.zeros((steps + 1, self.count))
        normal_stress = np.zeros((steps + 1, self.count))
        coefficient = np.zeros((steps + 1, self.count))
        # the state of rate and state, which the solution reports
        states = np.zeros((steps + 1, self.count))
        time = np.arange(steps + 1) * step

        # the pore pressure p, which lowers the compression of the normal
        # stress, from its initial value p_0 = -tp_iniPressure on
        pressure = 0.0
        if thermal is not None:
            pressure = -thermal.pore_pressure
            pore = np.full((steps + 1, self.count), thermal.pore_pressure)
            temperature = np.full((steps + 1, self.count), thermal.t0)
        shear, normal = loading(0.0)
        state = friction.initial_state(shear, normal + pressure)
        rate[0] = friction.rate(shear, normal + pressure, state, 0.0, self.eta)
        rate_k[:, 0] = np.fft.rfft(rate[0])
        traction[0] = shear - self.eta * rate[0]
        normal_stress[0] = normal + pressure
        coefficient[0] = friction.friction(state, rate[0], 0.0, normal + pressure)
        reported = isinstance(friction, _RateAndState)
        if reported:
            states[0] = state

        for n in range(steps):
            t = time[n + 1]
            past, history = self._memory(rate_k, slip_k, n, step)
            shear, normal = loading(t)
            # what rounding leaves of the traction, which the slip rate takes
            # divided by eta, bounds how closely the iteration can settle
            settled = 1e-12 * max(1.0, np.max(np.abs(shear)) / self.eta)
            if thermal is not None:
                thermal.begin()
                heated = np.abs(traction[n] * rate[n])

            def heating(guess, transfer):
                # the mean heating of the step, with the traction of the end
                # that the slip rate leaves
                ending = np.abs((shear + transfer - self.eta * guess) * guess)
                return 0.5 * (heated + ending)

            guess = rate[n]
            for _ in range(iterations):
                new_slip = slip[n] + 0.5 * step * (rate[n] + guess)
                new_state = friction.advance(state, rate[n], guess, step, t, normal)
                transfer_k = self._transfer(
                    np.fft.rfft(new_slip), np.fft.rfft(guess), past, history
                )
                transfer = np.fft.irfft(transfer_k, self.count)
                if thermal is not None:
                    _, rise = thermal.changes(heating(guess, transfer))
                    pressure = rise - thermal.pore_pressure
                update = friction.rate(
                    shear + transfer, normal + pressure, new_state, t, self.eta, guess
                )
                change = np.max(np.abs(update - guess))
                guess = update
                if change <= max(settled, 1e-12 * np.max(np.abs(guess))):
                    break
            else:
                raise fault1d.Unsupported(f"the step at t = {t:.6g} does not settle")
            rate[n + 1] = guess
            slip[n + 1] = slip[n] + 0.5 * step * (rate[n] + guess)
            state = friction.advance(state, rate[n], guess, step, t, normal)
            rate_k[:, n + 1] = np.fft.rfft(rate[n + 1])
            slip_k[:, n + 1] = np.fft.rfft(slip[n + 1])
            traction[n + 1] = shear + transfer - self.eta * rate[n + 1]
            if thermal is not None:
                current = heating(guess, transfer)
                warming, rise = thermal.changes(current)
                pressure = rise - thermal.pore_pressure
                thermal.record(current)
                pore[n + 1] = -pressure
                temperature[n + 1] = thermal.t0 + warming
            normal_stress[n + 1] = normal + pressure
            coefficient[n + 1] = friction.friction(
                state, rate[n + 1], t, normal + pressure
            )
            if reported:
                states[n + 1] = state
        speed = np.abs(rate)
        accumulated = np.zeros(rate.shape)
        accumulated[1:] = np.cumsum(0.5 * step * (speed[1:] + speed[:-1]), axis=0)
        solution = {
            "time": time,
            "slip-rate": rate,
            "slip": slip,
            "traction": traction,
            "normal-stress": normal_stress,
            "friction": coefficient,
            "accumulated-slip": accumulated,
            "peak-slip-rate": np.maximum.accumulate(speed, axis=0),
            "rupture-time": _first_crossings(time, speed, fault1d.RUPTURE_SLIP_RATE),
        }
        if reported:
            solution["state"] = states
        if isinstance(friction, _ImposedSlipRate):
            # the traction the fault would carry if it were welded
            solution["locked-traction"] = traction + self.eta * rate
        if thermal is not None:
            solution["pore-pressure"] = pore
            solution["temperature"] = temperature
        solution["dynamic-stress-time"] = friction.weakened(solution)
        return solution

    def _sources(self, x, end):
        """The distances from which the fault and its images reach the point
        at x before `end`, each with the sign of the displacement it sends
        there (see the description of the module)."""
        if self.width is None:
            if x == 0.0:
                raise fault1d.Unsupported("the point lies on the fault")
            return [(abs(x), math.copysign(1.0, x))]
        offset = x % self.width
        if min(offset, self.width - offset) <= 1e-9 * self.width:
            raise fault1d.Unsupported("the point lies on the fault")
        sources = []
        n = 0
        reach = self.fastest * end
        while min(offset, self.width - offset) + n * self.width < reach:
            for distance, sign in (
                (offset + n * self.width, 1.0),
                ((n + 1) * self.width - offset, -1.0),
            ):
                if distance < reach:
                    sources.append((distance, sign))
            n += 1
        return sources

    def radiated(self, solution, x, y):
        """The particle velocity (x, y, z) and the stress (ordered as
        fault1d.VOIGT) that the slip of `solution` leaves at the point (x, y),
        per time of the solution, as the receivers record them: without the
        initial stress. Off the fault, the field is smooth along y, as the
        modes of the front decay away from it, so that its Fourier series
        holds at any y."""
        time = solution["time"]
        step = time[1] - time[0]
        rate_k = np.fft.rfft(solution["slip-rate"], axis=1).T
        slip_k = np.fft.rfft(solution["slip"], axis=1).T
        a = self.speed * self.k
        shape = rate_k.shape
        propagator, propagator_end = np.zeros(shape), np.zeros(shape)
        image_rate, image_rate_end = np.zeros(shape), np.zeros(shape)
        image_slip, image_slip_end = np.zeros(shape), np.zeros(shape)
        velocity_k = np.zeros(shape, dtype=complex)
        displacement_k = np.zeros(shape, dtype=complex)
        arrivals = np.zeros(shape, dtype=complex)
        for distance, sign in self._sources(x, time[-1]):
            delay = distance / self.speed
            if not delay > step:
                raise fault1d.Unsupported(
                    "the point lies closer to the fault than a wave travels in a step"
                )
            weights, rising = _propagator_weights(a, delay, step, time.size)
            propagator += sign * weights
            propagator_end += sign * rising
            rate_weights, slip_weights, rate_end, slip_end = _image_weights(
                a, delay, step, time.size
            )
            image_rate += rate_weights
            image_rate_end += rate_end
            image_slip += slip_weights
            image_slip_end += slip_end
            arrived = _delayed(rate_k, delay, step)
            velocity_k += 0.5 * sign * arrived
            displacement_k += 0.5 * sign * _delayed(slip_k, delay, step)
            arrivals += arrived
        velocity_k -= 0.5 * _convolve(propagator, propagator_end, rate_k)
        displacement_k -= 0.5 * _convolve(propagator, propagator_end, slip_k)
        shear_k = -(0.5 * self.mu / self.speed) * (
            arrivals
            + _convolve(image_rate, image_rate_end, rate_k)
            + _convolve(image_slip, image_slip_end, slip_k)
        )

        # the real Fourier series at y, of the modes of numpy.fft.rfft
        factor = np.full(self.k.size, 2.0)
        factor[0] = 1.0
        if self.count % 2 == 0:
            factor[-1] = 1.0
        factor = factor * np.exp(1j * self.k * (y - self.y[0])) / self.count

        def at_y(values):
            return (factor @ values).real

        velocity = np.zeros((time.size, 3))
        velocity[:, 2] = at_y(velocity_k)
        stress = np.zeros((time.size, 6))
        stress[:, fault1d.VOIGT.index("yz")] = self.mu * at_y(
            1j * self.k[:, None] * displacement_k
        )
        stress[:, fault1d.VOIGT.index("xz")] = at_y(shear_k)
        return velocity, stress


def _inplane_transfer_kernel(argument, ratio):
    """C_II(T), the transfer kernel of an in-plane slip in T = c_s |k| t, at
    the ratio r = c_p / c_s of the wave speeds, from C_II(0) = 1 to 0:

        1 + J1(T) / (2 kappa) - [Ji0(T) - Ji0(r T) / r^2 + T J0(T)
                                 - T J0(r T) / r - T^2 (C(T) - C(r T))] / kappa,

    kappa = 1 - 1 / r^2 and Ji0 the integral of J0 from 0, with the antiplane
    kernel C, whose Laplace transform in p = s / (c_s |k|) is
    (1 / p) [1 - Q(p) / (2 kappa)],
    Q(p) = (2 + p^2)^2 / (p^2 sqrt(p^2 + 1)) - (4 / r) sqrt(p^2 + r^2) / p^2 - p.
    The terms in T and T^2 cancel to an oscillation of the order sqrt(T)."""
    from scipy import special

    integral, _ = special.itj0y0(argument)
    integral_p, _ = special.itj0y0(ratio * argument)
    bracket = (
        integral
        - integral_p / (ratio * ratio)
        + argument * special.j0(argument)
        - argument * special.j0(ratio * argument) / ratio
        - argument
        * argument
        * (_transfer_kernel(argument) - _transfer_kernel(ratio * argument))
    )
    kappa = 1.0 - 1.0 / (ratio * ratio)
    return 1.0 + special.j1(argument) / (2.0 * kappa) - bracket / kappa


def _near_field_weights(k, speeds, distance, step, count, kinds, chunk=48):
    """The weights of the near fields of an in-plane slip at `distance`
    against the hat functions of the nodes j step, j = 0 .. count - 1, whole
    and their rising halves, per wave number k >= 0, with the wave speeds
    (c_s, c_p), by kind, for the `kinds` among "Q" and "E": from the arrival
    of the P wave on, of

        Q_d(t') = t' [F(w_P) - F(w_S)] - t' / d before the S wave arrives,
        F' = [w J0(k w) + k d^2 J1(k w)] / (w^2 + d^2)^(3/2),

    the kernel of (beta e^(-beta d) - alpha e^(-alpha d)) / s^2, for "Q",
    and of E_d(t') = -d t' [F(w_P) - F(w_S)], F' = w J0(k w) /
    (w^2 + d^2)^(3/2), the kernel of (e^(-beta d) - e^(-alpha d)) / s^2, for
    "E" (see the description of the module); w_P = sqrt(c_p^2 t'^2 - d^2),
    and w_S = sqrt(c_s^2 t'^2 - d^2) once the S wave has arrived, 0 before.
    Between the arrivals, the integration runs in R = sqrt(t'^2 - tau_P^2),
    after them in R = sqrt(t'^2 - tau_S^2), on pieces of at most a radian of
    the fastest mode; F between the values of w the pieces take, in order,
    which lie within a fifth of a radian of that mode of each other, where
    three points integrate F' as six do, to 1e-13.
    """
    from scipy import sparse, special

    c_s, c_p = speeds
    tau_p, tau_s = distance / c_p, distance / c_s
    weights = {
        kind: (np.zeros((k.size, count)), np.zeros((k.size, count))) for kind in kinds
    }
    first = int(math.floor(tau_p / step))
    if first >= count - 1:
        return weights
    fastest = max(float(np.max(k)) * c_p, 1.0 / step)
    pieces = []
    for window, (start, stop, delay) in enumerate(
        ((tau_p, tau_s, tau_p), (tau_s, math.inf, tau_s))
    ):
        nodes = np.arange(max(first, int(math.floor(start / step))), count - 1)
        lower = np.maximum(nodes * step, start)
        upper = np.minimum((nodes + 1) * step, stop)
        inside = upper > lower
        nodes, lower, upper = nodes[inside], lower[inside], upper[inside]
        if nodes.size == 0:
            continue
        radius_low = np.sqrt(np.maximum(lower**2 - delay**2, 0.0))
        radius_high = np.sqrt(upper**2 - delay**2)
        parts = np.maximum(np.ceil(fastest * (radius_high - radius_low)), 1).astype(int)
        owner = np.repeat(np.arange(nodes.size), parts)
        within = np.arange(owner.size) - np.repeat(np.cumsum(parts) - parts, parts)
        width = ((radius_high - radius_low) / parts)[owner]
        radius = (
            (radius_low[owner] + within * width)[:, None] + _NODES * width[:, None]
        ).ravel()
        weight = (_WEIGHTS * width[:, None]).ravel()
        time = np.sqrt(radius * radius + delay * delay)
        if window == 0:
            w_p, w_s = c_p * radius, np.zeros_like(radius)
        else:
            w_p = c_p * np.sqrt(np.maximum(time * time - tau_p * tau_p, 0.0))
            w_s = c_s * radius
        pieces.append(
            (
                np.repeat(nodes[owner], _NODES.size),
                time,
                weight * radius / time,
                w_p,
                w_s,
                np.full(radius.size, window == 0),
            )
        )
    node, time, measure, w_p, w_s, between = (
        np.concatenate(part) for part in zip(*pieces)
    )
    # the sums of a kernel at the points of the pieces into the weights of
    # their nodes, the whole hats and their rising halves
    rising = measure * (time - node * step) / step
    falling = measure * ((node + 1) * step - time) / step
    point = np.arange(node.size)
    whole = sparse.csr_matrix(
        (
            np.concatenate((rising, falling)),
            (np.concatenate((node + 1, node)), np.concatenate((point, point))),
        ),
        shape=(count, node.size),
    )
    rise = sparse.csr_matrix((rising, (node + 1, point)), shape=(count, node.size))
    # F at every w the pieces take, from w = 0 on
    values = np.unique(np.concatenate((w_p, w_s, [0.0])))
    half = 0.5 * np.diff(values)
    abscissae, rule = np.polynomial.legendre.leggauss(3)
    u = (0.5 * (values[1:] + values[:-1]))[:, None] + half[:, None] * abscissae
    rule = half[:, None] * rule
    index_p = np.searchsorted(values, w_p)
    index_s = np.searchsorted(values, w_s)
    square = distance * distance
    scale = rule / (u * u + square) ** 1.5
    for part in range(0, k.size, chunk):
        modes = k[part : part + chunk][:, None, None]
        x = modes * u[None]
        bessel = u * special.j0(x)
        for kind in kinds:
            slope = bessel + modes * square * special.j1(x) if kind == "Q" else bessel
            primitive = np.zeros((modes.shape[0], values.size))
            primitive[:, 1:] = np.cumsum(np.einsum("kgn,gn->kg", slope, scale), axis=1)
            difference = primitive[:, index_p] - primitive[:, index_s]
            if kind == "Q":
                kernel = time * difference - np.where(between, time / distance, 0.0)
            else:
                kernel = -distance * time * difference
            full, halves = weights[kind]
            full[part : part + chunk] += (whole @ kernel.T).T
            halves[part : part + chunk] += (rise @ kernel.T).T
    return weights


class Inplane(Antiplane):
    """The spectral boundary integral method of an in-plane rupture on the
    plane x = 0, which slips along y (see the description of the module), in
    the arguments of Antiplane."""

    axis = "y"
    direction = np.array((1.0, 0.0))

    def __init__(self, medium, length, width=None, count=512, offset=0.0):
        super().__init__(medium, length, width, count, offset)
        self.lam = float(medium.stiffness[0, 1])
        self.fastest = math.sqrt((self.lam + 2.0 * self.mu) / medium.rho)
        self.ratio = self.fastest / self.speed
        self.kappa = 1.0 - 1.0 / (self.ratio * self.ratio)

    def _weights(self, step, steps):
        from scipy import special

        a = self.speed * self.k[:, None]
        rising, falling = _hat_halves(
            lambda t: _inplane_transfer_kernel(a * t[None, :], self.ratio),
            step,
            steps + 1,
        )
        self.direct, self.direct_end = rising + falling, rising
        self.images = []
        if self.width is None:
            return
        if not self.width / self.fastest > step:
            raise fault1d.Unsupported(
                "the step is not shorter than the time a wave takes across the period"
            )
        speeds = (self.speed, self.fastest)
        factor = (4.0 * self.k * self.k * self.speed**3)[:, None]
        n = 1
        while n * self.width / self.fastest < steps * step:
            distance = n * self.width
            delay = distance / self.speed
            ((rate_weights, rate_end),) = _delayed_weights(
                a[:, 0],
                delay,
                step,
                steps + 1,
                ((lambda a, x: -a[:, None, None] * special.j1(x), None),),
            )
            slip_weights, slip_end = _near_field_weights(
                self.k, speeds, distance, step, steps + 1, ("Q",)
            )["Q"]
            self.images.append(
                (
                    delay,
                    distance / self.fastest,
                    rate_weights,
                    factor * slip_weights,
                    rate_end,
                    factor * slip_end,
                )
            )
            n += 1

    def _memory(self, rate_k, slip_k, n, step):
        """As Antiplane._memory; the slip of an image transfers from the
        arrival of its P wave on, its slip rate from that of its S wave."""
        t = (n + 1) * step
        past = np.einsum("km,km->k", self.direct[:, 1 : n + 1], rate_k[:, n:0:-1])
        past = past + self.direct_end[:, n + 1] * rate_k[:, 0]
        history = np.zeros(self.k.size, dtype=complex)
        for delay, early, rate_weights, slip_weights, rate_end, slip_end in self.images:
            if early >= t:
                break
            convolution = (
                np.einsum("km,km->k", slip_weights[:, 1 : n + 1], slip_k[:, n:0:-1])
                + slip_end[:, n + 1] * slip_k[:, 0]
            )
            if delay < t:
                # V at t - delay, linear between the steps
                position = (t - delay) / step
                lower = min(int(math.floor(position)), n)
                fraction = position - lower
                upper = min(lower + 1, n)
                convolution = convolution + (1.0 - fraction) * rate_k[:, lower]
                convolution = convolution + fraction * rate_k[:, upper]
                convolution = convolution + (
                    np.einsum("km,km->k", rate_weights[:, 1 : n + 1], rate_k[:, n:0:-1])
                    + rate_end[:, n + 1] * rate_k[:, 0]
                )
            history = history - (self.mu / self.speed) * convolution
        return past, history

    def _transfer(self, slip_k, rate_k, past, history):
        damping = self.mu * self.kappa * self.k
        return -damping * (slip_k - (self.direct[:, 0] * rate_k + past)) + history

    def radiated(self, solution, x, y):
        """The particle velocity (x, y, z) and the stress (ordered as
        fault1d.VOIGT) that the slip of `solution` leaves at the point (x, y),
        per time of the solution, as the receivers record them: without the
        initial stress (see the description of the module)."""
        from scipy import special

        time = solution["time"]
        step = time[1] - time[0]
        count = time.size
        rate_k = np.fft.rfft(solution["slip-rate"], axis=1).T
        slip_k = np.fft.rfft(solution["slip"], axis=1).T
        k = self.k[:, None]
        c_s, c_p, mu, lam = self.speed, self.fastest, self.mu, self.lam
        a, b = c_s * self.k, c_p * self.k
        shape = rate_k.shape

        def zero():
            return [np.zeros(shape), np.zeros(shape)]

        # the weights of all the faults sum up, with the sign of what is odd in x
        shear_wave, pressure_wave, near_e, near_q, bessel, slope = (
            zero() for _ in range(6)
        )
        velocity_y = np.zeros(shape, dtype=complex)
        displacement_y = np.zeros(shape, dtype=complex)
        divergence = np.zeros(shape, dtype=complex)
        shear = np.zeros(shape, dtype=complex)

        def add(target, weights, sign=1.0):
            target[0] += sign * weights[0]
            target[1] += sign * weights[1]

        for distance, sign in self._sources(x, time[-1]):
            delay, early = distance / c_s, distance / c_p
            if not early > step:
                raise fault1d.Unsupported(
                    "the point lies closer to the fault than a wave travels in a step"
                )
            add(shear_wave, _propagator_weights(a, delay, step, count), sign)
            add(pressure_wave, _propagator_weights(b, early, step, count), sign)
            near_field = _near_field_weights(
                self.k, (c_s, c_p), distance, step, count, ("E", "Q")
            )
            add(near_e, near_field["E"], sign)
            add(near_q, near_field["Q"])
            kernels = (
                (
                    lambda a, x: special.j0(x),
                    lambda weight, radius, t: weight * radius / t,
                ),
                (lambda a, x: -a[:, None, None] * special.j1(x), None),
            )
            for target, weights in zip(
                (bessel, slope), _delayed_weights(a, delay, step, count, kernels)
            ):
                add(target, weights)
            arrived = _delayed(rate_k, delay, step)
            velocity_y += 0.5 * sign * arrived
            displacement_y += 0.5 * sign * _delayed(slip_k, delay, step)
            divergence += sign * _delayed(slip_k, early, step)
            shear += arrived
        near = k * k * c_s * c_s
        velocity_y += near * _convolve(*near_e, rate_k) - 0.5 * _convolve(
            *shear_wave, rate_k
        )
        displacement_y += near * _convolve(*near_e, slip_k) - 0.5 * _convolve(
            *shear_wave, slip_k
        )
        velocity_x = (
            1j
            * k
            * (
                c_s * c_s * _convolve(*near_q, rate_k)
                - 0.5 * c_s * _convolve(*bessel, rate_k)
            )
        )
        divergence = (1j * k * c_s * c_s / (c_p * c_p)) * (
            divergence - _convolve(*pressure_wave, slip_k)
        )
        shear = (
            -0.5
            * mu
            * (
                (shear + _convolve(*slope, rate_k)) / c_s
                + 4.0 * near * _convolve(*near_q, slip_k)
            )
        )

        # the real Fourier series at y, of the modes of numpy.fft.rfft
        factor = np.full(self.k.size, 2.0)
        factor[0] = 1.0
        if self.count % 2 == 0:
            factor[-1] = 1.0
        factor = factor * np.exp(1j * self.k * (y - self.y[0])) / self.count

        def at_y(values):
            return (factor @ values).real

        velocity = np.zeros((count, 3))
        velocity[:, 0] = at_y(velocity_x)
        velocity[:, 1] = at_y(velocity_y)
        normal_y = lam * divergence + 2.0 * mu * 1j * k * displacement_y
        stress = np.zeros((count, 6))
        stress[:, fault1d.VOIGT.index("xx")] = at_y(
            2.0 * (lam + mu) * divergence - normal_y
        )
        stress[:, fault1d.VOIGT.index("yy")] = at_y(normal_y)
        stress[:, fault1d.VOIGT.index("zz")] = at_y(lam * divergence)
        stress[:, fault1d.VOIGT.index("xy")] = at_y(shear)
        return velocity, stress


class _Loading:
    """The shear traction along `axis`, z for an antiplane and y for an
    in-plane rupture, and the normal stress on the welded fault, per point:
    the initial stress and the nucleations, as fault1d.Loading."""

    def __init__(self, fault, axis="z"):
        parameters = fault[None]
        size = np.shape(
            next(value for name, value in fault.items() if name is not None)
        )
        other = "y" if axis == "z" else "z"

        def field(name):
            values = fault.get(name)
            return np.zeros(size) if values is None else values

        tractions = ("T_n", "T_s", "T_d", "Tnuc_n", "Tnuc_s", "Tnuc_d")
        if any(name in fault for name in tractions):
            raise fault1d.Unsupported(
                "tractions in fault coordinates are not supported"
            )
        if np.any(field(f"s_x{other}") != 0.0):
            raise fault1d.Unsupported(f"the rupture takes no traction along {other}")
        self.normal = field("s_xx")
        self.shear = field(f"s_x{axis}")
        self.sources = []
        count = int(parameters.get("nucleationcount", 1))
        for index in range(count):
            suffix = "" if index == 0 else str(index + 1)
            prefix = f"nuc{suffix}_"
            if not any(f"{prefix}{c}" in fault for c in ("xx", "xy", "xz")):
                continue
            if np.any(field(f"{prefix}x{other}") != 0.0):
                raise fault1d.Unsupported(
                    f"the rupture takes no traction along {other}"
                )
            rise = fault.get(f"nuc{suffix}_rise_time")
            if rise is None:
                rise = float(parameters.get(f"t{suffix}_0", 0.0))
            onset = fault.get(f"nuc{suffix}_onset")
            if onset is None:
                onset = float(parameters.get(f"s{suffix}_0", 0.0))
            self.sources.append(
                (field(f"{prefix}xx"), field(f"{prefix}x{axis}"), onset, rise)
            )

    def __call__(self, t):
        normal, shear = self.normal, self.shear
        for dn, dz, onset, rise in self.sources:
            fraction = _smooth_step(t - onset, rise)
            normal = normal + fraction * dn
            shear = shear + fraction * dz
        return shear, normal


def _given(fault, *names, default=None):
    """A parameter of the fault: per point where the easi model gives it,
    otherwise the value of the parameter file, whose names are
    case-insensitive."""
    for name in names:
        if name in fault:
            return fault[name]
    for name in names:
        if name.lower() in fault[None]:
            return float(fault[None][name.lower()])
    if default is None:
        raise fault1d.Unsupported(f"the case gives none of {', '.join(names)}")
    return default


class _ThermalPressurization:
    """fault1d.ThermalPressurization per point: the temperature and the pore
    pressure at the centre of the shear zone, from the frictional heating
    tau V of each step, constant over it, against the kernels of fault1d,
    integrated exactly over the steps. The heating of a step is the mean of
    tau V at its ends, as the state of rate and state takes the mean slip
    rate. tp_halfWidthShearZone and tp_hydraulicDiffusivity may vary along
    the fault."""

    def __init__(self, fault):
        self.a_th = float(_given(fault, "tp_thermalDiffusivity", "alpha_th"))
        self.rho_c = float(_given(fault, "tp_heatCapacity", "rho_c"))
        self.lam = float(_given(fault, "tp_undrainedTPResponse", "tp_lambda"))
        self.t0 = float(_given(fault, "tp_iniTemp", "iniTemp"))
        # the initial pore pressure as SeisSol states it, -p_0
        self.pore_pressure = float(_given(fault, "tp_iniPressure", "iniPressure"))
        self.a_hy = np.asarray(_given(fault, "tp_hydraulicDiffusivity"), dtype=float)
        self.width = np.asarray(_given(fault, "tp_halfWidthShearZone"), dtype=float)
        if np.any(
            np.abs(self.a_hy - self.a_th) <= 1e-9 * np.maximum(self.a_hy, self.a_th)
        ):
            raise fault1d.Unsupported("equal hydraulic and thermal diffusivity")
        # along the fault, or one for all of it
        self.uniform = np.ptp(self.a_hy) == 0.0 and np.ptp(self.width) == 0.0
        if self.uniform:
            self.a_hy, self.width = self.a_hy.flat[0], self.width.flat[0]

    def _kernels(self, low, high):
        """fault1d.ThermalPressurization._kernels, per point."""
        width = self.width * self.width

        def root(diffusivity, lag):
            return np.sqrt(width + 2.0 * diffusivity * lag)

        scale = 1.0 / math.sqrt(2.0 * math.pi)
        thermal = (
            2.0 * scale * (high - low) / (root(self.a_th, high) + root(self.a_th, low))
        )
        pressure = (
            scale
            * (
                (root(self.a_hy, high) - root(self.a_hy, low))
                - (root(self.a_th, high) - root(self.a_th, low))
            )
            / (self.a_hy - self.a_th)
        )
        return thermal, pressure

    def prepare(self, step, count, points):
        """Tabulate the kernels of the end of a step for `count` steps of
        width `step`: the step j < n, seen from t_(n + 1), lies at the lags
        [(n - j) step, (n - j + 1) step], and the step n itself at [0, step]."""
        distance = np.arange(1, count + 1, dtype=float)
        if not self.uniform:
            distance = distance[:, None]
        # the most distant step first, in the order the heating is recorded
        thermal, pressure = self._kernels(distance * step, (distance + 1.0) * step)
        self.history = (thermal[::-1].copy(), pressure[::-1].copy())
        self.local = self._kernels(0.0, step)
        self.heating = np.zeros((count, points))
        self.recorded = 0

    def begin(self):
        """The part of the temperature and the pressure change at the end of
        the next step that the recorded steps leave."""
        count = self.recorded
        past = self.heating[:count]
        size = len(self.history[0])
        if self.uniform:
            self.past = tuple(kernel[size - count :] @ past for kernel in self.history)
        else:
            self.past = tuple(
                np.einsum("jn,jn->n", kernel[size - count :], past)
                for kernel in self.history
            )

    def changes(self, current):
        """Temperature and pressure change at the end of the next step, with
        the heating `current` within it."""
        temperature = self.past[0] + current * self.local[0]
        change = self.past[1] + current * self.local[1]
        return temperature / self.rho_c, self.lam * change / self.rho_c

    def record(self, heating):
        """Close a step with its heating."""
        self.heating[self.recorded] = heating
        self.recorded += 1


def imposed_slip(fault):
    """The y and z components of the slip that an imposed slip rate gives,
    per point, from (strike_slip, dip_slip) in the fault frame."""
    frame = fault1d.fault_frame(fault1d.Inputs(fault[None], {}))
    return frame.T @ np.array((fault["strike_slip"], fault["dip_slip"]))


def in_plane(fault):
    """Whether the fault slips along y, as an in-plane rupture, rather than
    along z: where its shear traction, or the slip it imposes, lies along y."""
    law = int(fault[None].get("fl", 16))
    if law in (33, 34):
        along_y, along_z = imposed_slip(fault)
        return bool(np.max(np.abs(along_y)) > np.max(np.abs(along_z)))
    names = [name for name in fault if isinstance(name, str)]
    along_y = [
        name for name in names if name == "s_xy" or re.fullmatch(r"nuc\d*_xy", name)
    ]
    return any(np.any(fault[name] != 0.0) for name in along_y)


def _friction_law(fault, axis="z"):
    """The friction law of the parameter file, per point; an imposed slip
    rate slips along `axis`."""
    law = int(fault[None].get("fl", 16))
    if law in (2, 16, 6, 1058):
        return _SlipWeakening(fault, 16 if law == 2 else law)
    if law in (3, 4, 7, 103):
        return _RateAndState(fault, law)
    if law in (33, 34):
        return _ImposedSlipRate(fault, law, axis)
    raise fault1d.Unsupported(
        "the spectral method takes slip weakening, rate and state and imposed "
        f"slip rates, not FL={law}"
    )


class _ImposedSlipRate:
    """The imposed slip rates FL=33 (regularized Yoffe) and FL=34 (Gaussian),
    as fault1d.ImposedSlipRate, per point: the slip (strike_slip, dip_slip)
    at the rate of the time since rupture_onset, whatever the traction, along
    `axis`. It has no state, and no friction coefficient."""

    def __init__(self, fault, law, axis="z"):
        self.law = law
        along = dict(zip("yz", imposed_slip(fault)))
        other = "y" if axis == "z" else "z"
        if np.any(np.abs(along[other]) > 1e-12 * np.max(np.abs(along[axis]))):
            raise fault1d.Unsupported(f"the rupture slips along {axis} only")
        self.slip = along[axis]
        self.onset = fault["rupture_onset"]
        if law == 33:
            self.tau_s = fault["tau_S"]
            self.tau_r = np.maximum(fault["tau_R"], self.tau_s)
        else:
            self.rise = fault["rupture_rise_time"]

    def initial_state(self, shear, normal):
        return np.zeros(np.shape(shear))

    def advance(self, state, rate, guess, step, t, normal):
        return state

    def rate(self, shear, normal, state, t, eta, start=None):
        local = t - self.onset
        if self.law == 33:
            return self.slip * _yoffe_rate(local, self.tau_s, self.tau_r)
        return self.slip * _smooth_step_rate(local, self.rise)

    def friction(self, state, rate, t, normal):
        return np.full(np.shape(rate), np.nan)

    def weakened(self, solution):
        return np.full(solution["slip-rate"].shape[1], np.nan)


class _SlipWeakening:
    """Linear slip weakening (FL=16, or FL=2), with or without a forced
    rupture time, its proxy of thermal pressurization (FL=1058) and its
    Prakash-Clifton regularization (FL=6), as fault1d.SlipWeakening, per
    point. Its state is the accumulated slip; with FL=6, also the strength
    that the slip rate sees, which relaxes towards the strength of slip
    weakening at the rate (|V| + v*) / L, and that strength itself.

    With lsw_healingThreshold V_h > 0, the fault heals instantaneously as
    SeisSol heals it: once its slip rate has exceeded V_h, a slip rate below
    V_h restores the static friction coefficient mu_s, which the strength of
    the next step takes, without undoing the accumulated slip; above V_h
    again, the friction coefficient is that of the accumulated slip. The
    state then holds as well the largest slip rate so far, whether the slip
    rate of its time heals the fault, and whether the strength of its time
    is the healed one."""

    def __init__(self, fault, law):
        self.law = law
        parameters = fault[None]
        threshold = float(parameters.get("lsw_healingthreshold", -1.0))
        # the slip rate is never below a threshold of zero or less
        self.healing = threshold if threshold > 0.0 else None
        self.mu_s, self.mu_d, self.d_c = fault["mu_s"], fault["mu_d"], fault["d_c"]
        self.cohesion = fault.get("cohesion", 0.0)
        self.forced = None
        if "forced_rupture_time" in fault:
            if law == 1058:
                raise fault1d.Unsupported(
                    "a forced rupture time with FL=1058 is not supported"
                )
            self.forced = (
                fault["forced_rupture_time"],
                float(parameters.get("t_0", 0.0)),
            )
        self.exponent = float(_given(fault, "TpProxyExponent", default=1.0 / 3.0))
        if law == 6:
            self.v_star = float(_given(fault, "pc_vstar", "v_star"))
            self.length = float(_given(fault, "pc_prakashLength", "L"))

    def _coefficient(self, t, accumulated, healed=None):
        """The friction coefficient of slip weakening, mu_s where `healed`."""
        if self.law == 1058:
            friction = self.mu_d + (self.mu_s - self.mu_d) / (
                1.0 + accumulated / self.d_c
            ) ** (self.exponent)
        else:
            weakening = np.minimum(accumulated / self.d_c, 1.0)
            if self.forced is not None:
                onset, rise = self.forced
                if rise > 0.0:
                    weakening = np.maximum(
                        weakening, np.clip((t - onset) / rise, 0.0, 1.0)
                    )
                else:
                    weakening = np.where(t >= onset, 1.0, weakening)
            friction = self.mu_s - (self.mu_s - self.mu_d) * weakening
        if healed is not None:
            friction = np.where(healed, self.mu_s, friction)
        return friction

    def strength(self, t, accumulated, normal, healed=None):
        """The strength of slip weakening, with mu_s where `healed`."""
        friction = self._coefficient(t, accumulated, healed)
        return -self.cohesion - np.minimum(normal, 0.0) * friction

    def friction(self, state, rate, t, normal):
        """The friction coefficient that SeisSol reports at the time of the
        state, which the strength of the next step takes: with healing, mu_s
        where the slip rate of the time heals the fault."""
        weakening = self._weakening(state)
        accumulated = weakening if self.law != 6 else weakening[0]
        healed = None if self.healing is None else state[-2] > 0.5
        return self._coefficient(t, accumulated, healed)

    def weakened(self, solution):
        """When the accumulated slip reaches d_c, per point."""
        return _first_crossings(
            solution["time"], solution["accumulated-slip"], self.d_c
        )

    def initial_state(self, shear, normal):
        accumulated = np.zeros(np.shape(shear))
        rows = [accumulated]
        if self.law == 6:
            strength = self.strength(0.0, accumulated, normal)
            rows += [strength, strength]
        if self.healing is not None:
            # no slip rate so far, and nothing healed
            rows += [accumulated] * 3
        return np.array(rows) if len(rows) > 1 else accumulated

    def _weakening(self, state):
        """The state of slip weakening, without that of the healing."""
        if self.healing is None:
            return state
        return state[0] if self.law != 6 else state[:3]

    def advance(self, state, rate, guess, step, t, normal):
        """The state at the end of a step, from that at its start, the slip
        rates at both ends, and the time and the normal stress at its end.
        With FL=6, the strength of slip weakening varies linearly over the
        step, and the regularized one follows it as it does then. Whether
        the strength of the end is healed, the slip rate of the start
        decides, as in SeisSol that of the sub-step before."""
        healed = None
        if self.healing is not None:
            healed = state[-2] > 0.5
        weakening = self._weakening(state)
        if self.law != 6:
            advanced = weakening + 0.5 * step * (np.abs(rate) + np.abs(guess))
        else:
            accumulated = weakening[0] + 0.5 * step * (np.abs(rate) + np.abs(guess))
            strength = self.strength(t, accumulated, normal, healed)
            mean = 0.5 * (np.abs(rate) + np.abs(guess))
            relaxation = (mean + self.v_star) * step / self.length
            weight = -np.expm1(-relaxation) / relaxation
            regularized = strength - (strength - weakening[2]) * weight
            regularized = regularized + (weakening[1] - weakening[2]) * np.exp(
                -relaxation
            )
            advanced = np.array([accumulated, regularized, strength])
        if self.healing is None:
            return advanced
        peak = state[-3]
        heals = (peak > self.healing) & (np.abs(guess) < self.healing)
        return np.vstack(
            [np.atleast_2d(advanced), np.maximum(peak, np.abs(guess)), heals, healed]
        )

    def rate(self, shear, normal, state, t, eta, start=None):
        """The slip rate that the traction of the welded fault leaves."""
        weakening = self._weakening(state)
        if self.law == 6:
            strength = weakening[1]
        else:
            healed = None if self.healing is None else state[-1] > 0.5
            strength = self.strength(t, weakening, normal, healed)
        return np.sign(shear) * np.maximum(np.abs(shear) - strength, 0.0) / eta


def _log_sinh(x):
    """log(sinh(x)) for x > 0, without overflow."""
    return x + np.log(-np.expm1(-2.0 * x) / 2.0)


def _arsinh_exp(log_z):
    """asinh(exp(log_z)), also where exp(log_z) overflows."""
    large = log_z > 20.0
    safe = np.where(large, 0.0, log_z)
    return np.where(
        large,
        log_z + math.log(2.0) + np.exp(-2.0 * np.where(large, log_z, 20.0)) / 4.0,
        np.arcsinh(np.exp(safe)),
    )


def _arsinh_exp_slope(log_z):
    """The derivative of _arsinh_exp."""
    return 1.0 / np.sqrt(1.0 + np.exp(-2.0 * np.clip(log_z, -350.0, 350.0)))


def _solve_slip_rate(friction, shear, sigma, eta, start=None):
    """The slip rate V >= 0 that solves shear - eta V = sigma f(V), per point
    of shear >= 0, by Newton in u = ln V within a bracket, from `start` where
    it lies inside; friction(u) gives f and df/du at V = exp(u). f grows with
    V, so that there is at most one root, and none where sigma f(0) holds the
    shear: the fault is locked there, and the slip rate 0."""
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):

        def residual(u):
            value, slope = friction(u)
            rate = np.exp(u)
            return sigma * value + eta * rate - shear, sigma * slope + eta * rate

        rate = np.zeros(np.shape(shear))
        moving = shear > 0.0
        moving &= residual(np.full(rate.shape, -np.inf))[0] < 0.0
        if not moving.any():
            return rate
        upper = np.where(moving, np.log(np.where(moving, shear, 1.0) / eta), 0.0)
        lower = upper - 1.0
        for _ in range(64):
            high = moving & (residual(lower)[0] > 0.0)
            if not high.any():
                break
            lower = np.where(high, lower - 2.0 * (upper - lower), lower)
        u = 0.5 * (lower + upper)
        if start is not None:
            guess = np.log(np.abs(start))
            u = np.where((lower < guess) & (guess < upper), guess, u)
        for _ in range(200):
            value, slope = residual(u)
            upper = np.where(value > 0.0, u, upper)
            lower = np.where(value < 0.0, u, lower)
            candidate = np.where(value == 0.0, u, u - value / slope)
            outside = (value != 0.0) & ~((lower < candidate) & (candidate < upper))
            candidate = np.where(outside, 0.5 * (lower + upper), candidate)
            change = np.abs(candidate - u)
            u = candidate
            if np.all(change[moving] <= 1e-15 * np.maximum(1.0, np.abs(u[moving]))):
                break
        else:
            raise fault1d.Unsupported("the slip rate does not settle")
        rate[moving] = np.exp(u[moving])
    return rate


class _RateAndState:
    """The rate-and-state laws FL=3 (aging), FL=4 (slip), FL=7 (severe
    velocity weakening) and FL=103 (fast velocity weakening), as
    fault1d.RateAndState, per point. Its state is psi, which a step advances
    in closed form at the mean slip rate V of the step: psi relaxes towards
    L / V with aging law, ln(V psi / L) towards 0 with slip law, and psi
    towards the steady state of V with FL=103, all at the rate V / L, and
    with FL=7, psi towards V L / V_0 at the rate V_0 / L."""

    def __init__(self, fault, law):
        self.law = law
        size = np.shape(fault["s_xx"])

        def field(*names):
            return np.broadcast_to(np.asarray(_given(fault, *names), dtype=float), size)

        # the names easi takes, of which the parameter file may give some
        self.a = field("rs_a")
        self.b = field("rs_b", "RS_b")
        self.length = field("rs_sl0", "RS_sl0")
        self.v0 = field("rs_sr0")
        self.f0 = field("rs_f0", "RS_f0")
        if law == 103:
            self.mu_w = field("rs_muw", "RS_muw", "mu_w")
            self.v_w = field("rs_srW")
        self.initial_rate = math.hypot(
            float(_given(fault, "rs_iniSlipRate1")),
            float(_given(fault, "rs_iniSlipRate2")),
        )

    def initial_state(self, shear, normal):
        """The state that holds the traction at the initial slip rate."""
        friction = np.abs(shear) / -normal
        rate = self.initial_rate
        if self.law == 7:
            weakening = (self.f0 + self.a * rate / (rate + self.v0) - friction) / self.b
            slipping = weakening < 1.0
            # where no state slips at the initial slip rate: locked, and healed
            healed = self.length * rate / self.v0
            safe = np.where(slipping, weakening, 0.0)
            return np.where(slipping, self.length * safe / (1.0 - safe), healed)
        psi = self.a * (
            math.log(2.0) + _log_sinh(friction / self.a) - np.log(rate / self.v0)
        )
        if self.law == 103:
            return psi
        return self.length / self.v0 * np.exp((psi - self.f0) / self.b)

    def advance(self, state, rate, guess, step, t=None, normal=None):
        """The state at the end of a step, from that at its start and the
        slip rates at both ends."""
        mean = 0.5 * (np.abs(rate) + np.abs(guess))
        if self.law == 7:
            # the steady state is large where V_0 / L is small, and the
            # relaxation towards it tiny: their product, and not their sum
            steady = mean * self.length / self.v0
            return state + (steady - state) * -np.expm1(-self.v0 / self.length * step)
        relaxation = mean * step / self.length
        with np.errstate(divide="ignore", invalid="ignore"):
            if self.law == 3:
                # psi e^-x + (L / V) (1 - e^-x), x = V t / L, and t where V = 0
                weight = np.where(
                    relaxation > 0.0,
                    -np.expm1(-relaxation)
                    / np.where(relaxation > 0.0, relaxation, 1.0),
                    1.0,
                )
                return state * np.exp(-relaxation) + step * weight
            if self.law == 4:
                logarithm = np.log(mean * state / self.length)
                return np.where(
                    mean > 0.0,
                    state * np.exp(np.expm1(-relaxation) * logarithm),
                    state,
                )
            steady = self._steady_state(np.where(mean > 0.0, mean, 1.0))
            return np.where(
                mean > 0.0, state - (state - steady) * -np.expm1(-relaxation), state
            )

    def friction(self, state, rate, t=None, normal=None):
        """The friction coefficient of the slip rate under the state."""
        speed = np.abs(rate)
        if self.law == 7:
            return (
                self.f0
                + self.a * speed / (speed + self.v0)
                - self.b * state / (state + self.length)
            )
        big_psi = state
        if self.law != 103:
            big_psi = self.f0 + self.b * np.log(self.v0 * state / self.length)
        with np.errstate(divide="ignore"):
            log_z = big_psi / self.a - np.log(2.0 * self.v0) + np.log(speed)
        return self.a * _arsinh_exp(log_z)

    def weakened(self, solution):
        """As fault1d.RateAndState.weakened, per point: when the friction
        coefficient of FL=103 falls to mu_w + 0.05 (f_0 - mu_w), not before
        the rupture time; NaN for the other laws."""
        time = solution["time"]
        weakened = np.full(solution["slip-rate"].shape[1], np.nan)
        if self.law != 103:
            return weakened
        threshold = self.mu_w + 0.05 * (self.f0 - self.mu_w)
        for point, start in enumerate(solution["rupture-time"]):
            if np.isnan(start):
                continue
            after = time >= start
            crossing = fault1d._first_crossing(
                time[after], -solution["friction"][after, point], -threshold[point]
            )
            if crossing is not None:
                weakened[point] = crossing
        return weakened

    def _steady_state(self, rate):
        """The steady state of FL=103 at the slip rate."""
        ratio = rate / self.v0
        f_lv = self.f0 - (self.b - self.a) * np.log(ratio)
        mu_ss = self.mu_w + (f_lv - self.mu_w) / (1.0 + (rate / self.v_w) ** 8) ** 0.125
        return self.a * (math.log(2.0) - np.log(ratio) + _log_sinh(mu_ss / self.a))

    def rate(self, shear, normal, state, t, eta, start=None):
        """The slip rate that the traction of the welded fault leaves, under
        the state."""
        sigma = -normal
        if np.any(sigma <= 0.0):
            raise fault1d.Unsupported("the normal stress is not compressive")
        if self.law == 7:
            weakening = self.f0 - self.b * state / (state + self.length)

            def friction(u):
                rate = np.exp(u)
                value = weakening + self.a * rate / (rate + self.v0)
                return value, self.a * rate * self.v0 / (rate + self.v0) ** 2

        else:
            big_psi = state
            if self.law != 103:
                big_psi = self.f0 + self.b * np.log(self.v0 * state / self.length)
            offset = big_psi / self.a - np.log(2.0 * self.v0)

            def friction(u):
                log_z = offset + u
                return self.a * _arsinh_exp(log_z), self.a * _arsinh_exp_slope(log_z)

        magnitude = _solve_slip_rate(friction, np.abs(shear), sigma, eta, start)
        return np.sign(shear) * magnitude


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------


def _klein_gordon(medium, k, half_width, slip, end, step, probe):
    """The traction mu dU/dx (0, t) of the mode U(x, t) cos(k y) of an
    antiplane slip slip(t) cos(k y) on the faults x = n W of a solid periodic
    across them, and U, dU/dt and dU/dx at the grid point nearest to
    x = probe, from U_tt = c^2 (U_xx - k^2 U) on 0 < x < W / 2, with
    U = slip / 2 on the fault and U = 0 halfway to the next, by finite
    differences of the second order at `step`, and at the Courant number 1,
    or just below, at which they carry the waves along x without error.
    Returns the traction, the position of the grid point and its U, dU/dt
    and dU/dx, per step."""
    mu = float(medium.stiffness[3, 3])
    speed = math.sqrt(mu / medium.rho)
    points = int(math.floor(half_width / (speed * step)))
    spacing = half_width / points
    courant = (speed * step / spacing) ** 2
    index = int(round(probe / spacing))
    previous = np.zeros(points + 1)
    current = np.zeros(points + 1)
    traction = [0.0]
    field = [(0.0, 0.0, 0.0)]
    for n in range(1, int(round(end / step)) + 2):
        following = np.empty_like(current)
        following[1:-1] = (
            2.0 * current[1:-1]
            - previous[1:-1]
            + courant * (current[2:] - 2.0 * current[1:-1] + current[:-2])
            - (speed * step * k) ** 2 * current[1:-1]
        )
        following[0] = 0.5 * slip(n * step)
        following[-1] = 0.0
        # the rate of the previous step, centred
        rate = (following[index] - previous[index]) / (2.0 * step)
        field[-1] = (field[-1][0], rate, field[-1][2])
        previous, current = current, following
        traction.append(
            mu * (-3.0 * current[0] + 4.0 * current[1] - current[2]) / (2.0 * spacing)
        )
        gradient = (current[index + 1] - current[index - 1]) / (2.0 * spacing)
        field.append((current[index], 0.0, gradient))
    return np.array(traction[:-1]), index * spacing, np.array(field[:-1])


def _plane_strain(medium, k, half_width, slip, end, step, probe, spacing):
    """The in-plane mode u_y = U(x, t) cos(k y), u_x = X(x, t) sin(k y) of
    a slip slip(t) cos(k y) on the faults x = n W of a solid periodic across
    them, from

        rho U_tt = mu U_xx - (lambda + 2 mu) k^2 U + (lambda + mu) k X_x,
        rho X_tt = (lambda + 2 mu) X_xx - mu k^2 X - (lambda + mu) k U_x,

    on 0 < x < W / 2, with U = slip / 2 and sigma_xx = 0, i.e.
    (lambda + 2 mu) X_x = lambda k U, on the fault, and U = 0 and X_x = 0
    halfway to the next, by finite differences of the second order on the
    grid of `spacing`, at the Courant number 1/2 of the P wave. Returns, per
    step, the traction mu (k X + U_x)(0, t) and, at the grid point nearest
    to x = probe, U, X, U_t, X_t, U_x and X_x; and the position of the
    point."""
    mu = float(medium.stiffness[3, 3])
    lam = float(medium.stiffness[0, 1])
    modulus = lam + 2.0 * mu
    points = int(round(half_width / spacing))
    h = half_width / points
    index = int(round(probe / h))
    substeps = int(math.ceil(step / (0.5 * h / math.sqrt(modulus / medium.rho))))
    dt = step / substeps
    along = [np.zeros(points + 1) for _ in range(2)]
    across = [np.zeros(points + 1) for _ in range(2)]
    traction = [0.0]
    field = [np.zeros(6)]
    total = int(round(end / step))
    # one substep beyond the end, for the rates of the last sample
    for n in range(1, total * substeps + 2):
        previous_u, u = along
        previous_x, x = across
        # sigma_xx = 0 on the fault; X_x = 0, and U odd, halfway to the next
        x_out = np.concatenate(
            ([x[1] - 2.0 * h * lam * k * u[0] / modulus], x, [x[-2]])
        )
        u_out = np.concatenate(([0.0], u, [-u[-2]]))
        slope_u = (u_out[2:] - u_out[:-2]) / (2.0 * h)
        slope_u[0] = (-3.0 * u[0] + 4.0 * u[1] - u[2]) / (2.0 * h)
        slope_x = (x_out[2:] - x_out[:-2]) / (2.0 * h)
        curvature_u = (u_out[2:] - 2.0 * u + u_out[:-2]) / (h * h)
        curvature_x = (x_out[2:] - 2.0 * x + x_out[:-2]) / (h * h)
        next_u = (
            2.0 * u
            - previous_u
            + dt
            * dt
            / medium.rho
            * (mu * curvature_u - modulus * k * k * u + (lam + mu) * k * slope_x)
        )
        next_x = (
            2.0 * x
            - previous_x
            + dt
            * dt
            / medium.rho
            * (modulus * curvature_x - mu * k * k * x - (lam + mu) * k * slope_u)
        )
        next_u[0] = 0.5 * slip(n * dt)
        next_u[-1] = 0.0
        if (n - 1) % substeps == 0:
            # the rates of the sample of the step before, centred
            field[-1][2] = (next_u[index] - previous_u[index]) / (2.0 * dt)
            field[-1][3] = (next_x[index] - previous_x[index]) / (2.0 * dt)
        along, across = [u, next_u], [x, next_x]
        if n % substeps == 0:
            slope = (-3.0 * next_u[0] + 4.0 * next_u[1] - next_u[2]) / (2.0 * h)
            traction.append(mu * (k * next_x[0] + slope))
            field.append(
                np.array(
                    (
                        next_u[index],
                        next_x[index],
                        0.0,
                        0.0,
                        (next_u[index + 1] - next_u[index - 1]) / (2.0 * h),
                        (next_x[index + 1] - next_x[index - 1]) / (2.0 * h),
                    )
                )
            )
    return np.array(traction[: total + 1]), np.array(field[: total + 1]), index * h


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
    # the plane waves of that slip at x = 0.3 W: the waves of the fault, at
    # x + n W, carry +V / 2, those of the next, at (n + 1) W - x, -V / 2, and
    # all of them the stress -eta V
    x = 0.3 * width
    velocity, stress = solver.radiated(series, x, 0.0)
    time = series["time"]
    exact_velocity = np.zeros(time.size)
    exact_stress = np.zeros(time.size)
    for n in range(4):
        for distance, sign in ((x + n * width, 1.0), ((n + 1) * width - x, -1.0)):
            arrived = np.array([rate(t - distance / solver.speed) for t in time])
            exact_velocity += 0.5 * sign * arrived
            exact_stress -= solver.eta * arrived
    check(
        "the plane waves of a uniform slip across periodic faults",
        max(
            float(np.max(np.abs(velocity[:, 2] - exact_velocity)))
            / float(np.max(np.abs(exact_velocity))),
            float(np.max(np.abs(stress[:, 4] - exact_stress)))
            / float(np.max(np.abs(exact_stress))),
            float(np.max(np.abs(np.delete(velocity, 2, axis=1)))),
            float(np.max(np.abs(np.delete(stress, (3, 4), axis=1)))),
        ),
        3e-3,
    )

    # the same slip, weakened to 0.6 by a forced rupture at t = 0 over t_0
    # from mu_s = 0.7, to which the fault heals once the echoes slow it below
    # lsw_healingThreshold = 0.5 m/s, and under which it stays locked: the
    # traction stays below 0.7 * 120 MPa. The slip rate of a step heals the
    # strength of the next, as in SeisSol, so that the fault stops up to two
    # steps after the time the slip rate falls below the threshold.
    healing = {
        **loading,
        "mu_s": 0.7,
        "mu_d": 0.6,
        "d_c": 1.0,
        "forced_rupture_time": 0.0,
    }
    fault = {name: np.full(8, value) for name, value in healing.items()}
    fault[None] = {"t_0": 0.02, "lsw_healingthreshold": 0.5}
    solver = Antiplane(medium, 800.0, width, 8)
    series = solver.solve(fault, end, step)

    def free(t):
        if t <= 0.0:
            return 0.0
        shear = 70.0e6 + 10.0e6 * fault1d.smooth_step(t, 0.02)
        strength = (0.7 - 0.1 * min(t / 0.02, 1.0)) * 120.0e6
        echoes = sum(free(t - n * delay) for n in range(1, int(t / delay) + 1))
        return max(0.0, (shear - strength) / solver.eta - 2.0 * echoes)

    time = series["time"]
    unhealed = np.array([free(t) for t in time])
    # the first time it falls below the threshold, between the steps
    below = int(np.argmax((np.maximum.accumulate(unhealed) > 0.5) & (unhealed < 0.5)))
    early, late = time[below - 1], time[below]
    for _ in range(40):
        middle = 0.5 * (early + late)
        early, late = (middle, late) if free(middle) >= 0.5 else (early, middle)
    healed = np.where(time < late, unhealed, 0.0)
    exact = np.concatenate(([0.0], np.cumsum(0.5 * step * (healed[1:] + healed[:-1]))))
    check(
        "the slip of a uniform slip that heals across periodic faults",
        float(np.max(np.abs(series["slip"][:, 0] - exact)) / exact[-1]),
        1e-3,
    )
    stop = int(np.argmax(time >= late)) + np.argmax(
        series["slip-rate"][time >= late, 0] <= 0.0
    )
    check(
        "the steps from the healing to the stop of the slip",
        (time[stop] - late) / step,
        2.0,
    )
    check(
        "the slip rate and the friction coefficient but mu_s after the stop",
        float(
            np.max(np.abs(series["slip-rate"][stop:, 0]))
            + np.max(np.abs(series["friction"][stop:, 0] - 0.7))
        ),
        0.0,
    )

    # a mode of a prescribed slip, against the finite differences of its wave
    # equation across the faults: the transfer kernel and the images, and the
    # propagators of the field off the fault
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
    exact, probe, field = _klein_gordon(
        medium, solver.k[mode], 0.5 * width, slip, end, step, 0.15 * width
    )
    check(
        f"the traction of the mode k = {solver.k[mode]:.4g} / m across periodic faults",
        float(np.max(np.abs(computed - exact)) / np.max(np.abs(exact))),
        1e-4,
    )
    # the same mode at the points, cos(k (y - y_0)) of the series of rfft, and
    # its field off the fault, which is odd in x, but for mu du/dx; once
    # across the periodic faults and once in a solid that extends to
    # infinity, where the finite differences reach far enough that what their
    # boundary reflects does not return to the point within the time
    angle = solver.k[mode] * (solver.y - solver.y[0])
    prescribed = {
        "time": time,
        "slip-rate": np.outer(rate_k[mode].real, np.cos(angle)),
        "slip": np.outer(slip_k[mode].real, np.cos(angle)),
    }
    infinite = Antiplane(medium, 800.0, None, 16)
    _, far, alone = _klein_gordon(
        medium, solver.k[mode], 2.0 * solver.speed * end, slip, end, step, probe
    )
    for method, x, values, side, y in (
        (solver, probe, field, 1.0, 30.0),
        (solver, probe, field, -1.0, -110.0),
        (infinite, far, alone, -1.0, 30.0),
    ):
        velocity, stress = method.radiated(prescribed, side * x, y)
        phase = solver.k[mode] * (y - solver.y[0])
        expected = (
            (velocity[:, 2], side * values[:, 1] * math.cos(phase)),
            (stress[:, 4], solver.mu * values[:, 2] * math.cos(phase)),
            (
                stress[:, 3],
                -side * solver.mu * solver.k[mode] * values[:, 0] * math.sin(phase),
            ),
        )
        solid = "across periodic faults" if method.width else "alone"
        check(
            f"the field of the mode k = {solver.k[mode]:.4g} / m at "
            f"({side * x:.4g}, {y:.4g}) {solid}",
            max(
                float(np.max(np.abs(computed - exact)) / np.max(np.abs(exact)))
                for computed, exact in expected
            ),
            1e-4,
        )

    # the in-plane slip of the same mode and of a mode of eight times the
    # wave number, whose P waves reach the next fault before the S waves, and
    # its field off the fault, against the finite differences of the plane
    # strain of the mode on a grid of 2.5 cm, to whose error the deviations
    # converge
    lam, mu = float(medium.stiffness[0, 1]), float(medium.stiffness[3, 3])
    short = 0.06
    inplane_time = time[time <= short + 0.5 * step]
    for mode, probe_y in ((5, 30.0), (40, 33.0)):
        inplane = Inplane(medium, 800.0, width, 2 * mode + 2)
        k = inplane.k[mode]
        rate_k = np.zeros((inplane.k.size, inplane_time.size), dtype=complex)
        slip_k = np.zeros_like(rate_k)
        rate_k[mode] = [slip_rate(t) for t in inplane_time]
        slip_k[mode] = [slip(t) for t in inplane_time]
        computed = inplane.transfer(rate_k, slip_k, step)[mode].real
        computed = computed - inplane.eta * rate_k[mode].real
        exact, values, probe = _plane_strain(
            medium, k, 0.5 * width, slip, short, step, 0.15 * width, 0.025
        )
        check(
            f"the traction of the in-plane mode k = {k:.4g} / m across periodic faults",
            float(np.max(np.abs(computed - exact)) / np.max(np.abs(exact))),
            2e-4,
        )
        angle = k * (inplane.y - inplane.y[0])
        prescribed = {
            "time": inplane_time,
            "slip-rate": np.outer(rate_k[mode].real, np.cos(angle)),
            "slip": np.outer(slip_k[mode].real, np.cos(angle)),
        }
        u, x, rate_u, rate_x, slope_u, slope_x = values.T
        for side, y in ((1.0, probe_y), (-1.0, probe_y - 140.0)):
            velocity, stress = inplane.radiated(prescribed, side * probe, y)
            phase = k * (y - inplane.y[0])
            cosine, sine = math.cos(phase), math.sin(phase)
            # u_y is odd in x and u_x even
            expected = (
                (velocity[:, 0], rate_x * sine),
                (velocity[:, 1], side * rate_u * cosine),
                (stress[:, 5], mu * (k * x + slope_u) * cosine),
                (
                    stress[:, 0],
                    side * ((lam + 2.0 * mu) * slope_x - lam * k * u) * sine,
                ),
                (
                    stress[:, 1],
                    side * (lam * slope_x - (lam + 2.0 * mu) * k * u) * sine,
                ),
                (stress[:, 2], side * lam * (slope_x - k * u) * sine),
            )
            check(
                f"the in-plane field of the mode k = {k:.4g} / m at "
                f"({side * probe:.4g}, {y:.4g}) across periodic faults",
                max(
                    float(np.max(np.abs(computed - exact)) / np.max(np.abs(exact)))
                    for computed, exact in expected
                ),
                1e-2 if mode == 40 else 3e-4,
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

    # a uniform slip of every other law on a fault between two solids that
    # extend to infinity: the 0D problem of fault1d, which integrates its own
    # way, with what SeisSol reports from the Gauss points of the fault
    def uniform(law, rupture, given, method=Antiplane):
        problem = fault1d.Problem(rupture, given, medium, medium)
        exact = problem.solve(0.06, 1e-5)
        solver = method(medium, 800.0, None, 8)
        fault = {name: np.full(8, value) for name, value in given.items()}
        fault[None] = {key.lower(): value for key, value in rupture.items()}
        computed = solver.solve(fault, 0.06, 3e-5)
        deviation = {}
        for name in ("slip-rate", "friction", "accumulated-slip", "peak-slip-rate"):
            if name not in exact:
                continue
            reference = exact[name][::3]
            error = np.max(np.abs(computed[name][:, 0] - reference))
            deviation[name] = float(error / np.max(np.abs(reference)))
        for name in ("pore-pressure", "temperature"):
            if name in computed:
                # the change, which is small against the initial value
                values = computed[name][:, 0] - computed[name][0, 0]
                reference = exact[name][::3] - exact[name][0]
                error = np.max(np.abs(values - reference))
                deviation[name] = float(error / np.max(np.abs(reference)))
        if "state" in computed:
            values, reference = computed["state"][:, 0], exact["state"][::3]
            if law in (3, 4):
                # psi of aging and slip law spans orders of magnitude
                error = np.max(np.abs(np.log(values) - np.log(reference)))
            else:
                error = np.max(np.abs(values - reference)) / np.max(np.abs(reference))
            deviation["state"] = float(error)
        label = f"FL={law}"
        if int(rupture.get("thermalPress", 0)):
            label += " with thermal pressurization"
        if method is Inplane:
            label += " in plane"
        check(f"a uniform slip under {label}", deviation.pop("slip-rate"), 1e-5)
        check(
            f"what the Gauss points report of a uniform slip under {label}",
            max(deviation.values()),
            3e-5,
        )
        # the times of the events, which both take between the steps of the
        # method, where a slip rate that sets in with a kink, as that of FL=6
        # and FL=7 does, crosses fault1d.RUPTURE_SLIP_RATE elsewhere than
        # between the finer steps of fault1d
        sampled = {name: value[::3] for name, value in exact.items() if np.ndim(value)}
        sampled["rupture-time"] = fault1d._first_crossing(
            sampled["time"], sampled["slip-rate"], fault1d.RUPTURE_SLIP_RATE
        )
        events = ["rupture-time"]
        if not problem.imposed:
            sampled["dynamic-stress-time"] = problem.friction.weakened(sampled)
            events.append("dynamic-stress-time")
        deviation = 0.0
        for name in events:
            value, reference = computed[name][0], sampled[name]
            if reference is None:
                deviation = max(deviation, 0.0 if np.isnan(value) else math.inf)
            else:
                deviation = max(deviation, abs(value - reference) / 3e-5)
        check(
            f"the rupture and dynamic stress times of a uniform slip under {label}, "
            "in steps",
            deviation,
            1e-2,
        )

    # the imposed slip rates of fault1d-imposed.yaml, along z
    imposed = {"strike_slip": 0.0, "dip_slip": -1.0, "rupture_onset": 1.5625e-2}
    imposed.update({"tau_S": 1.5625e-2, "tau_R": 6.25e-2, "rupture_rise_time": 6.25e-2})
    for law in (33, 34):
        uniform(law, {"FL": law, "XRef": -1.0e5}, imposed)
    # in plane, along y, where the uniform slip transfers nothing either
    uniform(
        33,
        {"FL": 33, "XRef": -1.0e5},
        {**imposed, "strike_slip": -1.0, "dip_slip": 0.0},
        Inplane,
    )
    given = {"s_xx": -120.0e6, "s_xy": 72.0e6, "nuc_xy": 14.0e6}
    given.update({"mu_s": 0.677, "mu_d": 0.525, "d_c": 0.08})
    uniform(16, {"FL": 16, "t_0": 0.02, "XRef": -1.0e5}, given, Inplane)
    for law, extra in (
        (6, {"pc_vstar": 0.1, "pc_prakashLength": 0.01}),
        (1058, {}),
    ):
        given = {"s_xx": -120.0e6, "s_xz": 72.0e6, "nuc_xz": 14.0e6}
        given.update({"mu_s": 0.677, "mu_d": 0.525, "d_c": 0.08})
        if law == 1058:
            given.update({"d_c": 0.01, "s_xz": 75.0e6})
        rupture = {"FL": law, "t_0": 0.02, "XRef": -1.0e5, **extra}
        uniform(law, rupture, given)
    parameters = {
        "rs_f0": 0.6,
        "rs_sr0": 1e-6,
        "rs_b": 0.012,
        "rs_iniSlipRate1": 1e-16,
        "rs_iniSlipRate2": 0.0,
        "t_0": 0.02,
        "XRef": -1.0e5,
    }
    loading = {"s_xx": -120.0e6, "s_xz": 75.0e6, "nuc_xz": 40.0e6}
    # thermal pressurization, of fault1d-rs-fvw-tp.yaml and SEISSOL_TP, with an
    # initial pore pressure
    thermal = {
        "thermalPress": 1,
        "tp_iniTemp": 483.15,
        "tp_iniPressure": -20.0e6,
        "tp_thermalDiffusivity": 1.0e-6,
        "tp_heatCapacity": 2.7e6,
        "tp_undrainedTPResponse": 0.1e6,
    }
    shear_zone = {"tp_halfWidthShearZone": 0.005, "tp_hydraulicDiffusivity": 4.0e-4}
    for law, extra, changes, pressurized in (
        (3, {"rs_a": 0.008, "rs_sl0": 0.02}, {}, False),
        (4, {"rs_a": 0.008, "rs_sl0": 0.1}, {}, False),
        (7, {"rs_a": 0.008, "rs_sl0": 0.02}, {"s_xz": 70.8e6}, False),
        (103, {"rs_a": 0.01, "rs_sl0": 0.2, "rs_srW": 0.1}, {"s_xz": 40.0e6}, False),
        (103, {"rs_a": 0.01, "rs_sl0": 0.2, "rs_srW": 0.1}, {"s_xz": 34.0e6}, True),
    ):
        given = {**loading, **changes, **extra}
        rupture = {**parameters, "FL": law}
        if law == 103:
            given["nuc_xz"] = 45.0e6
            rupture.update({"rs_b": 0.014, "rs_muW": 0.2})
        if pressurized:
            given.update(shear_zone)
            rupture.update(thermal)
        uniform(law, rupture, given)
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
