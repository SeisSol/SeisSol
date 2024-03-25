from math import exp
import numpy as np


def SmoothStep(time, Tnuc):
    positive_ids = time > 0
    result = np.zeros_like(time)
    result[positive_ids] = np.exp(
        (time[positive_ids] - Tnuc[positive_ids]) ** 2
        / (time[positive_ids] * (time[positive_ids] - 2.0 * Tnuc[positive_ids]))
    )
    result[time >= Tnuc] = 1.0
    return result


def GaussianSTF(time, Tnuc, dt):
    Gnuc = np.where(
        (time > 0) & (time < Tnuc),
        SmoothStep(time, Tnuc) - SmoothStep(time - dt, Tnuc),
        0,
    )
    return Gnuc
