from math import exp


def SmoothStep(time, Tnuc):
    if time <= 0:
        return 0
    elif time < Tnuc:
        return exp((time - Tnuc) ** 2 / (time * (time - 2.0 * Tnuc)))
    else:
        return 1.0


def GaussianSTF(time, Tnuc, dt):
    if time > 0 and time < Tnuc:
        Gnuc = SmoothStep(time, Tnuc) - SmoothStep(time - dt, Tnuc)
    else:
        Gnuc = 0
    return Gnuc
