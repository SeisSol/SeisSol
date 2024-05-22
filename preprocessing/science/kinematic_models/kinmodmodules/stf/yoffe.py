from math import sqrt, asin, atan, pi


def C1(t, ts, tr):
    """ C1 to C6 are analytical functions
    used for building the regularized Yoffe function
    """
    return (
        (0.5 * t + 0.25 * tr) * sqrt(t * (tr - t))
        + (t * tr - tr * tr) * asin(sqrt(t / tr))
        - 0.75 * tr * tr * atan(sqrt((tr - t) / t))
    )


def C2(t, ts, tr):
    return 0.375 * pi * tr * tr


def C3(t, ts, tr):
    return (
        (ts - t - 0.5 * tr) * sqrt((t - ts) * (tr - t + ts))
        + tr * (2 * tr - 2 * t + 2 * ts) * asin(sqrt((t - ts) / tr))
        + 1.5 * tr * tr * atan(sqrt((tr - t + ts) / (t - ts)))
    )


def C4(t, ts, tr):
    # fixed 2 typos in the second term
    return (
        (-ts + 0.5 * t + 0.25 * tr) * sqrt((t - 2.0 * ts) * (tr - t + 2.0 * ts))
        - tr * (tr - t + 2.0 * ts) * asin(sqrt((t - 2.0 * ts) / tr))
        - 0.75 * tr * tr * atan(sqrt((tr - t + 2.0 * ts) / (t - 2.0 * ts)))
    )


def C5(t, ts, tr):
    return 0.5 * pi * tr * (t - tr)


def C6(t, ts, tr):
    return 0.5 * pi * tr * (2.0 * ts - t + tr)


def regularizedYoffe(t, ts, tr):
    """
    Implementation of the regularized Yoffe function
    defined in Appendix of Tinti et al. (2005), https://doi.org/10.1785/0120040177
    """
    assert ts <= tr, f"ts:{ts} > tr:{tr}"
    K = 2.0 / (pi * tr * ts * ts)
    if tr > 2.0 * ts:
        if t <= 0:
            return 0
        elif t <= ts:
            return K * (C1(t, ts, tr) + C2(t, ts, tr))
        elif t <= 2.0 * ts:
            return K * (C1(t, ts, tr) - C2(t, ts, tr) + C3(t, ts, tr))
        elif t < tr:
            return K * (C1(t, ts, tr) + C3(t, ts, tr) + C4(t, ts, tr))
        elif t < tr + ts:
            return K * (C3(t, ts, tr) + C4(t, ts, tr) + C5(t, ts, tr))
        elif t < tr + 2.0 * ts:
            return K * (C4(t, ts, tr) + C6(t, ts, tr))
        else:
            return 0
    else:
        if t <= 0:
            return 0
        elif t <= ts:
            return K * (C1(t, ts, tr) + C2(t, ts, tr))
        elif t < tr:
            return K * (C1(t, ts, tr) - C2(t, ts, tr) + C3(t, ts, tr))
        elif t <= 2.0 * ts:
            return K * (C5(t, ts, tr) + C3(t, ts, tr) - C2(t, ts, tr))
        elif t < tr + ts:
            return K * (C3(t, ts, tr) + C4(t, ts, tr) + C5(t, ts, tr))
        elif t < tr + 2.0 * ts:
            return K * (C4(t, ts, tr) + C6(t, ts, tr))
        else:
            return 0
