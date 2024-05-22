import numpy as np


def asymmetric_cosine(t, trise, tfall=None):
    """
    Initialize a source time function with asymmetric cosine, normalized to 1
    modified from instaseis
    https://github.com/krischer/instaseis/blob/master/instaseis/source.py#L174
    """
    # initialize
    if not isinstance(tfall, np.ndarray):
        tfall = trise
    # build slices
    slrise = np.logical_and(t > 0.0, t <= trise)
    slfall = np.logical_and(t > trise, t <= trise + tfall)

    # compute stf
    asc = np.zeros_like(t)
    asc[slrise] = 1.0 - np.cos(np.pi * t[slrise] / trise[slrise])
    asc[slfall] = 1.0 - np.cos(
        np.pi * (t[slfall] - trise[slfall] + tfall[slfall]) / tfall[slfall]
    )

    # normalize
    asc /= trise + tfall
    return asc
