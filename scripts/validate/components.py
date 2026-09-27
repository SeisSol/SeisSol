# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The vectors and tensors whose components the ``compare-*.py`` scripts compare together.

A component is measured relative to the largest reference norm among the components of
its vector or tensor, not relative to its own: a component that is zero or nearly so,
like the shear stress in water, the tangential displacement at the free surface of a
fluid or the normal traction change on a fault that the setup keeps symmetric, holds
rounding noise only, and relative to that noise any other rounding would look like a
change of order one.
"""

import re

# output quantities that are components of one vector or tensor
COMPONENTS = {
    **{name: "stress" for name in ("s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz")},
    **{
        name: "strain rate"
        for name in ("epsxx", "epsyy", "epszz", "epsxy", "epsyz", "epsxz")
    },
    **{name: "velocity" for name in ("v1", "v2", "v3")},
    **{name: "displacement" for name in ("u1", "u2", "u3")},
    **{name: "fluid velocity" for name in ("v1_f", "v2_f", "v3_f")},
    **{name: "traction" for name in ("T_s", "T_d", "P_n")},
    **{name: "initial traction" for name in ("Ts0", "Td0", "Pn0")},
    **{name: "slip" for name in ("Sls", "Sld")},
    **{name: "slip rate" for name in ("SRs", "SRd")},
}


def component_group(name: str) -> str:
    """The vector or tensor a quantity is a component of, per simulation.

    A fused simulation appends its index to every quantity: directly in a receiver
    (s_xx3, v13), after a dash on the fault (SRs-4) and in the mesh outputs (v1-2); the
    components of one simulation form one group. A quantity that is no component is a
    group of its own.
    """
    for end in range(len(name), 0, -1):
        base, suffix = name[:end], name[end:]
        if base in COMPONENTS and re.fullmatch(r"(-?\d+)?", suffix):
            return COMPONENTS[base] + suffix
    return name


def group_scales(norms: dict) -> dict:
    """The largest norm of each group, from the reference norm of every quantity."""
    scales = {}
    for name, norm in norms.items():
        group = component_group(name)
        scales[group] = max(scales.get(group, 0.0), norm)
    return scales
