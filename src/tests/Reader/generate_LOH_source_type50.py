#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2023-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff


import numpy as np


def main():
    fout = open("fsrm_source1.dat", "w")

    # Note that we use a moment tensor definition that differs with the standards in Seismology (e.g. Komatitsch and Tromp, 1999).
    # the moment tensor is opposite to the usual definition

    # Note that strike, dip and rake allow rotating the seismic moment tensor.
    # Typically, the moment tensor does not need to be rotated and strike = dip = rake = 0.

    dt = 0.02
    T = 0.1
    vtime = np.arange(0, 4, dt)
    momentRate = 1 / T**2 * vtime * np.exp(-vtime / T)

    header = f""" seismic moment tensor
    0.00000E+00 -0.10000E+19 0.00000E+00
    -0.10000E+19 0.00000E+00 0.00000E+00
    0.00000E+00 0.00000E+00 0.00000E+00
    number of point sources
                        1
    x,y,z,strike, dip, rake, area, onset time.
        0.00000     0.00000     2000.00000     0.00000    0.00000     0.00000  1.0     0.00000
    source time function: delta, total sample point
        {dt}                  {len(vtime)}
    samples
    """
    fout.write(header)

    # Now write the moment rate:
    # check that integral is 1
    # print(np.trapz(momentRate, vtime))

    np.savetxt(fout, momentRate, fmt="%.18e")
    fout.close()
    print("done writing LOH1_source_50.dat")


if __name__ == "__main__":
    main()
