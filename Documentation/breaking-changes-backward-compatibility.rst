..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Breaking changes in backward compatibility
==========================================

To keep up-to-date with changes in compute-centers and geoscientists' needs, breaking changes sometimes needed to be introduced in SeisSol.

All breaking changes for version 0.9.0 and later are listed here.

Energy Output
~~~~~~~~~~~~~
(since 0.9.0, `#531 <https://github.com/SeisSol/SeisSol/pull/531>`_, April 2022)

Since we merged GitHub pull request `#531 <https://github.com/SeisSol/SeisSol/pull/531>`_ (April 2022), the seismic moment time history output,
from which the moment rate can be post-processed, is integrated into the energy output (see :ref:`energy_output`).
Therefore, the parameters `magnitude_output_on`, `energy_rate_output_on` and `energy_rate_printtimeinterval` have been removed from the `DynamicRupture` namelist in the main parameter file.

C++ dynamic rupture implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(since 1.0.0, `#625 <https://github.com/SeisSol/SeisSol/pull/625>`_, September 2022)

While porting dynamic rupture to C++, we changed a few parameter names to make things more consistent.
The new dynamic rupture implementation has been merged in September 2022 (GitHub pull request `#625 <https://github.com/SeisSol/SeisSol/pull/625>`_).
The linear slip weakening friction laws FL=2 (nucleation by stress increase) and FL=16 (forced time rupture nucleation) have been merged (the new friction law is FL=16; but FL=2 will be accepted as input after `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_).
Because of this change, FL=16 now requires nucleation stress or tractions to be specified in the fault-specific yaml file.

Parameter file (`parameters.par`):

+---------------+-----------------------------------------------------------------------------------------------+
| old           | new                                                                                           |
+===============+===============================================================================================+
| `0d0`         | `0.0` (only until `#1173 <https://github.com/SeisSol/SeisSol/pull/1173>`_)                    |
+---------------+-----------------------------------------------------------------------------------------------+
| `v_star`      | `pc_vStar` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)               |
+---------------+-----------------------------------------------------------------------------------------------+
| `L`           | `pc_prakashLength` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)       |
+---------------+-----------------------------------------------------------------------------------------------+
| `mu_w`        | `rs_muW` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)                 |
+---------------+-----------------------------------------------------------------------------------------------+
| `alpha_th`    | `tp_thermalDiffusivity` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)  |
+---------------+-----------------------------------------------------------------------------------------------+
| `rho_c`       | `tp_heatCapacity` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)        |
+---------------+-----------------------------------------------------------------------------------------------+
| `tp_lambda`   | `tp_undrainedTPResponse` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_) |
+---------------+-----------------------------------------------------------------------------------------------+
| `initemp`     | `tp_iniTemp` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)             |
+---------------+-----------------------------------------------------------------------------------------------+
| `inipressure` | `tp_iniPressure` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)         |
+---------------+-----------------------------------------------------------------------------------------------+

Fault-specific yaml file (`fault.yaml`):

+-----------------------------+----------------------------+
| old                         | new                        |
+=============================+============================+
| `RS_sl0`                    |  `rs_sl0`                  |
+-----------------------------+----------------------------+
| `alpha_hy`                  |  `tp_hydraulicDiffusivity` |
+-----------------------------+----------------------------+
| `TP_half_width_shear_zone`  |  `tp_halfWidthShearZone`   |
+-----------------------------+----------------------------+
| `Ts0`                       |  `T_s`                     |
+-----------------------------+----------------------------+
| `Td0`                       |  `T_d`                     |
+-----------------------------+----------------------------+
| `Pn0`                       |  `T_n`                     |
+-----------------------------+----------------------------+

FORTRAN Removal
~~~~~~~~~~~~~~~

(since 1.1.0, `#829 <https://github.com/SeisSol/SeisSol/pull/829>`_, April 2023)

All FORTRAN had been removed; and several parameters have been deprecated.
However, all previous configuration files continue working.

Relative Paths
~~~~~~~~~~~~~~

(since 1.2.0, `#1156 <https://github.com/SeisSol/SeisSol/pull/1156>`_, August 2024)

SeisSol now looks for all paths in the parameter file first relatively to the parameter file,
and only then relative to the execution directory. This rule does not apply for the output directory.

Before this change, also all other paths were only viewed relative to the execution directory.

For absolute paths, nothing has changed.

DR Traction Computation
~~~~~~~~~~~~~~~~~~~~~~~

(since 1.1.0, `#895 <https://github.com/SeisSol/SeisSol/pull/895>`_, September 2023)

There was an adjustment in the Rate-and-State friction law computation; thus the results differ slightly between 1.0.1, and 1.1.0 onwards for the CPU computation.
The adjustment was propagated to the GPU implementation after 1.3.0.

The respective `commit <https://github.com/SeisSol/SeisSol/commit/73b284b7a8a2323170766f3ab594312a31f514c1>_`
updated the CPU implementation; the GPU implementation was updated by
`#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_.
