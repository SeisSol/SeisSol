Breaking changes in backward compatibility
==========================================

To keep up-to-date with changes in compute-centers and geoscientists' needs, breaking changes sometimes needed to be introduced in SeisSol.
They are listed here.

Energy Output
~~~~~~~~~~~~~

Since we merged GitHub pull request `#531 <https://github.com/SeisSol/SeisSol/pull/531>`_ (April 2022), the seismic moment time history output, 
from which the moment rate can be post-processed, is integrated into the energy output  (see :ref:`energy_output`).
Therefore the parameters `magnitude_output_on`, `energy_rate_output_on` and `energy_rate_printtimeinterval` have been removed from the `DynamicRupture` namelist in the main parameter file.

C++ dynamic rupture implementation (dr/cpp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While porting dynamic rupture to C++, we changed a few parameter names to make things more consistent.
The new dynamic rupture implementation has been merged in September 2022 (GitHub pull request `#625 <https://github.com/SeisSol/SeisSol/pull/625>`_).
The linear slip weakening friction laws FL=2 (nucleation by stress increase) and FL=16 (forced time rupture nucleation) have been merged (the new friction law is FL=16).
Because of this change, FL=16 now requires nucleation stress or tractions to be specified in the fault-specific yaml file.

Parameter file (`parameters.par`):

+---------------+----------------------------------------------------------------------------+
| old           | new                                                                        |
+===============+============================================================================+
| `0d0`         | `0.0` (only until `#1173 <https://github.com/SeisSol/SeisSol/pull/1173>`_) |
+---------------+----------------------------------------------------------------------------+
| `v_star`      | `pc_vStar`                                                                 |
+---------------+----------------------------------------------------------------------------+
| `L`           | `pc_prakashLength`                                                         |
+---------------+----------------------------------------------------------------------------+
| `mu_w`        | `rs_muW`                                                                   |
+---------------+----------------------------------------------------------------------------+
| `alpha_th`    | `tp_thermalDiffusivity`                                                    |
+---------------+----------------------------------------------------------------------------+
| `rho_c`       | `tp_heatCapacity`                                                          |
+---------------+----------------------------------------------------------------------------+
| `tp_lambda`   | `tp_undrainedTPResponse`                                                   |
+---------------+----------------------------------------------------------------------------+
| `initemp`     | `tp_iniTemp`                                                               |
+---------------+----------------------------------------------------------------------------+
| `inipressure` | `tp_iniPressure`                                                           |
+---------------+----------------------------------------------------------------------------+

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
