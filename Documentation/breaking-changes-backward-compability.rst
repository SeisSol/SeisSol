Breaking changes in backwards compability
=========================================

SeisSol has been around since the early 2000s.
During this time, the code has seen a lot of enhancements and refactorings to always keep up-to-date with what compute-centers offer and geoscientists need.
Unfortunately, development of such a large software stack over such a long period of time, does not go without braking conventions from time to time.
In this page, we list breaking changes.

Convert parameter files for new DR release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While porting Dynamic Rupture to C++, we refactored a few parameter names to make things more consistent.
The new Dynamic Rupture implementation has been accecpted in August 2022 (commit XXXXXXX and later)

Parameter file (`parameters.par`):

+---------------+--------------------------+
| old           | new                      |
+===============+==========================+
| `0d0`         | `0.0`                    |
| `1d-16`       | `1e-16`                  |
| `v_star`      | `pc_vStar`               |
| `L`           | `pc_prakashLength`       |
| `mu_w`        | `rs_muW`                 |
| `alpha_th`    | `TP_thermalDiffusivity`  |
| `rho_c`       | `TP_heatCapacity`        |
| `tp_lambda`   | `TP_undrainedTPResponse` |
| `initemp`     | `TP_iniTemp`             |
| `inipressure` | `TP_iniPressure`         |
+---------------+--------------------------+

Fault material file (`fault.yaml`):

+-----------------------------+----------------------------+
| old                         | new                        |
+=============================+============================+
| `RS_sl0`                    |  `rs_sl0`                  |
| `alpha_hy`                  |  `tp_hydraulicDiffusivity` |
| `TP_half_width_shear_zone`  |  `tp_halfWidthShearZone`   |
+-----------------------------+----------------------------+

