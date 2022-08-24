Breaking changes in backwards compability
=========================================

SeisSol has been around since the early 2000s.
During this time, the code has undergone a lot of enhancements and refactorings to always keep up-to-date with what compute-centers offer and geoscientists need.
Unfortunately, development of such a large software stack over such a long period of time, does not go without breaking conventions from time to time.
In this page, we list breaking changes.

Energy Output
~~~~~~~~~~~~~
Todo(TU, LK): I realized, I have to remove some lines regarding the energy output from the parameter file, after the refactoring of the energy output file.

Convert parameter files for latest DR implementation (dr/cpp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While porting Dynamic Rupture to C++, we changed a few parameter names to make things more consistent.
The new Dynamic Rupture implementation has been accecpted in September 2022 (github pull request `#625 <https://github.com/SeisSol/SeisSol/pull/625>`_).

Parameter file (`parameters.par`):

+---------------+--------------------------+
| old           | new                      |
+===============+==========================+
| `0d0`         | `0.0`                    |
+---------------+--------------------------+
| `1d-16`       | `1e-16`                  |
+---------------+--------------------------+
| `v_star`      | `pc_vStar`               |
+---------------+--------------------------+
| `L`           | `pc_prakashLength`       |
+---------------+--------------------------+
| `mu_w`        | `rs_muW`                 |
+---------------+--------------------------+
| `alpha_th`    | `tp_thermalDiffusivity`  |
+---------------+--------------------------+
| `rho_c`       | `tp_heatCapacity`        |
+---------------+--------------------------+
| `tp_lambda`   | `tp_undrainedTPResponse` |
+---------------+--------------------------+
| `initemp`     | `tp_iniTemp`             |
+---------------+--------------------------+
| `inipressure` | `tp_iniPressure`         |
+---------------+--------------------------+

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

Note that the identifiers in the `yaml` file are case sensitive.
Maybe you also have to change e.g. `RS_b` to `rs_b`.
