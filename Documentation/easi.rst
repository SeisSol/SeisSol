..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _easi:

easi
====

easi is a library for the Easy Initialization of model parameters in
three (or less) dimensional domains. easi offers the possibility to
parameterize the simulation without having to recompile SeisSol. Thanks
to easi, all user can run their simulations with the same executable (no
more hardcoded fortran initialization routines).

Writing easi models
-------------------

easi uses configuration files written in Yet Another Markup Language (YAML).
For example, a simple model with constant material parameters could be
described in the following way:

.. code-block:: YAML

  !ConstantMap
  map:
    lambda: 3.2044e+010
    mu:     3.2038e+010
    rho:    2670.
    Qp:     69.3
    Qs:     155.9

Complex models may be easily built, too.
For example, built-in with easi come `linear <https://easyinit.readthedocs.io/en/latest/maps.html#affinemap>`__ or `polynomial <https://easyinit.readthedocs.io/en/latest/maps.html#polynomialmap>`__ models
(e.g. to build depth-dependent models).
Furthermore, one may sample parameters from large two- or three-dimensional `uniform grids <https://easyinit.readthedocs.io/en/latest/maps.html#asagi>`__.
Lastly, code may be supplied written in a `C-like language <https://easyinit.readthedocs.io/en/latest/maps.html#functionmap>`__, which is compiled at run-time.

Please refer to `easi's documentation <https://easyinit.readthedocs.io/>`__
and to `easi's examples <https://github.com/SeisSol/easi/tree/master/examples>`__ for further details.

Invoking SeisSol
----------------

It is recommended to place easi models in the folder that also contains
the parameter file.
Within the parameter-file, add the parameter ``MaterialFileName`` to
the equations block, e.g.

.. code-block:: Fortran

  &equations
  MaterialFileName = 'material.yaml'
  /

When using :doc:`dynamic-rupture`, add the parameter ``ModelFileName`` to
the DynamicRupture block, e.g.

.. code-block:: Fortran

  &DynamicRupture
  ModelFileName = 'fault.yaml'
  /

Rheological model parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following parameters need to be set by easi.
The columms E, A, V, and P denote if the respective parameter is required
when using an (isotropic) elastic, anisotropic (elastic), viscoelastic, and viscoplastic rheological model.


.. |checkmark| unicode:: U+2713

.. list-table::
   :widths: 25 10 5 5 5 5 50
   :header-rows: 1

   * - Parameter
     - Unit
     - E
     - A
     - V
     - P
     - Description
   * - rho
     - :math:`\frac{kg}{m^3}`
     - |checkmark|
     - |checkmark|
     - |checkmark|
     - |checkmark|
     - Density.
   * - mu, lambda
     - Pa
     - |checkmark|
     -
     - |checkmark|
     - |checkmark|
     - Lamé parameters.
   * - c11, ..., c66 [#]_
     - Pa
     -
     - |checkmark|
     -
     -
     - stiffness tensor.
   * - Qp, Qs
     -
     -
     -
     - |checkmark|
     -
     - P-wave and S-wave quality factors.
   * - bulkFriction
     -
     -
     -
     -
     - |checkmark|
     - Bulk friction coefficient.
   * - plastCo
     - Pa
     -
     -
     -
     - |checkmark|
     - Plastic cohesion.
   * - s_xx, s_yy, s_zz, s_xy, s_yz, s_xz
     - Pa
     -
     -
     -
     - |checkmark|
     - Initial stress tensor.


.. [#] See :ref:`anisotropic` for more details.

Fault parameters (dynamic rupture)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following parameters need to be set by easi.
The column FL denotes for which friction law the respective parameter is required.
Please note that there are two ways to specify the initial stress on the fault:
You may either specify a stress tensor (s_xx, s_yy, s_zz, s_xy, s_yz, s_xz),
which has to be given for the same cartesian coordinate system as the mesh,
or you may specify a traction vector (T_n, T_s, T_d),
which has to be given in a fault local coordinate system.
You must not specify both.

.. list-table::
   :widths: 25 10 10 55
   :header-rows: 1

   * - Parameter
     - Unit
     - FL
     - Description
   * - s_xx, s_yy, s_zz, s_xy, s_yz, s_xz
     - Pa
     - all (excludes initial traction)
     - Initial stress tensor.
   * - T_n, T_s, T_d
     - Pa
     - all (excludes initial stress)
     - Initial traction in n=normal, s=strike, d=dip direction.
   * - cohesion
     - Pa
     - 2, 6
     - Magnitude of cohesive force.
   * - mu_s, mu_d
     -
     - 2, 6
     - Linear slip weakening: Static and dynamic friction coefficient.
   * - d_c
     - m
     - 2, 6
     - Linear slip weakening: Critical distance.
   * - forced_rupture_time
     - s
     - 16
     - Time of forced rupture.
   * - rs_a, rs_srW, RS_sl0
     -
     - 101, 103
     - Rate-and-state friction parameter.
   * - nuc_{xx, yy, zz, xy, yz, xz} or Tnuc_{n, s, d}
     - Pa
     - 2, 3, 4, 103
     - Nucleation stress or tractions.

Debugging easi script
---------------------


| Most easi components return easy to track error, for example
| ``test.yaml: yaml-cpp: error at line 6, column 9: illegal map value``
| Yet implajit function maps are more complex to debug. The following
  example:
| ``27.1: syntax error, unexpected '}', expecting ;``
| indicates that an error occurred in the 27th line of the function, but
  does not indicate which file and which function.
| Hopefully this will be improved in the future.


An example illustrating some subtleties of easi error logs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let suppose that we try to retrieve s_zz located at (x,y,z)=(0,0,0) in group 1 from the following easi file:

.. code-block:: YAML

    [s_zz,s_yy,s_yz,s_xx,s_xz,s_xy,d_c,mu_s]: !AffineMap
      matrix:
        xf: [0.4054811 , -0.91410343,  0.   ]
        yf: [-0.62424723, -0.2769057 ,  0.73050574]
        zf: [-0.6677578 , -0.29620627, -0.68290656]
      translation:
        xf: 348441.377459
        yf: 4760209.93637
        zf: 0.0
      components: !ASAGI
              file: norciax_210fault_nncia.nc
              parameters: [s_zz,s_yy,s_yz,s_xx,s_xz,s_xy,d_c,mu_s]
              var: data
              interpolation: nearest

and get the following error log:


.. code-block:: none

    terminate called after throwing an instance of 'std::runtime_error'
      what():  fault2.yaml@2: Could not find model for point [ 348441 4.76021e+06 0 ] in group 1.

How to interpret this error log?
The component at Line 2 is throwing the error (the AffineMap).
The AffineMap component is complaining that its output point is not accepted by any of its child components.
In this case, the point is outside the bounds of the ASAGI file.


Note that in the slightly different example below, without the AffineMap, easi will not verify that the point is outside the bounds of ASAGI file:

.. code-block:: YAML

    [s_zz,s_yy,s_yz,s_xx,s_xz,s_xy,d_c,mu_s]: !ASAGI
              file: norciax_210fault_nncia.nc
              parameters: [s_zz,s_yy,s_yz,s_xx,s_xz,s_xy,d_c,mu_s]
              var: data
              interpolation: nearest

In fact, in this case, ASAGI is directly queried and easi, therefore, does no verify that the point queried in inside the bounds of the ASAGI file.
If the point is out of bounds, ASAGI will pick the value of the nearest grid point and issue a warning:

.. code-block:: none

    Thu Jan 09 14:32:22, Warn:  ASAGI: Coordinate in dimension 2  is out of range. Fixing.


