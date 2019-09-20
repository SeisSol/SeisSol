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
The columms E, V, and P denote if the respective parameter is required
when using an elastic, viscoelastic, and viscoplastic rheological model.

.. |checkmark| unicode:: U+2713

.. list-table::
   :widths: 25 10 5 5 5 50
   :header-rows: 1

   * - Parameter
     - Unit
     - E
     - V
     - P
     - Description
   * - rho
     - :math:`\frac{kg}{m^3}`
     - |checkmark|
     - |checkmark|
     - |checkmark|
     - Density.
   * - mu, lambda
     - Pa
     - |checkmark|
     - |checkmark|
     - |checkmark|
     - Lamé parameters.
   * - Qp, Qs
     -
     - 
     - |checkmark|
     -
     - P-wave and S-wave quality factors.
   * - bulkFriction
     -
     - 
     -
     - |checkmark|
     - Bulk friction coefficient.
   * - plastCo
     - Pa
     - 
     -
     - |checkmark|
     - Plastic cohesion.
   * - s_xx, s_yy, s_zz, s_xy, s_yz, s_xz
     - Pa
     - 
     - 
     - |checkmark|
     - Initial stress tensor.

Fault parameters (dynamic rupture)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following parameters need to be set by easi.
The columm FL denotes for which friction law the respective parameter is required.
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
     - 2
     - Time of forced rupture.
   * - rs_a, rs_srW, RS_sl0
     - 
     - 101, 103
     - Rate-and-state friction parameter.
   * - nuc_{xx, yy, zz, xy, yz, xz}
     - Pa
     - 103
     - Nucleation stress.

Debugging easi script
---------------------

| Most easi components return easy to track error, for example
| ``test.yaml: yaml-cpp: error at line 6, column 9: illegal map value``
| Yet implajit function map are more complex to debug. The following
  example:
| ``27.1: syntax error, unexpected '}', expecting ;``
| indicates that an error occur in the 27th line of the function, but
  does not indicate which file and which function.
| Hopefully this will be improved in the future.
