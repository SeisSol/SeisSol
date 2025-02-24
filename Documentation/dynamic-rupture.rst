..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Dynamic rupture
===============

SeisSol is verified for a wide range of dynamic rupture problems
(Pelties et al. 2014): it is possible to use branched faults, dipping
fault geometries and laboratory-derived constitutive laws such as the
rate-and-state friction law.

Sign conventions
~~~~~~~~~~~~~~~~

The following definitions regarding fault and stress orientation are
conventions. The initial stresses need to be given in the global
coordinate system. The Matlab script
'/preprocessing/science/stressrotation3D' is available for correct
rotation of initial stress values on arbitrary orientated faults.

Definitions
^^^^^^^^^^^

-  Traction = (normal vector)\* (stress tensor)

   -  Shear traction in strike direction right (left) lateral: positive
      (negative)
   -  Shear traction along dip (anti-dip) direction: positive (negative)

-  Normal stress = negative in compression
-  Cohesion = negative (acts on normal stress, encapsulates
   the effect of pore pressurization)
-  Reference point = defines dip direction, starting point of the normal
   vector n pointing towards the fault, at +y
-  Note: SCEC benchmark descriptions define normal stress as positive in
   compression, opposite to SeisSol 3D convention.

Fault geometry
^^^^^^^^^^^^^^

The right-handed coordinate system is required. Please, make sure to follow
this convention during model/mesh generation!

.. ~ TODO: what's the point of these arrows?
.. ~ z y free-surface North ↑ ↗ ↑ ↗ ↑ ↗ ↑ ↗ ↑ → → → x = ↑ → → → East -z depth

-  Arbitrary fault orientation possible, **except** a fault in the
   xy-plane, i.e. parallel to the free surface

-  Fault Output is in fault coordinate system

   -  n = ( nx, ny, nz ).T normal vector, from + y → - y
   -  s = ( ny, - nx, 0).T \* sqrt ( nx^2 + ny^2) along-strike vector
   -  d = (-sy*nz, sx*\ nz, sy\ *nx-ny*\ sx).T / \|\| d \|\| along-dip
      vector, pointing in -z direction (downwards)

.. _example-1:-initial-stresses-for-a-60-dipping-fault:

Example 1: Initial stresses for a 60° dipping fault
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following the SCEC benchmark description of
`TPV10 <http://strike.scec.org/cvws/download/TPV10_11_Description_v7.pdf>`__
we need to transfer the following initial stresses given in fault
coordinates in global coordinates for SeisSol:

**INPUT** (from tpv10)

::

   normal stress = - 7378 Pa/m * (meter downdip-distance)
   shear stress = 0.55 * (normal stress) = + 4057.9 Pa/m * (meter downdip-distance)

Assuming the fault plane in xz-plane, dipping 60° in -y-direction we can
use the script '/preprocessing/science/stressrotation3D' to rotate the
stresses 30° clockwise at the x-axis (considering we look towards +x) to
get the full stress tensor in global coordinates:

**OUTPUT** (to be used in SeisSol's parameter file)

::

   s_yy = -  2019.26 Pa/m * (meter downdip-distance)
   s_zz = -  5358.74 Pa/m * (meter downdip-distance)
   s_yz = + 5223.72 Pa/m * (meter downdip-distance)
   s_xx = s_xy = s_xz = 0.0 Pa

.. _example-2:-initial-stresses-for-a-branched-fault:

Example 2: Initial stresses for a branched fault
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Normal and shear stresses in fault coordinate system for strike-slip
fault aligned with the x- or y-axis are already in the global
coordinates. But if we have additional branches the initial stresses
need to be rotated again. Following for example the SCEC benchmark
description for `TPV14 <http://strike.scec.org/cvws/download/TPV14_15_Description_v08.pdf>`__,
the stresses for the right-lateral strike-slip fault in xz-plane are:

::

   s_yy = - 120 MPa
   s_xy = + 70 MPa
   s_xy (nucleation patch) = + 81.6 MPa
   s_xx = s_zz = s_xz = s_yz = 0.0 Pa

For the 30° branch in -y direction we rotate the given stresses (in that case, they are the same as on the main fault) by 330° counterclockwise
(or 30° clockwise) around the z-axis using the script
'/preprocessing/science/stressrotation3D'.

**OUTPUT** (to be used in SeisSol's parameter file)

::

   s_xx = + 30.62 MPa
   s_yy = - 150.62 MPa
   s_xy =  - 16.96 MPa
   s_zz = s_xz = s_yz = 0.0 Pa

In case of a left-lateral strike slip fault (for example TPV15) the
stresses on the main fault are:

::

   s_yy = - 120 MPa
   s_xy =  - 70 MPa
   s_xy (nucleation patch) = - 81.6 MPa
   s_xx = s_zz = s_xz = s_yz = 0.0 Pa


Projecting the state variable increment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A sudden decrease of slip rate to very small values, as for example occurring at a rupture front about to be arrested, may cause numerical issues and pollute the solution - in the worst case even leading to re-activation of rupture.
Such numerical issues are easy to diagnose in the fault output visualization: a checkerboard pattern with unphysical values for the slip rate and slip in the affected region will be visible.
This is a sign of mesh resolution (h-refinement) being locally not sufficient.
SeisSol mitigates such artifacts by projecting the state variable (e.g. cumulative slip for linear slip weakening) increment on the 2D fault interface basis function after evaluation of the friction law.
This implementation lowers the required h-refinement across the fault in comparison to earlier implementations (cf. Pelties et al. , JGR 2012; GMD 2014)


Visualisation: SlipRateOutputType (default =1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the fault output will be showing regions affected by numerical problems. However, the user may choose to smooth out such artifacts for visualization purposes. Switching ``SlipRateOutputType`` in the ``DynamicRupture`` namelist from the default value to 0, will evaluate the slip-rate from the difference between the velocity on both sides of the fault, rather than evaluating the slip-rate from the fault tractions and the failure criterion directly.
Note that this fix only applies to the output, but does not suppress numerical problems in the simulation.
Also note that that ``SlipRateOutputType=0`` is slightly less accurate than the default ``SlipRateOutputType=1`` without numerical problems.

Friction laws
~~~~~~~~~~~~~

Linear slip-weakening friction (:code:`FL=6`, :code:`FL=16`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The linear slip-weakening friction is widely used for dynamic rupture simulations.

The fault strength is determined by

.. math::

  \tau = -C - \min\left(0, \sigma_n\right) \left( \mu_s - \frac{\mu_s - \mu_d}{d_c} \min\left(S, d_c\right)\right),

where :math:`S(t) = \int_0^t |V(s)| ds` is the accumulated fault slip, and the other variables are parameters of the friction, detailed below.

Friction parameters:

+------------------+----------------------------------------+-------------------------------+
| symbol           | quantity                               | SeisSol name                  |
+==================+========================================+===============================+
| :math:`\mu_s(x)` | static friction coefficient            | :code:`mu_s`                  |
+------------------+----------------------------------------+-------------------------------+
| :math:`\mu_d(x)` | dynamic friction coefficient           | :code:`mu_d`                  |
+------------------+----------------------------------------+-------------------------------+
| :math:`d_c(x)`   | slip-weakening critical distance       | :code:`d_c`                   |
+------------------+----------------------------------------+-------------------------------+
| :math:`C(x)`     | cohesion                               | :code:`cohesion`              |
+------------------+----------------------------------------+-------------------------------+
| :math:`T(x)`     | forced rupture time                    | :code:`forced_rupture_time`   |
+------------------+----------------------------------------+-------------------------------+
| :math:`v_0`      | threshold velocity                     | :code:`lsw_healingThreshold`  |
+------------------+----------------------------------------+-------------------------------+

Friction law :code:`16` implements linear slip-weakening with a forced rupture time.
If you are only interested in linear slip weakening friction without forced rupture time, do not supply the parameter `forced_rupture_time` in the fault `yaml` file.
Friction law :code:`6` uses Prakash-Clifton regularization for bimaterial faults.
For friction law :code:`16`, we resample the slip rate in every step to suppress spurious oscillations.
In the case of Prakash-Clifton regularization, we do not resample the slip rate.
If the slip rate :math:`V` drops below the threshold velocity :math:`v_0`, we reset the friction parameter :math:`\mu = \mu_s`.
Also, we reset the state variable :math:`S = 0`.
The threshold :math:`v_0` is set to :math:`-1.0` by default, such that healing is disabled.


Examples of input files for the friction laws :code:`6` and :code:`16` are availbable in the :ref:`cookbook<cookbook overview>`.

Linear slip weakening can be seen as a special case of rate-and-state friction with

.. math::
  \begin{aligned}
    f(V, \psi) &= C - \left( \mu_s - \frac{\mu_s - \mu_d}{d_c}\right) \min\left(\psi, d_c\right), \\
    g(V, \psi) &= V.
  \end{aligned}

Now the state variable stores the accumulated slip.

TP proxy slip-weakening friction (:code:`FL=1058`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The TP proxy slip-weakening friction has been proposed by Herrera et al. (2024), GJI, to approximate thermal pressurization in a computationally efficient way.
The fault strength is determined by

.. math::

  \tau = -C - \min\left(0, \sigma_n\right) \left( \mu_d + \frac{{(\mu_s - \mu_d)}}{{\left(1 + \frac{S}{d_c}\right)^{\alpha}}} \right),


All variables are the same as defined in previous section for :code:`FL=16`.
The friction law also supports forced rupture time.
You can modify the default value 1/3 of :math:`\alpha`, by adjusting the :code:`TpProxyExponent` parameter in the main parameter file (namelist: :code:`DynamicRupture`).

Rate-and-state friction
^^^^^^^^^^^^^^^^^^^^^^^
Rate-and-state friction laws allow modeling the frictional shear strength variations as a function of slip rate and of the evolving properties of the contact population (Dieterich, 1979, 1981; Ruina, 1983).
In SeisSol, we currently support 3 types of rate-and-state friction laws, which differ by the set of ordinary differential equations describing the evolution of the state variable.
The type of rate-and-state friction law is set by the FL variable in the DynamicRupture namelist (parameters.par):
Friction law :code:`3` implements the aging law, friction law :code:`4` implements the slip law, and friction law :code:`103` implements a slip law with strong rate-weakening.
More details about these friction laws can be found in the `SCEC benchmarks descriptions <https://strike.scec.org/cvws/benchmark_descriptions.html>`_ (TPV101 to 105) or in Pelties et al. (2013, `GMD <https://gmd.copernicus.org/articles/7/847/2014/>`_).

Some parameters are considered homogeneous across the fault and defined in the main parameter file.
Others can spatially vary (:code:`rs_a`, :code:`RS_sl0` for FL=3,4 and 103 and :code:`rs_srW` for FL=103) and are defined in the fault yaml file.
Examples of input files for the `aging law <https://github.com/SeisSol/Examples/tree/master/tpv101>`_
and for the `rate and state friction with strong velocity weakening <https://github.com/SeisSol/Examples/tree/master/tpv104>`_
are available at the given links.

All rate-and-state friction laws are described by the following system of differential algebraic equations, which depend on the state variable :math:`\psi` and the slip velocity :math:`V`.

.. math::

  \begin{aligned}
    \tau &= \sigma_n f(V,\psi) \\
    \frac{\partial\psi}{\partial t} &= g(V,\psi)
  \end{aligned}

Aging law (:code:`FL=3`)
-------------------------
Reference benchmarks: TVP101 and TPV102

Friction parameters:

+------------------+----------------------------------------+-------------------------------+
| symbol           | quantity                               | seisSol name                  |
+==================+========================================+===============================+
| :math:`a(x)`     | frictional evolution coefficient       | :code:`rs_a`                  |
+------------------+----------------------------------------+-------------------------------+
| :math:`b`        | frictional state coefficient           | :code:`rs_b`                  |
+------------------+----------------------------------------+-------------------------------+
| :math:`L(x)`     | characteristic slip scale              | :code:`rs_sl0`                |
+------------------+----------------------------------------+-------------------------------+
| :math:`V_0`      | reference slip velocity                | :code:`rs_sr0`                |
+------------------+----------------------------------------+-------------------------------+
| :math:`f_0`      | reference friction coefficient         | :code:`rs_f0`                 |
+------------------+----------------------------------------+-------------------------------+

.. math::
  \begin{aligned}
    f(V, \psi) &= a \sinh^{-1}\left[\frac{V}{2V_0} \exp\left( \frac{f_0 + b \ln(V_0 \psi / L)}{a}\right) \right] \\
    g(V, \psi) &= 1 - \frac{V \psi}{L}
  \end{aligned}

Slip law (:code:`FL=4`)
-----------------------
The slip law has the same parameters as the Aging Law.

.. math::
  \begin{aligned}
    f(V, \psi) &= a \sinh^{-1}\left[\frac{V}{2V_0} \exp\left( \frac{f_0 + b \ln(V_0 \psi / L)}{a}\right) \right] \\
    g(V, \psi) &= -V\frac{\psi}{L}\ln \left(V \frac{\psi}{L} \right)
  \end{aligned}

Severe velocity weakening (:code:`FL=7`)
----------------------------------------
No reference benchmark.

The Severe Velocity Weakening Law has the same parameters as the Aging Law does.

Strong velocity weakening (:code:`FL=103`)
------------------------------------------
Reference TPV103 and TPV104

In addition to the Aging and the Slip Law, strong velocity weakening requires two more parameters:

+------------------+----------------------------------------+-------------------------------+
| symbol           | quantity                               | seisSol name                  |
+==================+========================================+===============================+
| :math:`V_w(x)`   | weakening slip velocity                | :code:`rs_srW`                |
+------------------+----------------------------------------+-------------------------------+
| :math:`\mu_w`    | weakening friction coefficient         | :code:`rs_muW`                |
+------------------+----------------------------------------+-------------------------------+

.. math::
  \begin{aligned}
    f(V, \psi) &= a \sinh^{-1}\left[\frac{V}{2V_0} \exp\left(\frac{\psi}{a}\right) \right] \\
    g(V, \psi) &= - \frac{V}{L} \left(\psi - a \ln\left[ \frac{2V_0}{V} \sinh\left( \frac{\mu_{ss}(V)}{a} \right) \right] \right)
  \end{aligned}

with

.. math::
  \begin{aligned}
    \mu_{ss}(V) = \mu_w + \frac{f_0 - (b-a) \ln\left( \frac{V}{V_0} \right) - \mu_W}{\left( 1 + \left[ \frac{V}{V_W}\right]^8\right)^{1/8}}
  \end{aligned}.


Note that from the merge of pull request `#306 <https://github.com/SeisSol/SeisSol/pull/306>`__ of March 17th, 2021 to the merge of pull request `#752 <https://github.com/SeisSol/SeisSol/pull/752>`__ of December 22nd, 2022, the state variable was enforced positive in this friction law.
This enforcement aimed at avoiding the state variable getting negative because of Gibbs effects when projecting the state increment onto the modal basis functions (resampling matrix).
Since then, we realized that the state variable can get negative due to other factors, and, therefore, reverted this change.

Thermal Pressurization
~~~~~~~~~~~~~~~~~~~~~~

Seissol can account for thermal pressurization (TP) of pore fluids.
As deformation occurs within the fault gauge, frictional heating increases the temperature of the rock matrix and pore fluids.
The pore fluids then pressurize, which weakens the fault.
The evolution of the pore fluid pressure and temperature is governed by the diffusion of heat and fluid.
TP can be activated using ``thermalPress`` in the ``DynamicRupture`` namelist.
TP can be enabled with all rate-and-state friction laws (FL=3,4 and 103).
The TP parameters for which no spatial dependence has been implemented are defined directly in the ``DynamicRupture`` namelist:

.. code-block:: Fortran

  &DynamicRupture
  thermalPress = 1                     ! Thermal pressurization 0: inactive; 1: active
  tp_iniTemp = 483.15                  ! Initial temperature [K]
  tp_iniPressure = -80.0e6             ! Initial pore pressure; have to be added to normal stress in your initial stress yaml file [Pa]
  tp_thermalDiffusivity = 1.0e-6       ! Thermal diffusivity [m^2/s]
  tp_heatCapacity = 2.7e6              ! Specific heat [Pa/K]
  tp_undrainedTPResponse = 0.1e6       ! Pore pressure change per unit temperature [Pa/K]

Two additional thermal pressurization parameters are space-dependent and therefore have to be specified in the dynamic rupture yaml file:

.. code-block:: YAML

  !ConstantMap
  map:
    tp_hydraulicDiffusivity: 1e-4   # Hydraulic diffusivity [m^2/s]
    tp_halfWidthShearZone: 0.01     # Half width of shearing zone [m]

TP generates 2 additional on-fault outputs: Pore pressure and temperature (see fault output).

