Dynamic Rupture
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
-  Cohesion = negative (acts on shear stress components, encapsulates
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


Projected linear slip-weakening
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A sudden decrease of slip rate to very small values, as for example occurring at a rupture front about to be arrested, may cause numerical issues and pollute the solution - in the worst case even leading to re-activation of rupture. 
Such numerical issues are easy to diagnose in the fault output visualization: a checkerboard pattern with unphysical values for the slip rate and slip in the affected region will be visible. 
This is a sign of mesh resolution (h-refinement) being locally not sufficient.
SeisSol mitigates such artifacts (currently only for linear slip-weakening friction) by projecting the slip increment on the 2D fault interface basis function after evaluation of the friction law. 
This implementation lowers the required h-refinement for linear-slip weakening friction across the fault in comparison to earlier implementations (cf. Pelties et al. , JGR 2012; GMD 2014)
Note, that the corresponding implementation for rate-and-state friction laws is currently pending - thus, a sudden decrease to small slip rates may cause numerical artifacts. 


Visualisation: SlipRateOutputType (default =1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the fault output will be showing regions affected by numerical problems. However, the user may choose to smooth out such artifacts for visualization purposes. Switching ``SlipRateOutputType`` in the ``DynamicRupture`` namelist from the default value to 0, will evaluate the slip-rate from the difference between the velocity on both sides of the fault, rather than evaluating the slip-rate from the fault tractions and the failure criterion directly. 
Note that this fix only applies to the output, but does not suppress numerical problems in the simulation.
Also note that that ``SlipRateOutputType=0`` is slightly less accurate than the default ``SlipRateOutputType=1`` without numerical problems. 

Friction laws
~~~~~~~~~~~~~

Linear-Slip Weakening Friction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rate-and-state Friction
^^^^^^^^^^^^^^^^^^^^^^^
Rate-and-state friction laws allow modeling the frictional shear strength variations as a function of slip rate and of the evolving properties of the contact population (Dieterich, 1979, 1981; Ruina, 1983).
In SeisSol, we currently support 3 types of rate-and-state friction laws, which differ by the set of ordinary differential equations describing the evolution of the state variable.
The type of rate-and-state friction law used is set by the FL variable in the DynamicRupture namelist (parameters.par):  

 
.. code-block:: Fortran

    3: ageing Law, 
    4: slip law, 
    103: rate-and-state with strong velocity weakening

More details about these friction law can be found in the `SCEC benchmarks descriptions <https://strike.scec.org/cvws/benchmark_descriptions.html>`_
(tpv101 to 105) or in Pelties et al. (2013, `GMD <https://gmd.copernicus.org/articles/7/847/2014/>`_).
Some parameters are considered homogeneous across the fault and defined in the main parameter file.
Others can spatially vary (rs_a, RS_sl0 for FL=3,4 and 103 and rs_srW for FL=103) and are defined in the fault yaml file.
Examples of input files for the `ageing law <https://github.com/SeisSol/Examples/tree/master/tpv101>`_
and for the `rate and state friction with strong velocity weakening <https://github.com/SeisSol/Examples/tree/master/tpv104>`_
are available at the given links.


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
  thermalPress = 1                  ! Thermal pressurization 0: inactive; 1: active
  IniTemp = 483.15                  ! Initial temperature [K]
  IniPressure = -80.0e6             ! Initial pore pressure; have to be added to normal stress in your initial stress yaml file [Pa]
  alpha_th = 1.0e-6                 ! Thermal diffusivity [m^2/s]
  rho_c = 2.7e6                     ! Specific heat [Pa/K]
  TP_lambda = 0.1e6                 ! Pore pressure change per unit temperature [Pa/K]

Two additional thermal pressurization parameters are space-dependent and therefore have to be specified in the dynamic rupture yaml file:

.. code-block:: YAML

  !ConstantMap
  map:
    alpha_hy: 1e-4                  ! Hydraulic diffusivity [m^2/s]
    TP_half_width_shear_zone: 0.01  ! Half width of shearing zone [m]

TP generates 2 additional on-fault outputs: Pore pressure and temperature (see fault output).

Slip-rate imposed on a DR boundary condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This friction law allows imposing slip-rate on a dynamic rupture boundary (potentially any kinematic models, but the current implementation is limited, see below).
The FL id for this friction law is 33.
The advantage of this approach compared to a multi point-sources representation is that the fault slip is not condensed to points. 
Therefore the discontinuity of the displacement across the fault can be accurately accounted for, and more generally the wavefield is accurate in the near-field.

The current implementation allows imposing a slip distribution on the DR Boundary using the same arbitrary smooth-step SR function everywhere on the fault.
The slip distribution is imposed simultaneously everywhere on the fault, smoothly over a time ``t_0``, where ``t_0`` is a parameter of the ``DynamicRupture`` namelist.
This is an expensive way of getting the final stress distribution from a given slip distribution.
The slip distribution is defined using easi by the ``strike_slip`` and ``dip_slip`` variables. 
Warning: the direction of positive ``strike_slip`` and ``dip_slip`` is based on the convention of Seissol (e.g. positive strike_slip for right-lateral faulting). 
Below is an example of an input file for defining the slip distribution:

.. code-block:: YAML

    !Switch
    [strike_slip, dip_slip]: !Any
      components:
       - !AxisAlignedCuboidalDomainFilter
          limits:
            x: [-1000, 0]
            y: [-1e10, 1e10]
            z: [-5000, -3000]
          components: !ConstantMap
            map:
              strike_slip:   0.01
              dip_slip: 0.1
       - !ConstantMap
          map:
            strike_slip: 0.05
            dip_slip: 0.05

