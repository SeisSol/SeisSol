Scaling
=======

When working with SI units in earthquake scenarios numbers might get very large.
Rescaling the equations might be advantageous, e.g. when working with single precision arithmetic.
In this section, we show how to properly scale the elastic wave equation.

The elastic wave equation in velocity-stress form with source terms is given by

.. math::

   \begin{aligned}
     \frac{\partial\sigma_{ij}}{\partial t}
           - \lambda \delta_{ij}\frac{\partial u_k}{\partial x_k}
           - \mu \left(
               \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}
            \right) &= s_i \\
     \rho\frac{\partial u_i}{\partial t} - \frac{\partial \sigma_{ij}}{\partial x_j} &= f_i
   \end{aligned}

We define scaled quantities

.. math::

   \bar{x}_i = \frac{x_i}{L}, \quad \bar{u}_i = \frac{u_i}{u_c}, \quad
       \bar{\sigma}_{ij} = \frac{\sigma_{ij}}{\sigma_c}

Inserting these into the elastic wave equation gives

.. math::

   \begin{aligned}
       \sigma_c\frac{\partial \bar{\sigma}_{ij}}{\partial t}
         - u_cL^{-1} \lambda \delta_{ij}\frac{\partial \bar{u}_k}{\partial\bar{x}_k}
         - u_c L^{-1} \mu \left(
            \frac{\partial \bar{u}_i}{\partial\bar{x}_j} +
            \frac{\partial \bar{u}_j}{\partial\bar{x}_i}\right) &= s_i \\
       \rho u_c\frac{\partial \bar{u}_i}{\partial t} - L^{-1}\sigma_c
         \frac{\partial \bar{\sigma}_{ij}}{\partial \bar{x}_j} &= f_i
   \end{aligned}

Multiplying the first equation with :math:`\sigma_c^{-1}`, the second with :math:`L\sigma_c^{-1}`,
and defining

.. math::

   \bar{\rho} = L\sigma_c^{-1} u_c\rho, \quad
   \bar{\lambda} = \sigma_c^{-1}u_cL^{-1}\lambda, \quad
   \bar{\mu} = \sigma_c^{-1}u_cL^{-1}\mu, \\
   \bar{f}_i = L\sigma_c^{-1} f_i, \quad
   \bar{s}_i = \sigma_c^{-1} s_i

leads to

.. math::

   \begin{aligned}
     \frac{\partial \bar{\sigma}_{ij}}{\partial t}
         - \bar{\lambda} \delta_{ij}\frac{\partial \bar{u}_k}{\partial\bar{x}_k }
         - \bar{\mu} \left(
            \frac{\partial \bar{u}_i}{\partial\bar{x}_j} +
            \frac{\partial \bar{u}_j}{\partial\bar{x}_i}\right) &= \bar{s}_i \\
     \bar{\rho}\frac{\partial \bar{u}_i}{\partial t } -
      \frac{\partial \bar{\sigma}_{ij}}{\partial \bar{x}_j } &= \bar{f}_i\end{aligned}

Example
-------
Set

.. math:: L = 10^3, \quad u_c = 1, \quad \sigma_c = 10^6

Then the spatial dimension of the mesh is [km], velocities are in [m/s], and stresses are in [MPa].
Parameters and source terms are scaled with

.. math::

   \bar{\rho} = 10^{-3}\rho, \quad
   \bar{\lambda} = 10^{-9}\lambda, \quad
   \bar{\mu} = 10^{-9}\mu, \quad
   \bar{f}_i = 10^{-3} f_i, \quad
   \bar{s}_i = 10^{-6} s_i

The wave speeds are given in [km/s]:

.. math::

   \bar{c}_s = \sqrt{\frac{\bar{\mu}}{\bar{\rho}}} = 10^{-3}c_s, \quad
       \bar{c}_p = \sqrt{\frac{\bar{\lambda} + 2\bar{\mu}}{\bar{\rho}}} = 10^{-3}c_p,
