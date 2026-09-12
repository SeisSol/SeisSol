..
  SPDX-FileCopyrightText: 2022 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _energy_output:

Energy output
=============

Introduction
------------

The energy output computes the energy of the simulation. It is divided into multiple parts:

Energy in the water layer
~~~~~~~~~~~~~~~~~~~~~~~~~

The water layer is modelled as an acoustic medium over the domain :math:`\Omega_a`. The constitutive behaviour is governed by the relation :math:`p = -K \nabla \cdot \mathbf{u}`, with bulk modulus :math:`K`, pressure :math:`p`, and displacement :math:`\mathbf{u}`.

**Gravitational energy** is associated with the deformed sea surface :math:`\Gamma_\mathrm{free}` (the top boundary of :math:`\Omega_a`):

.. math::

   W_\mathrm{grav} = \int_{\Gamma_\mathrm{free}} \frac{1}{2} \rho\, g\, \eta^2 \,\mathrm{d}S

with :math:`\rho` the density, :math:`g` the gravitational acceleration and :math:`\eta` the sea-surface elevation. This formulation follows from the linearised free-surface boundary condition with gravitational restoring force [LottoDunham2015]_.

**Acoustic energy:**

.. math::

   W_\mathrm{ac} = \int_{\Omega_a} \frac{1}{2K}\, p^2 \,\mathrm{d}\mathbf{x}

with :math:`p` the acoustic pressure and :math:`K` the bulk modulus (compressibility).

**Acoustic kinetic energy:**

.. math::

   W_\mathrm{ac,kin} = \int_{\Omega_a} \frac{1}{2} \rho\, v_i\, v_i \,\mathrm{d}\mathbf{x}

with :math:`\rho` the density and :math:`v_i` the velocity.

**Total acoustic energy:**

.. math::

   W_\mathrm{ac,tot} = W_\mathrm{grav} + W_\mathrm{ac} + W_\mathrm{ac,kin}

**Dissipation.** The acoustic medium is modelled as inviscid and non-dissipative:

.. math::

   \dot{D}_\mathrm{ac} = 0

**Energy balance:**

.. math::

   \frac{\mathrm{d}}{\mathrm{d}t} W_\mathrm{ac,tot} = P_\mathrm{ext,ac} - \Phi_{\Gamma}

where :math:`P_\mathrm{ext,ac}` is the power of external sources in the acoustic domain and :math:`\Phi_{\Gamma}` is the energy flux through the acoustic-elastic coupling interface :math:`\Gamma`. At this interface, continuity of normal velocity and normal traction ensures that energy leaving the water layer enters the solid domain and vice versa:

.. math::

   \Phi_{\Gamma} = \int_{\Gamma} p\, v_i\, n_i \,\mathrm{d}S

with :math:`n_i` the outward normal of :math:`\Omega_a`. In a closed coupled system without external sources, the total energy :math:`W_\mathrm{ac,tot} + W_\mathrm{kin} + W_e` (summed over both domains) is conserved; any drift in the discrete balance indicates numerical dissipation or energy injection by the scheme.

Energy in the Earth
~~~~~~~~~~~~~~~~~~~

For comprehensive derivations of the energy expressions for elastic, anelastic, anisotropic and poroelastic media used throughout this section, see [Carcione2014]_. The Earth domain is denoted :math:`\Omega_e`.

Kinetic energy
^^^^^^^^^^^^^^

For single-phase solid material models (elastic, anisotropic, viscoelastic), the kinetic energy of the solid is:

.. math::

   W_\mathrm{kin} = \int_{\Omega_e} \frac{1}{2} \rho\, v_i\, v_i \,\mathrm{d}\mathbf{x}

with :math:`\rho` the density and :math:`v_i` the velocity. The poroelastic case is two-phase (solid skeleton + pore fluid) and requires additional terms for solid-fluid coupling and relative fluid motion; the corresponding expression is given in the poroelastic subsection below.

Elastic energy
^^^^^^^^^^^^^^

**General form.** The elastic (strain) energy is defined as:

.. math::

   W_e = \int_{\Omega_e} \frac{1}{2} \epsilon_{ij}\, \sigma_{ij} \,\mathrm{d}\mathbf{x}

with :math:`\epsilon_{ij}` the strain tensor and :math:`\sigma_{ij}` the stress tensor.

**Isotropic elastic materials.** For isotropic materials with constitutive law :math:`\sigma_{ij} = \lambda \delta_{ij} \epsilon_{kk} + 2\mu \epsilon_{ij}`, the strain can be expressed purely in terms of stress by inverting the constitutive relation:

.. math::

   \epsilon_{ij} = \frac{1}{2\mu}\left(\sigma_{ij} - \frac{\lambda}{3\lambda+2\mu}\,\sigma_{kk}\,\delta_{ij}\right)

Substituting into the energy expression yields:

.. math::

   W_e = \int_{\Omega_e} \frac{1}{4\mu} \left(\sigma_{ij}\, \sigma_{ij} - \frac{\lambda}{3\lambda+2\mu}\, \sigma_{kk}^2\right) \mathrm{d}\mathbf{x}

with :math:`\lambda` and :math:`\mu` the Lamé parameters. Written out explicitly (:math:`\sigma_{kk}^2` is shorthand for :math:`\sigma_{kk}\,\sigma_{ll} = (\mathrm{tr}\,\boldsymbol{\sigma})^2`):

.. math::

   W_e = \int_{\Omega_e} \frac{1}{4\mu}\Big(\sigma_{xx}^2 + \sigma_{yy}^2 + \sigma_{zz}^2 + 2\sigma_{xy}^2 + 2\sigma_{xz}^2 + 2\sigma_{yz}^2 - \frac{\lambda}{3\lambda+2\mu}\left(\sigma_{xx}+\sigma_{yy}+\sigma_{zz}\right)^2\Big)\,\mathrm{d}\mathbf{x}

**Anisotropic elastic materials.** For anisotropic materials with constitutive law :math:`\sigma_{ij} = c_{ijkl}\,\epsilon_{kl}`, the strain is obtained by introducing the compliance tensor :math:`s_{ijkl}`, defined as the inverse of the stiffness tensor on the space of symmetric second-order tensors:

.. math::

   c_{ijmn}\, s_{mnkl} = \tfrac{1}{2}\left(\delta_{ik}\,\delta_{jl} + \delta_{il}\,\delta_{jk}\right)

so that :math:`\epsilon_{ij} = s_{ijkl}\,\sigma_{kl}`. The elastic energy becomes:

.. math::

   W_e = \int_{\Omega_e} \frac{1}{2}\, s_{ijkl}\, \sigma_{ij}\, \sigma_{kl} \,\mathrm{d}\mathbf{x}

In Voigt notation, with the stress vector :math:`\boldsymbol{\sigma} = (\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz}, \sigma_{xy})^T`, the :math:`6 \times 6` Voigt stiffness matrix :math:`\mathbf{C}` (with components :math:`C_{IJ}` corresponding to :math:`c_{ijkl}` under the standard index pairing) and its inverse :math:`\mathbf{S} = \mathbf{C}^{-1}`, the energy reads:

.. math::

   W_e = \int_{\Omega_e} \frac{1}{2}\, \boldsymbol{\sigma}^T \mathbf{S}\, \boldsymbol{\sigma} \,\mathrm{d}\mathbf{x}

This relation uses the engineering convention in which the off-diagonal entries of the Voigt strain vector are :math:`2\epsilon_{ij}` (engineering shear strains); only under this convention does the matrix identity :math:`\mathbf{S} = \mathbf{C}^{-1}` reproduce the tensor relation :math:`\epsilon_{ij} = s_{ijkl}\,\sigma_{kl}`.

.. warning::

   The Voigt vector above is ordered :math:`(xx, yy, zz, yz, xz, xy)`, which is
   also the ordering the material parameters :math:`c_{IJ}` follow. SeisSol's
   *quantity* vector, and hence the stress moments the energy output works with,
   is ordered :math:`(xx, yy, zz, xy, yz, xz)` -- the last three are permuted.
   The compliance therefore has to be assembled in the quantity ordering before
   it is contracted with the stress moments. Getting this wrong is silent for any
   material with an isotropic shear block, and only shows up for genuinely
   anisotropic ones.

The isotropic result is recovered when :math:`c_{ijkl} = \lambda\,\delta_{ij}\,\delta_{kl} + \mu\,(\delta_{ik}\,\delta_{jl} + \delta_{il}\,\delta_{jk})`, the standard isotropic elasticity tensor in the Lamé parameters :math:`\lambda` and :math:`\mu`.

**Dissipation.** Both isotropic and anisotropic elastic materials are non-dissipative:

.. math::

   \dot{D} = 0

**Energy balance.** The total mechanical energy is exactly conserved:

.. math::

   \frac{\mathrm{d}}{\mathrm{d}t}\left(W_\mathrm{kin} + W_e\right) = P_\mathrm{ext}

where :math:`P_\mathrm{ext}` is the power of external sources (boundary tractions, body forces, fault slip). Without external sources:

.. math::

   W_\mathrm{kin}(t) + W_e(t) = W_\mathrm{kin}(0) + W_e(0) = \mathrm{const}

This conservation property makes the elastic energy balance a useful diagnostic for numerical accuracy: any drift in :math:`W_\mathrm{kin} + W_e` indicates numerical dissipation or energy injection by the scheme.


Viscoelastic energy
^^^^^^^^^^^^^^^^^^^

SeisSol implements the generalised Maxwell body with *strain-rate-like* anelastic
variables, which is worth stating explicitly because much of the literature uses
stress-like memory variables instead. The system solved is

.. math::

   \dot{\sigma}_{ij} = c_{ijkl}^U\,\dot{\epsilon}_{kl}
       - \sum_{l=1}^{L} d_{ijkl}^{(l)}\,\vartheta_{kl}^{(l)}

.. math::

   \dot{\vartheta}_{ij}^{(l)} = \omega_l\left(\dot{\epsilon}_{ij} - \vartheta_{ij}^{(l)}\right)

with the unrelaxed stiffness :math:`c^U` (this is what ``material.lambda`` and
``material.mu`` hold once ``fitAttenuation`` has run), the relaxation frequencies
:math:`\omega_l` (``material.omega``), and the modulus defect :math:`d^{(l)}` of
mechanism :math:`l`. The latter is stored in ``material.theta`` as

.. math::

   \theta^{(l)} = \bigl(-(\delta\lambda_l + 2\,\delta\mu_l),\; -\delta\lambda_l,\; -2\,\delta\mu_l\bigr)

so that :math:`d^{(l)}` is the isotropic stiffness built from
:math:`(\delta\lambda_l, \delta\mu_l)`, and the source tensor assembled by
``getTransposedSourceCoefficientTensor`` is exactly :math:`-d^{(l)}`. Note that
:math:`\delta\lambda_l` and :math:`\delta\mu_l` are fitted **independently** (from
two separate least-squares solves in ``fitAttenuation``); in particular
:math:`d^{(l)} \neq \omega_l\, c^U`, so formulations that assume a single scalar
weight per mechanism do not apply here.

The relaxed stiffness is

.. math::

   c^R = c^U - \sum_{l=1}^{L} d^{(l)},
   \qquad
   \lambda^R = \lambda^U - \sum_l \delta\lambda_l,
   \qquad
   \mu^R = \mu^U - \sum_l \delta\mu_l.

**The anelastic variables are strain rates.** The relaxation equation forces
:math:`[\vartheta] = [\dot{\epsilon}]`. Writing :math:`\xi^{(l)}` for their time
integral and :math:`e^{(l)} := \epsilon - \xi^{(l)}` for the spring strain of
Maxwell branch :math:`l`, one gets :math:`\dot{\xi}^{(l)} = \omega_l\, e^{(l)}`,
hence

.. math::

   e_{ij}^{(l)} = \frac{\vartheta_{ij}^{(l)}}{\omega_l}.

This is what makes the whole energy budget computable without any history: every
quantity below is a function of the instantaneous state
:math:`(\sigma, \vartheta^{(1..L)})` alone.

**Recovering the strain.** The total strain is not a state variable, but follows
from the relaxed constitutive relation:

.. math::

   \sigma^R_{ij} := \sigma_{ij} - \sum_{l=1}^{L} \frac{1}{\omega_l}\,
       d_{ijkl}^{(l)}\,\vartheta_{kl}^{(l)} \;=\; c^R_{ijkl}\,\epsilon_{kl}.

**Stored energy.** The free energy held in the relaxed spring and in the
:math:`L` Maxwell springs:

.. math::

   W_\mathrm{stored} = \frac{1}{2}\int_{\Omega_e}\left(
       c^R_{ijkl}\,\epsilon_{ij}\,\epsilon_{kl}
     + \sum_{l=1}^{L} d^{(l)}_{ijkl}\,e^{(l)}_{ij}\,e^{(l)}_{kl}
   \right)\mathrm{d}\mathbf{x}

The first term is reported as ``elastic_energy``, the second as
``viscoelastic_energy``. For :math:`L = 0` the first term reduces to the purely
elastic expression, as it must.

Both belong to the same reporting group as the kinetic energy, so the terminal
line reads

.. code-block:: none

   Elastic energy: <total>, kinematic x%, potential y%, viscoelastic z%

with the total being :math:`W_\mathrm{kin} + W_\mathrm{stored}`. Reporting only
the relaxed spring as "potential" would understate the stored energy by the
branch contribution, which grows with the attenuation.

.. warning::

   :math:`\tfrac{1}{2}\,\sigma_{ij}\,s^U_{ijkl}\,\sigma_{kl}` -- the elastic
   formula evaluated with the unrelaxed compliance -- is **not** the stored
   energy, and neither is the same expression with the relaxed compliance. The
   :math:`\sigma^R` correction cannot be absorbed into a choice of moduli; both
   variants are wrong by several percent, with a sign that depends on the loading
   history.

**Dissipation rate.** The power absorbed by the dashpots:

.. math::

   \dot{D} = \int_{\Omega_e} \sum_{l=1}^{L} \omega_l\;
       d^{(l)}_{ijkl}\,e^{(l)}_{ij}\,e^{(l)}_{kl}\,\mathrm{d}\mathbf{x}
   \;=\; \int_{\Omega_e} \sum_{l=1}^{L} \frac{1}{\omega_l}\;
       d^{(l)}_{ijkl}\,\vartheta^{(l)}_{ij}\,\vartheta^{(l)}_{kl}\,\mathrm{d}\mathbf{x}
   \;\geq\; 0

This is reported as ``viscous_dissipation_rate``, in watts. It is an
*instantaneous rate*, not a cumulative energy: the time integral

.. math::

   W_\mathrm{diss}(t) = \int_0^t \dot{D}(t')\,\mathrm{d}t'

is left to postprocessing, since accumulating it in the solver would require
carrying an extra state through the time-stepping kernels.

**Volumetric/deviatoric evaluation.** All of the above are isotropic quadratic
forms, so the implementation evaluates them in a volumetric/deviatoric split and
never inverts a :math:`6\times 6` matrix. With
:math:`\Delta K_l = \delta\lambda_l + \tfrac{2}{3}\delta\mu_l`:

.. math::

   \mathrm{tr}\,\sigma^R = \mathrm{tr}\,\sigma
       - \sum_l \frac{3\,\Delta K_l}{\omega_l}\,\mathrm{tr}\,\vartheta^{(l)},
   \qquad
   \mathrm{dev}\,\sigma^R = \mathrm{dev}\,\sigma
       - \sum_l \frac{2\,\delta\mu_l}{\omega_l}\,\mathrm{dev}\,\vartheta^{(l)}

.. math::

   d^{(l)}_{ijkl}\,a_{ij}\,a_{kl}
     = \Delta K_l\,(\mathrm{tr}\,a)^2 + 2\,\delta\mu_l\,\|\mathrm{dev}\,a\|^2

**Admissibility.** The least-squares fit in ``fitAttenuation`` does not constrain
the signs of its coefficients, so a poor combination of :math:`Q` and frequency
band can yield a fit that reproduces the target :math:`Q` while being
thermodynamically inconsistent (a branch with :math:`\delta\mu_l \leq 0`, or
non-positive relaxed moduli). The reported dissipation would then be negative.
SeisSol checks this once during setup and warns; the remedy is a different central
frequency or frequency ratio.

**Energy balance:**

.. math::

   \frac{\mathrm{d}}{\mathrm{d}t}\left(W_\mathrm{kin} + W_\mathrm{stored}\right)
       = P_\mathrm{ext} - \dot{D}

Without external sources, the total energy decreases monotonically due to
internal relaxation.


Poroelastic energy
^^^^^^^^^^^^^^^^^^

For poroelastic materials following Biot's theory [Biot1962]_, the constitutive equations are:

.. math::

   \sigma_{ij} = c_{ijkl}^d\,\epsilon_{kl} - \alpha\, p\, \delta_{ij}

.. math::

   p = M\left(\zeta - \alpha\, \epsilon_{kk}\right)

with drained stiffness :math:`c_{ijkl}^d`, Biot coefficient :math:`\alpha`, Biot modulus :math:`M`, pore pressure :math:`p`, and fluid content increment :math:`\zeta`.

**Elastic (potential) energy.** The total stored energy has a skeleton and a fluid contribution:

.. math::

   W_p = \frac{1}{2}\int_{\Omega_e}\left(\sigma_{ij}\,\epsilon_{ij} + p\,\zeta\right)\mathrm{d}\mathbf{x}

In quadratic form over the state variables :math:`(\epsilon_{ij},\,\zeta)`, substituting the constitutive relations gives the compact form:

.. math::

   W_p = \frac{1}{2}\int_{\Omega_e}\left(c_{ijkl}^d\,\epsilon_{ij}\,\epsilon_{kl} + M\,\bigl(\zeta - \alpha\,\epsilon_{kk}\bigr)^2\right)\mathrm{d}\mathbf{x}

Expressed purely in stresses and pore pressure, using the Biot effective stress :math:`\sigma'_{ij} = \sigma_{ij} + \alpha\, p\, \delta_{ij}` and the drained compliance :math:`s_{ijkl}^d` (defined analogously to :math:`s_{ijkl}` above, with :math:`c_{ijkl}^d` in place of :math:`c_{ijkl}`):

.. math::

   W_p = \frac{1}{2}\int_{\Omega_e}\left(s_{ijkl}^d\,\sigma'_{ij}\,\sigma'_{kl} + \frac{p^2}{M}\right)\mathrm{d}\mathbf{x}

For isotropic drained material (Lamé parameters :math:`\lambda_d, \mu`):

.. math::

   s_{ijkl}^d\,\sigma'_{ij}\,\sigma'_{kl} = \frac{1}{2\mu}\left(\sigma'_{ij}\,\sigma'_{ij} - \frac{\lambda_d}{3\lambda_d+2\mu}\,{\sigma'_{kk}}^2\right)

**Kinetic energy.** With solid velocity :math:`v_i^s`, Darcy velocity :math:`w_i = \phi(\dot{u}_i^f - \dot{u}_i^s)`, total density :math:`\rho = (1-\phi)\rho_s + \phi\rho_f`, and tortuosity tensor :math:`T_{ij}`:

.. math::

   W_\mathrm{kin} = \frac{1}{2}\int_{\Omega_e}\left(\rho\, v_i^s\, v_i^s + 2\rho_f\, v_i^s\, w_i + \frac{\rho_f}{\phi}\, T_{ij}\, w_i\, w_j\right)\mathrm{d}\mathbf{x}

The three terms represent the kinetic energy of the skeleton, the solid-fluid coupling, and the fluid motion relative to the skeleton, respectively.

**Dissipation rate.** Viscous dissipation due to Darcy friction:

.. math::

   \dot{D} = \int_{\Omega_e} \frac{\eta_f}{\kappa}\, w_i\, w_i \,\mathrm{d}\mathbf{x}

with dynamic fluid viscosity :math:`\eta_f` and permeability :math:`\kappa`. For anisotropic permeability, replace :math:`\frac{\eta_f}{\kappa}\, w_i\, w_i` by :math:`\eta_f\, \kappa_{ij}^{-1}\, w_i\, w_j`.

**Energy balance:**

.. math::

   \frac{\mathrm{d}}{\mathrm{d}t}\left(W_\mathrm{kin} + W_p\right) = P_\mathrm{ext} - \dot{D}

Without external sources, the total energy decreases monotonically due to viscous fluid friction.


Earthquake source energy
~~~~~~~~~~~~~~~~~~~~~~~~~

**Total frictional work** done by the stress change (cf. eq. 3 in [MaArchuleta2006]_):

.. math::

   W_\mathrm{total} = -\int_{0}^{t_f} \int_{\Sigma} \Delta\boldsymbol{\sigma}(t) \cdot \Delta\dot{\mathbf{u}}(t) \,\mathrm{d}S\,\mathrm{d}t

with :math:`\Sigma` the fault surface, :math:`\Delta\boldsymbol{\sigma}(t) = \boldsymbol{\sigma}(t) - \boldsymbol{\sigma}(0)` the shear traction change, :math:`\Delta\dot{\mathbf{u}}(t)` the fault slip rate, and :math:`t_f` the end time of the simulation.

**Static frictional work** done by the stress change (cf. eq. 4 in [MaArchuleta2006]_):

.. math::

   W_\mathrm{static} = -\int_{\Sigma} \frac{1}{2} \Delta\boldsymbol{\sigma}(t_f) \cdot \Delta\mathbf{u}(t_f) \,\mathrm{d}S

**Radiated energy** can then be computed as (cf. eq. 5 in [MaArchuleta2006]_):

.. math::

   E_\mathrm{r} = W_\mathrm{total} - W_\mathrm{static}


Other quantities
~~~~~~~~~~~~~~~~

**Potency:**

.. math::

   \int_{\Sigma} \Delta u_\mathrm{acc}(t_f) \,\mathrm{d}S

with :math:`\Delta u_\mathrm{acc}` the accumulated fault slip (scalar).

**Seismic moment:**

.. math::

   \int_{\Sigma} \mu\, \Delta u_\mathrm{acc}(t_f) \,\mathrm{d}S

with :math:`\mu` the shear modulus (second Lamé parameter).

**Plastic moment:**

.. math::

   \int_{\Omega_e} \mu\, \eta_p \,\mathrm{d}\mathbf{x}

with :math:`\mu` the shear modulus and :math:`\eta_p` a scalar measure of accumulated plastic strain (off-fault material damage).


Summary of energy balances
~~~~~~~~~~~~~~~~~~~~~~~~~~

+---------------------+-------------------------------------+---------------------------------------------+
| Material model      | Total energy                        | Dissipation mechanism                       |
+=====================+=====================================+=============================================+
| Acoustic            | :math:`W_\mathrm{grav} +            | None internally; energy exchanged via       |
|                     | W_\mathrm{ac} +                     | interface flux :math:`\Phi_\Gamma`          |
|                     | W_\mathrm{ac,kin}`                  |                                             |
+---------------------+-------------------------------------+---------------------------------------------+
| Elastic (iso/aniso) | :math:`W_\mathrm{kin} + W_e`        | None (energy conserved)                     |
+---------------------+-------------------------------------+---------------------------------------------+
| Viscoelastic        | :math:`W_\mathrm{kin} +             | Internal relaxation of the Maxwell          |
|                     | W_\mathrm{stored}`                  | branches, ``viscous_dissipation_rate``      |
+---------------------+-------------------------------------+---------------------------------------------+
| Poroelastic         | :math:`W_\mathrm{kin} + W_p`        | Viscous Darcy friction in the pore fluid,   |
|                     |                                     | ``darcy_dissipation_rate``                  |
+---------------------+-------------------------------------+---------------------------------------------+

In all dissipative cases the energy balance reads:

.. math::

   \frac{\mathrm{d}}{\mathrm{d}t}\left(W_\mathrm{kin} + W_\mathrm{pot}\right) = P_\mathrm{ext} - \dot{D}

where :math:`W_\mathrm{pot}` is the respective potential/stored energy and :math:`\dot{D} \geq 0` the dissipation rate.

.. note::

   The dissipation columns are reported as instantaneous **rates** (unit W), not
   as cumulative energies. Closing the balance therefore means integrating the
   reported rate over time in postprocessing and comparing against the change in
   :math:`W_\mathrm{kin} + W_\mathrm{pot}`. The accuracy of that integration is
   limited by the output interval, so a coarse ``EnergyOutputInterval`` will not
   close the budget to machine precision even for a perfectly conservative
   scheme.

For coupled acoustic-elastic simulations, the global energy balance over both domains reads:

.. math::

   \frac{\mathrm{d}}{\mathrm{d}t}\left(W_\mathrm{ac,tot} + W_\mathrm{kin} + W_\mathrm{pot}\right) = P_\mathrm{ext} - \dot{D}

The interface flux :math:`\Phi_\Gamma` cancels in the sum since the energy leaving one domain enters the other.


Configuration
-------------

.. code-block:: Fortran

   &Output
   OutputFile = 'output/conv'
   EnergyOutput = 1
   EnergyTerminalOutput = 1
   EnergyTerminalPrecision = 6
   EnergyOutputInterval = 0.05
   ComputeVolumeEnergiesEveryOutput = 4 ! Compute volume energies only once every ComputeVolumeEnergiesEveryOutput * EnergyOutputInterval
   /

Energy output
~~~~~~~~~~~~~

Controlled via ``EnergyOutput``.

| 0 : no output
| 1 : csv output

For the example configuration, the output is written in the file "output/conv-energy.csv".

Terminal output
~~~~~~~~~~~~~~~

Controlled via ``EnergyTerminalOutput``.

| 0 : no output
| 1 : additional output to stdout

Additionally, the energy can be written to stdout.
This can be useful for debugging. To increase the precision of the terminal output, adjust the ``EnergyTerminalPrecision`` field to the desired accuracy (a value of about 15 is sufficient to check for machine accuracy).

Note that ``EnergyOutput`` needs to be enabled for the terminal output to work.

Output interval
~~~~~~~~~~~~~~~

The output interval is controlled by ``EnergyOutputInterval``.
If the output interval is not specified, the energy will be computed at the start of the simulation and at the end of the simulation.

Postprocessing and plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code below suggests a way to process and plot variables of the energy output file:

.. code-block:: Python

   import pandas as pd
   import numpy as np
   import matplotlib.pylab as plt

   df = pd.read_csv("prefix-energy.csv")
   df = df.pivot_table(index="time", columns="variable", values="measurement")
   df["seismic_moment_rate"] = np.gradient(df["seismic_moment"], df.index[1])
   df.plot(y="seismic_moment_rate", use_index=True)

   # if ComputeVolumeEnergiesEveryOutput > 1
   volume_output = df.dropna()
   volume_output.plot(y="elastic_energy", use_index=True)

   plt.show()


References
----------

.. [Biot1962] Biot, M. A. (1962). Mechanics of deformation and acoustic propagation in porous media. *Journal of Applied Physics*, 33(4), 1482–1498. doi: `10.1063/1.1728759 <https://doi.org/10.1063/1.1728759>`_

.. [Carcione2014] Carcione, J. M. (2014). *Wave Fields in Real Media: Theory and Numerical Simulation of Wave Propagation in Anisotropic, Anelastic, Porous and Electromagnetic Media* (3rd ed.). Elsevier. ISBN 978-0-08-099999-9.

.. [LottoDunham2015] Lotto, G. C. and Dunham, E. M. (2015). High-order finite difference modeling of tsunami generation in a compressible ocean from offshore earthquakes. *Computational Geosciences*, 19(2), 327–340. doi: `10.1007/s10596-015-9472-0 <https://doi.org/10.1007/s10596-015-9472-0>`_

.. [MaArchuleta2006] Ma, S. and Archuleta, R. J. (2006). Radiated seismic energy based on dynamic rupture models of faulting. *Journal of Geophysical Research: Solid Earth*, 111, B05315. doi: `10.1029/2005JB004055 <https://doi.org/10.1029/2005JB004055>`_
