..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _energy_output:

Energy output
==============

Introduction
--------------

The energy output computes the energy of the simulation. It is divided into multiple parts:

- Energy in the water layer (= acoustic medium, :math:`\Omega_a`)
    - Gravitational energy
        :math:`\int_{\Omega_a} \frac{1}{2} \rho g \eta^2 \,\mathbf{dx}`, with :math:`\rho` the density, :math:`g` the gravitational acceleration and :math:`\eta` the sea-surface elevation.
    - Acoustic energy
        :math:`\int_{\Omega_a} \frac{1}{2K} p^2 \,\mathbf{dx}`, with :math:`p` the acoustic pressure and :math:`K` the compressibility.
    - Acoustic kinetic energy
        :math:`\int_{\Omega_a} \frac{1}{2} \rho v^2 \,\mathbf{dx}`, with :math:`\rho` the density and :math:`v` the velocity.
- Energy in Earth :math:`\Omega_e`
    - Elastic kinetic energy
        :math:`\int_{\Omega_e} \frac{1}{2} \rho v^2 \,\mathbf{dx}`, with :math:`\rho` the density and :math:`v` the velocity.
    - Elastic energy
        :math:`\int_{\Omega_e} \frac{1}{2} \epsilon_{ij} \sigma_{kl} \,\mathbf{dx}`, with  :math:`\epsilon_{ij}` the strain tensor and :math:`\sigma_{ij}` the stress tensor. It reduces for isotropic materials to :math:`\int_{\Omega_e} \frac{1}{2\mu} (\sigma_{ij} \sigma_{ij} -\frac{\lambda}{3\lambda+2\mu} \sigma_{kk}^2)\,\mathbf{dx}`, with :math:`\lambda` and :math:`\mu` the Lame coefficients.
- Earthquake source energy
    - Total frictional work done by the stress change
        :math:`W_\mathrm{total} = -\int_{0}^{t_f} \int_{\Sigma} \Delta\mathbf{\sigma}(t) \cdot \Delta\mathbf{\dot{u}}(t) \,\mathbf{dx}dt`, with :math:`\Sigma` the fault surface, :math:`\Delta\mathbf{\sigma}(t) = \mathbf{\sigma}(t) - \mathbf{\sigma}(0)` the shear traction change, and :math:`\Delta\mathbf{\dot{u}}(t)` the fault slip rate, and :math:`t_f` the end time of the simulation (see eq. 3 in Ma and Archuleta, 2006).
    - Static frictional work done by the stress change
        :math:`W_\mathrm{static} = -\int_{\Sigma} \frac{1}{2} \mathbf{\Delta\sigma}(t_f) \cdot \mathbf{\Delta u}(t_f) \,\mathbf{dx}` (see eq. 4 in Ma and Archuleta, 2006).
    - Radiated energy can then be computed with:
        :math:`E_\mathrm{r} = W_\mathrm{total} - W_\mathrm{static}` (see eq. 5 in Ma and Archuleta, 2006).

- Potency
        :math:`\int_{\Sigma} \Delta u_\mathrm{acc}(t_f) \,\mathbf{dx}`, with :math:`\Delta u_\mathrm{acc}` the accumulated fault slip (scalar).
- Seismic moment
        :math:`\int_{\Sigma} \mu \Delta u_\mathrm{acc}(t_f) \,\mathbf{dx}`, with :math:`\mu` the second Lame coefficient.
- Plastic moment
    :math:`\int_{\Omega_e} \mu \eta  \,\mathbf{dx}`, with :math:`\mu` the second Lame coefficient and \eta a scalar quantity measuring the accumulated material damage.

Currently, the output is only supported for the elastic wave equation.

Configuration
--------------

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
~~~~~~~~~~~~~~
Controlled via ``EnergyOutput``.
| 0 : no output
| 1 : csv output

For the example configuration, the output is written in the file "output/conv-energy.csv".

Terminal output
~~~~~~~~~~~~~~~~
Controlled via ``EnergyTerminalOutput``.
| 0 : no output
| 1 : additional output to stdout

Additionally, the energy can be written to stdout.
This can be useful for debugging. To increase the precision of the terminal output, adjust the ``EnergyTerminalPrecision`` field to the desired accuracy (a value of about 15 is sufficient to check for machine accuracy).

Note that ``EnergyOutput`` needs to be enabled for the terminal output to work.

Output interval
~~~~~~~~~~~~~~~~
The output interval is controlled by ``EnergyOutputInterval``.
If the output interval is not specified, the energy will be computed at the start of the simulation and at the end of the simulation.


Postprocessing and plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code below suggests a way to process and plot variables of the energy output file:

.. code-block:: python

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
