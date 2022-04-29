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
    - Total frictional work 
        :math:`\int_{0}^{t_f} \int_{\Sigma} \frac{1}{2} \Delta\mathbf{\sigma}(t) \cdot \Delta\mathbf{u}(t) \,\mathbf{dx}dt` with :math:`\Sigma` the fault surface, :math:`\Delta\mathbf{\sigma}(t)` the shear traction change, and :math:`\Delta\mathbf{u}(t)` the fault slip, and :math:`t_f` the end time of the simulation.
    - Static frictional work 
        :math:`\int_{\Sigma} \frac{1}{2} \mathbf{\Delta\sigma}(t_f) \cdot \mathbf{\Delta u}(t_f) \,\mathbf{dx}`.
- Seismic moment
        :math:`\int_{\Sigma} \frac{1}{2} \mu \Delta u_\mathrm{acc}(t_f) \,\mathbf{dx}`, with :math:`\mu` the second Lame coefficient and :math:`\Delta u_\mathrm{acc}` the accumulated fault slip (scalar).
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
    EnergyOutputInterval = 0.05
    /

Energy output
~~~~~~~~~~~~~~
| 0 : no output
| 1 : csv output

For the example configuration, the output is written in the file "output/conv-energy.csv".

Terminal output
~~~~~~~~~~~~~~~~
| 0 : no output
| 1 : additional output to stdout

Additionally, the energy can be written to stdout.
This can be useful for debugging.

Output interval
~~~~~~~~~~~~~~~~
The output interval is controlled by EnergyOutputInterval.
If the output interval is not specified, the energy will be computed at the start of the simulation and at the end of the simulation.
