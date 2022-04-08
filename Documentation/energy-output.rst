Energy output
============

Introduction
------------

The energy output computes the energy of the simulation. It is divided into multiple parts:

- Energy in the water layer (= acoustic medium).
    - Gravitational energy
    - Acoustic energy
    - Acoustic kinetic energy
- Energy in Earth
    - Elastic kinetic energy
    - Elastic energy
- Earthquake source energy
    - Total frictional work
    - Static frictional work
- Plastic moment

Currently, the output is only supported for the elastic wave equation.

Configuration
~~~~~~~~~~~

.. code-block:: Fortran

    &Output
    OutputFile = 'output/conv'
    EnergyOutput = 1
    EnergyTerminalOutput = 1
    EnergyOutputInterval = 0.05
    /

Energy output
---------------
| 0 : no output
| 1 : csv output

For the example configuration, the output is written in the file "output/conv_energy.csv".

Terminal output
---------------
| 0 : no output
| 1 : additional output to stdout

Additionally, the energy can be written to stdout.
This can be useful for debugging.

Output interval
--------------
The output interval is controlled by EnergyOutputInterval.
If the output interval is not specified, the energy will be computed at the start of the simulation and at the end of the simulation.