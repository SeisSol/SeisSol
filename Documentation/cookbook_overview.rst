Overview
========

This documentation is a collection of useful dynamic simulation examples to help users build models from the scratch with little efforts. Each example is demonstrated carefully with geometry building, parameter setup and result visualization. Users are suggested to repeat each example in order to get a comprehensive idea of how to set up dynamic simulation models with SeisSol. The updates of source software and news of SeisSol should be posted in SeisSol GitHub (https://github.com/SeisSol/SeisSol) or Readdocs (https://seissol.readthedocs.io/en/latest).

SeisSol is a part of SCEC dynamic code validation project (Harris et al. 2018) (http://scecdata.usc.edu/cvws/). Here we show several SCEC benchmarks for beginners to quickly catch up with SeisSol workflow. Each benchmark example is composed of a short problem description, a section of *geometry, initial setups (stress, nucleation, friction, etc.)*, and *simulation results*.

+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| No.            | Fault type              | Difficulty     | Description                                                                                                         |
+================+=========================+================+=====================================================================================================================+
| TPV5           | strike-slip             | beginner       | slip-weakening and heterogeneous initial stress conditions                                                          |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV6           | strike-slip             | beginner       | bi-material fault and, slip-weakening and heterogeneous initial stress conditions                                   |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV12          | normal fault            | beginner       | linear elastic and initial stress conditions are dependent on depth                                                 |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV13          | normal fault            | beginner       | non-associative Drucker-Prager plastic with yielding in shear ad initial stress conditions are dependent on depth   |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV16          | strike-slip             | intermediate   | randomly-generated heterogeneous initial stress conditions                                                          |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV24          | branching strike-slip   | intermediate   | a rightward branch forming a 30 degree angle. There are linear elastic material properties                          |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV29          | strike-slip             | difficult      | stochastic roughness. Linear elastic material properties in a homogeneous half-space.                               |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| TPV104         | strike-slip             | difficult      | Rate-state friction, using a slip law with strong rate-weakening.                                                   |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| Point Source   | strike-slip             | intermediate   | benchmark of SISMOWINE WP2\_LOH1.                                                                                   |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+
| Kinematic      | reverse fault           | intermediate   | Kinematic source of 1994 Mw6.7 Northridge earthquake.                                                               |
+----------------+-------------------------+----------------+---------------------------------------------------------------------------------------------------------------------+

Table: Overall of examples suites.

Prerequisites
~~~~~~~~~~~~~

Before you begin any of the examples, you will need to install latest
SeisSol from (https://github.com/SeisSol/SeisSol). The instruction can be found at https://seissol.readthedocs.io/en/latest/compilation.html. All geometry and
tetrahedral meshing are generated using free software Gmsh (http://gmsh.info/) is needed.
If you do not wish to create your own mesh at this time, the meshes are
also provided as part of the example. The ParaView visualization package
(https://www.paraview.org/) may be used to view simulation results. You may use other visualization
software, but some adaption from what is described here will be
necessary. Furthermore, you can complete a subset of the example using
files provided (as described below), skipping the steps for which you do
not have the proper software packages installed.

Input file resources
~~~~~~~~~~~~~~~~~~~~

The files needed to work through the examples are provided. Users may
download freely from this repository (https://github.com/daisy20170101/SeisSol_Cookbook).

References
~~~~~~~~~~~~~~~~~~~~

Harris, R. A., Michael Barall, B. T. Aagaard, S. Ma, and K. B. O. Daniel Roten, Benchun Duan, Dunyu Liu, Bin Luo, Kangchen Bai, Jean-Paul Ampuero, Yoshihiro Kaneko, Alice-Agnes Gabriel, Kenneth Duru, Thomas Ulrich, Stephanie Wollherr, Zheqiang Shi, Eric Dunham, Sam Bydlon, Zhenguo Zhang, Xiaofei Chen, Surendra N. Somala, Christian Pelties, Josue Tago, Victor Manuel Cruz-Atienza, Jeremy Kozdon, Eric Daub, Khurram Aslam, Yuko Kase, Kyle Withers (2018), A Suite of Exercises for Verifying Dynamic Earthquake Rupture Codes, Seismol. Res. Lett., 89(3), 1146-1162, doi:10.1785/0220170222.

