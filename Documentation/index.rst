.. SeisSol documentation master file, created by
   sphinx-quickstart on Wed Nov  7 15:42:26 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======
SeisSol
=======

SeisSol is a software package for simulating wave propagation and dynamic
rupture based on the arbitrary high-order accurate derivative discontinuous
Galerkin method (ADER-DG).

Characteristics of the SeisSol simulation software are:

- use of arbitrarily high approximation order in time and space
- use of tetrahedral meshes to approximate complex 3D model geometries (faults & topography) and rapid model generation
- use of (an)isotropic elastic, viscoelastic and viscoplastic material to approximate realistic geological subsurface properties
- parallel geo-information input (ASAGI)
- to produce reliable and sufficiently accurate synthetic seismograms or other seismological data sets

.. toctree::
  :maxdepth: 2
  :caption: Introduction

  introduction
  compilation
  a-first-example
  acknowledge
  reproducible-research
  related-publications

.. toctree::
  :maxdepth: 2
  :caption: Structural models

  cad-models
  meshing-with-simmodeler
  meshing-with-pumgen
  gmsh

.. toctree::
  :maxdepth: 2
  :caption: Invoking SeisSol
  
  configuration
  parameter-file
  initial-condition
  left-lateral-right-lateral-normal-reverse.rst
  easi
  fault-tagging
  environment-variables
  gpus
  memory-requirements

.. toctree::
  :maxdepth: 2
  :caption: SeisSol on Supercomputers

  supermuc
  marconi

.. toctree::
  :maxdepth: 2
  :caption: Seismic source

  dynamic-rupture
  standard-rupture-format
  slip-rate-on-DR
  point-source-older-implementation

.. toctree::
  :maxdepth: 2
  :caption: Output
  
  fault-output
  free-surface-output
  off-fault-receivers
  postprocessing-and-visualization
  wave-field-output
  checkpointing
  energy-output
 
.. toctree::
  :maxdepth: 2
  :caption: Further documentation

  PUML-mesh-format
  asagi
  computing-time-vs-order-of-accuracy
  performance-measurement
  attenuation
  physical-models 
  scaling
  basic-code-structure
  known-issues

.. toctree::
  :maxdepth: 2
  :caption: Tutorials
  
  simmodelerCAD-workflow
  generating-a-cad-model-using-gocad-basic-tutorial
  generating-a-megathrust-geometry
  fully-coupled-mesh-tutorial
  remeshing-the-topography
  adapting-the-cad-model-resolution-using-gocad
  manually-fixing-an-intersection-in-gocad
  
.. toctree::
  :maxdepth: 2
  :caption: Cookbook

  cookbook_overview
  tpv5
  tpv6
  tpv12
  tpv13
  tpv16
  tpv24
  tpv29
  tpv34
  tpv104
  pointsource
  kinematic
  copyrights
