..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. Main SeisSol documentation file.
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

-----

We gratefully acknowledge the funding of the German Research Foundation (as part of project no. 391134334 - "CoCoReCS"), which massively contributed to creating all documentation, tutorials, example workflows and reproducible setups published on this website.

-----

.. toctree::
  :maxdepth: 2
  :caption: Introduction

  introduction
  acknowledge
  reproducible-research
  related-publications

.. toctree::
  :maxdepth: 2
  :caption: Installing SeisSol

  build-overview
  build-dependencies
  build-seissol
  gpus
  build-run
  a-first-example
  build-parameters
  build-archs
  build-problems

.. toctree::
  :maxdepth: 2
  :caption: Invoking SeisSol

  configuration
  parameter-file
  initial-condition
  local-timestepping
  left-lateral-right-lateral-normal-reverse.rst
  easi
  fault-tagging
  environment-variables
  memory-requirements

.. toctree::
  :maxdepth: 2
  :caption: SeisSol on Supercomputers

  behind_firewall
  supermuc
  shaheen
  frontera
  frontier
  marconi
  heisenbug
  leonardo
  lumi
  supermuc-ng-phase2

.. toctree::
  :maxdepth: 2
  :caption: Structural models and Meshing

  cad-models
  meshing-with-simmodeler
  meshing-with-pumgen
  gmsh
  asagi
  PUML-mesh-format

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

  io
  off-fault-receivers
  fault-output
  free-surface-output
  wave-field-output
  energy-output
  checkpointing
  postprocessing-and-visualization

.. toctree::
  :maxdepth: 2
  :caption: Further documentation

  sycl
  computing-time-vs-order-of-accuracy
  performance-measurement
  attenuation
  physical-models
  scaling
  basic-code-structure
  known-issues
  breaking-changes-backward-compatibility

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
