..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Reproducible research
======================

On this page we list datasets containing reproducible simulation scenarios realized with SeisSol and elementary benchmark examples. In these data sets, all required input files, etc., are provided, to promote reproducible research and open science, and to inspire new simulation scenarios.

Please contact us via `SeisSol Email list <mailto:seissol-maintainers@lists.lrz.de>`_, if you're interested to provide a scenario. While many data sets are hosted on Zenodo, we do not require specific formats or repositories.

We provide datasets containing all input files of simulations in recent SeisSol publications as Zenodo repositories to promote reproducible research and open science.
A non-exhaustive list is given below:

.. list-table::
   :widths: 20 20 20 20 20 20 20
   :header-rows: 1

   * - Description
     - Source
     - Faulting |br| mechanism [#f1]_
     - Friction |br| law [#f2]_
     - Number |br| of faults
     - Research Item
     - Data

   * - 2004 Mw 9.2 Sumatra |br| megathrust
     - dynamic
     - R
     - LSW
     - 4
     - `Uphoff et al. (2017) <https://doi.org/10.1145/3126908.3126948>`_
     - `zenodo.org/record/439946 <https://zenodo.org/record/439946>`_

   * - Collective Knowledge |br| (CK) workflow
     - n/a
     - n/a
     - n/a
     - n/a
     - Fursin et al. (2018)
     - `zenodo.org/record/2422877 <https://zenodo.org/record/2422877>`_

   * - 2017 Mw 5.5 Pohang |br| induced earthquake
     - dynamic
     - O
     - fvw-RS
     - 2
     - `Palgunadi et al. (2020) <https://doi.org/10.1785/0120200106>`_
     - `zenodo.org/record/3930819 <https://zenodo.org/record/3930819>`_

   * - 2016 Mw 7.8 Kaikōura |br| earthquake
     - dynamic
     - mixed SS & R
     - fvw-RS
     - 10
     - `Ulrich et al. (2019) <https://doi.org/10.1038/s41467-019-09125-w>`_
     - `zenodo.org/record/2538024 <https://zenodo.org/record/2538024>`_

   * - 2018 Mw 7.5 Palu |br| earthquake
     - dynamic
     - mixed SS & N
     - fvw-RS
     - 3
     - `Ulrich et al. (2019b) <https://doi.org/10.1007/s00024-019-02290-5>`_
     - `zenodo.org/record/3234664 <https://zenodo.org/record/3234664>`_

   * - 2010 Mw 7.1 Darfield |br| earthquake
     - dynamic
     - SS
     - LSW
     - 6
     - `Uphoff (2019) <http://mediatum.ub.tum.de/?id=1531433>`_
     - `zenodo.org/record/3565774 <https://zenodo.org/record/3565774>`_

   * - Acoustic-Elastic model |br| of Palu earthquake
     - dynamic
     - mixed SS & N
     - fvw-RS
     - 3
     - `Krenz et al. (2021) <https://doi.org/10.1145/3458817.3476173>`_
     - `zenodo.org/record/5159333  <https://zenodo.org/record/5159333>`_

   * - Megathrusts initialized |br| with a 2D STM model
     - dynamic
     - R
     - LSW
     - 1
     - `Wirp et al. (2021) <https://doi.org/10.3389/feart.2021.626844>`_
     - `zenodo.org/record/4686551 <https://zenodo.org/record/4686551>`_

   * - Benchmarks for |br| poroelasticity
     - point
     - n/a
     - n/a
     - n/a
     - `Wolf et al. (2021) <https://doi.org/10.1016/j.jcp.2021.110886>`_
     - `zenodo.org/record/5236133 <https://zenodo.org/record/5236133>`_

   * - Low-angle normal fault |br| scenarios
     - dynamic
     - N
     - fvw-RS
     - 1
     - `Biemiller et al. (2022) <https://doi.org/10.1029/2021GC010298>`_
     - `zenodo.org/record/6094294 <https://zenodo.org/record/6094294>`_

   * - Megathrusts scenarios
     - dynamic
     - R
     - LSW
     - 4
     - `Madden et al. (2022) <https://doi.org/10.1029/2021JB023382>`_
     - `zenodo.org/record/5914661 <https://zenodo.org/record/5914661>`_

   * - 2016 Mw 6.2 Amatrice |br| broadband model
     - dynamic
     - N
     - LSW
     - 1
     - `Taufiqurrahman et al. (2022) <https://doi.org/10.1002/essoar.10510965.1>`_
     - `zenodo.org/record/6386938 <https://zenodo.org/record/6386938>`_

   * - 2016 Mw 6.5 Norcia |br| earthquake
     - dynamic
     - N
     - LSW
     - 1
     - `Tinti et al. (2021) <https://doi.org/10.1016/j.epsl.2021.117237>`_
     - `github repository <https://github.com/git-taufiq/NorciaMultiFault>`_

   * - 2004 Mw 9.2 Sumatra |br| megathrust
     - dynamic
     - R
     - LSW
     - 4
     - `Ulrich et al. (2022) <https://doi.org/10.1038/s41561-021-00863-5>`_
     - `zenodo.org/record/5541271 <https://zenodo.org/record/5541271>`_

..
   * - description
     - `xxx et al. (xxx) <https://doi.org/>`_
     - `zenodo.org/record/ <https://zenodo.org/record/>`_



SeisSol setups for community benchmark are described in the cookbook  (see :ref:`cookbook overview<cookbook overview>`), and the input files are available at https://github.com/SeisSol/Examples.

.. include:: table_cookbook.rst


We provide the following small-scale examples, specifically designed for SeisSol training and tutorials, such as the  `CHEESE Advanced training on HPC for Computational Seismology <https://www.hlrs.de/training/2021-10-19-cheese/>`_ and `ICTP Advanced Workshop on Earthquake Fault Mechanics <https:We provide the following small-scale examples, specifically designed for SeisSol training and tutorials//indico.ictp.it/event/8715/overview>`_ .  These SeisSol training examples are part of the `SeisSol Docker container <https://github.com/SeisSol/Training>`_  which also includes related open-source tools (Gmsh and ParaView) and all required input files.

.. list-table::
   :widths: 20 20 20 20 20 20
   :header-rows: 1

   * - Description
     - Source
     - Faulting |br| mechanism [#f1]_
     - Friction |br| law [#f2]_
     - Number |br| of faults
     - Data

   * - `TPV13 <https://github.com/SeisSol/Examples/tree/master/tpv12_13>`_
     - dynamic
     - N
     - LSW
     - 1
     - `<https://github.com/SeisSol/Training/tree/main/tpv13>`_

   * - 2018 Mw 7.5 Palu earthquake |br| (reduced mesh-size)
     - dynamic
     - mixed SS & N
     - fvw-RS
     - 3
     - `<https://github.com/SeisSol/Training/tree/main/sulawesi>`_

   * - `Northridge <https://github.com/SeisSol/Examples/tree/master/Northridge>`_
     - kinematic
     - R
     - n/a
     - 1
     - `<https://github.com/SeisSol/Training/tree/main/northridge>`_


