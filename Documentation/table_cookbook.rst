..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. list-table::
   :widths: 20 20 20 20 20 20
   :header-rows: 1

   * - Description
     - Source
     - Faulting |br| mechanism [#f1]_
     - Friction |br| law [#f2]_
     - Number |br| of faults
     - Further details

   * - `TPV5 <https://github.com/SeisSol/Examples/tree/master/tpv5>`_
     - dynamic
     - SS
     - LSW
     - 1
     - 3 stress asperities, see :ref:`tpv5`

   * - `TPV6 <https://github.com/SeisSol/Examples/tree/master/tpv6_7>`_
     - dynamic
     - SS
     - LSW
     - 1
     - Bi-material fault, heterogeneous initial stress, see :ref:`tpv6`

   * - `TPV12 <https://github.com/SeisSol/Examples/tree/master/tpv12_13>`_
     - dynamic
     - N
     - LSW
     - 1
     - depth-dependent initial stress conditions, see :ref:`tpv12`

   * - `TPV13 <https://github.com/SeisSol/Examples/tree/master/tpv12_13>`_
     - dynamic
     - N
     - LSW
     - 1
     - Same as TPV12 with non-associative Drucker-Prager plastic |br| with yielding in shear, see :ref:`tpv-13`

   * - `TPV16 <https://github.com/SeisSol/Examples/tree/master/tpv16>`_
     - dynamic
     - SS
     - LSW
     - 1
     - Randomly-generated heterogeneous initial stress conditions, |br| see :ref:`tpv16`

   * - `TPV24 <https://github.com/SeisSol/Examples/tree/master/tpv24>`_
     - dynamic
     - SS
     - LSW
     - 2
     - Rightward branch forming a 30 degree angle, see :ref:`tpv24`

   * - `TPV29 <https://github.com/SeisSol/Examples/tree/master/tpv29>`_
     - dynamic
     - SS
     - LSW
     - 1
     - Stochastic roughness, see :ref:`tpv29`

   * - `TPV34 <https://github.com/SeisSol/Examples/tree/master/tpv34>`_
     - dynamic
     - SS
     - LSW
     - 1
     - Imperial Fault model with 3D velocity structure, see :ref:`tpv34`

   * - `TPV104 <https://github.com/SeisSol/Examples/tree/master/tpv104>`_
     - dynamic
     - SS
     - fvw-RS
     - 1
     - see :ref:`tpv104`

   * - `LOH.1 <https://github.com/SeisSol/Examples/tree/master/WP2_LOH1>`_
     - point
     - n/a
     - n/a
     - n/a
     - point-source benchmark, see :ref:`loh1`

   * - `Northridge <https://github.com/SeisSol/Examples/tree/master/Northridge>`_
     - kinematic
     - R
     - n/a
     - 1
     - see :ref:`northridge`


.. [#f1] SS: strike-slip, N: normal, R: reverse, O: oblique
.. [#f2] LSW: linear slip-weakening friction, fvw-RS: fast-velocity weakening rate-and-state friction

.. |br| raw:: html

     <br>
