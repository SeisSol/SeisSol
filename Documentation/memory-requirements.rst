..
  SPDX-FileCopyrightText: 2020-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Memory requirements
~~~~~~~~~~~~~~~~~~~

General
-------

Memory requirements per node are difficult to predict because the elements are not distributed equally when Local Time-Stepping (LTS) is enabled.
However, an adequate rule of thumb to estimate memory requirements is to multiply the number of elements with the number of degrees of freedom per element and then multiply this number with a factor of 10.
Therefore, a run using the viscoelastic wave equation with 100 million elements of order 5 requires about 1.4 terabytes of memory.


LTS weight balancing strategies
-------------------------------

By default, SeisSol uses a single-constraint node-weight mesh partitioning when LTS is enabled. The node-weights are evaluated as follows:

.. math::

   w_{k} = c_{k} R^{L - l_{k}}


where :math:`w_{k}` - node-weight of element :math:`k`; :math:`c_{k}` - cost of element :math:`k`, which depends whether 1) a cell is regular,
2) has :math:`n` dynamic rupture faces and 3) has :math:`m`  free surface with gravity faces; :math:`n, m \in [0, 4]`;  :math:`R` - LTS cluster update ratio;
:math:`L` - total number of LTS clusters; :math:`l_{k} \in [0, L)` - time cluster number, which element :math:`k` belongs to.

Because of the form of the node-weight function, we call this weight balancing strategy as *exponential*. The strategy reflects a computational
intensity of each element, taking element-wise update frequencies into account, and thus it aims to balance computations between MPI ranks.
However, the strategy may lead to memory imbalances, which can be a problem for systems with a limited amount of memory,
e.g. GPUs.


To address this issue, two other multi-constraint mesh partitioning strategies are available, namely: 1) *exponential-balanced* and 2) *encoded*.

*exponential-balanced* strategy can be described as follows:

.. math::

    w_{k} \in \mathbb{R}^{2} \mid
    w_{k} =
    \begin{bmatrix}
    c_{k} R^{L - l_{k}}\\
    1\\
    \end{bmatrix}

It features an additional constrain (the second component of :math:`w_{k}`) that
aims at equalizing the number of elements in each rank.

The *encoded* one:

.. math::

    w_{k} \in \mathbb{R}^{L} \mid
    w^{i}_{k} =
        \begin{cases}
            c_{k}, &  \text{if}\ i = l_{k} \\
            0, & \text{otherwise}
        \end{cases}
    & \text{and} \  i \in [0, L)

This strategy aims at distributing time clusters equally between all ranks
without taking into account element-wise update frequencies. The strategy may be
beneficial while working with LTS ratio 3 or 4.

A user can specify a particular partitioning strategy in *parameters.par* file:

.. code-block:: Fortran

    &Discretization
    ...
    ClusteredLTS = 2
    LtsWeightTypeId = 1  ! 0=exponential, 1=exponential-balanced, 2=encoded
    /


Note, the default (*exponential*) strategy is going to be used if *ClusteredLTS* is :math:`\geq 2` and
*LtsWeightTypeId* is not specified.
