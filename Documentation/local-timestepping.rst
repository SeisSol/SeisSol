..
  SPDX-FileCopyrightText: 2023-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Local time-stepping (LTS)
===================================

You can (and should) enable local time-stepping for your simulations.
Generally, this can lead to a better time to solution.
The following settings are relevant:

.. code-block:: Fortran

    &Discretization
    ...
    ClusteredLTS = 2

    LtsWiggleFactorMin = 0.51
    LtsWiggleFactorStepsize = 0.01
    LtsWiggleFactorEnforceMaximumDifference = 1

    LtsMaxNumberOfClusters = 20
    LtsAutoMergeClusters = 1
    LtsAllowedRelativePerformanceLossAutoMerge = 0.01
    LtsAutoMergeCostBaseline = 'bestWiggleFactor'
    /

To enable LTS, use the setting :code:`ClusteredLTS = 2`.
This is called rate-2 LTS.
Higher values are also supported, :code:`ClusteredLTS = N` leads to rate-N LTS (where N is a positive integer).
To disable LTS, use the value :code:`ClusteredLTS = 1`
For most simulations, we recommend to start with rate-2 LTS.

When using local time-stepping, SeisSol updates elements only as often as needed.
To do this, SeisSol computes a time step size for each element independently.
The elements are then grouped into so-called time-clusters, which are updated together.
The time step size of a cluster is the minimum of the time step sizes of the elements in the cluster.
Assuming rate-2 LTS (:code:`ClusteredLTS = 2`), the first cluster contains elements of time step sizes in the interval

.. math::

    [\lambda (\Delta t)^\text{min}, 2 \lambda (\Delta t)^\text{min}])

the second cluster elements of size

.. math::

    [2 \lambda (\Delta t)^\text{min}, 4 \lambda (\Delta t)^\text{min}])

and so on.
For LTS with higher rates, the base of the exponentiation changes.
The wiggle factor :math:`\lambda` is typically one.

Maximum Difference Property
----------------------------

SeisSol enforces some constraints on the clustering, for example, neighboring elements are always either in the same cluster,
or in a cluster which has at most a time step size difference of 2.
Elements that are connected by a dynamic rupture face have to be in the same time-cluster.
This is called the maximum difference property.


Wiggle factor (experimental)
----------------------------
This feature is only supported for rate-2 LTS (:code:`ClusteredLTS = 2`) at the moment.
The *LtsWiggleFactorMin* parameter sets the minimum allowed value for the wiggle factor :math:`0.5 < \lambda \leq 1`.
This wiggle factor can be used to reduce the overall number of time-steps in some cases and hence reduce the cost of the simulation.
Even though it seems counterproductive to reduce the global time step size, it can move the boundaries of the clusters such that
elements move from a cluster with a smaller time step size to clusters with a larger time step size.

SeisSol tries to find the optimal wiggle factor automatically by performing a grid search in the interval

.. math::

    \text{LtsWiggleFactorMin} \leq \lambda \leq 1

The step size of the grid search is controlled by the parameter *LtsWiggleFactorStepsize*.
The optimal wiggle factor is the one that minimizes the cost of updates per unit time.
When the setting *LtsWiggleFactorEnforceMaximumDifference* is set to one, SeisSol enforces the `Maximum Difference Property`_.
during the grid search. This leads to a better cost estimate, but computing this cost estimate can be costly for large scale simulations with many MPI ranks.
Hence, this is a trade-off between the time required to initialize the simulation and the time required to perform the actual simulation.
Typically this setting should be activated and only be deactivated when the initialization time becomes a bottleneck.

The wiggle factor was inspired by the implementation in (Breuer, Heinecke, 2022) [1]_

Enforcing maximum number of clusters (experimental)
----------------------------------------------------
You can set a maximum number of clusters by the parameter LtsMaxNumberOfClusters.
This can lead to better performance in some cases, especially on GPUs.
You can set a maximum number of clusters by setting :code:`LtsMaxNumberOfClusters=20`.
SeisSol also supports the automatic merging of clusters.
For this, you need to set *LtsAutoMergeClusters* to one.
The parameter *LtsAllowedRelativePerformanceLossAutoMerge* controls the allowed relative performance loss when merging clusters compared
to the baseline cost.
There are two different types of computing the baseline cost.
:code:`LtsAutoMergeCostBaseline = 'bestWiggleFactor'` allows merging clusters without losing the performance gained from the wiggle factor.
It computes the optimal wiggle factor without merging and then computes the best wiggle factor with auto-merging using the previous result as baseline.
This requires two iterations of finding the best wiggle factor but because the results of the most expensive operations are cached in the implementation, it is not much slower than :code:`LtsAutoMergeCostBaseline = 'maxWiggleFactor'`.
Alternatively, you can use :code:`LtsAutoMergeCostBaseline = 'maxWiggleFactor'`, which computes the cost without merging and wiggle factor and uses this as baseline cost.
The default and recommended choice is :code:`LtsAutoMergeCostBaseline = 'bestWiggleFactor'`.


These features should be considered experimental at this point in time.

.. [1] Breuer, A., & Heinecke, A. (2022). Next-Generation Local Time Stepping for the ADER-DG Finite Element Method. In 2022 IEEE International Parallel and Distributed Processing Symposium (IPDPS) (pp. 402-413). IEEE.
