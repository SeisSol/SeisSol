..
  SPDX-FileCopyrightText: 2023 SeisSol Group

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

Varying-Rate Clustering
-----------------------

Instead of just writing :code:`ClusteredLTS = N`,
you can also input a list of numbers, e.g. :code:`ClusteredLTS = '3 2'`,
which defines the rate between the first clusters.
In the example, we have rate 3 between the first and the second cluster;
and rate 2 between all subsequent clusters. I.e.

.. math::

    [\lambda (\Delta t)^\text{min}, 3 \lambda (\Delta t)^\text{min}])

and then

.. math::

    [3 \lambda (\Delta t)^\text{min}, 6 \lambda (\Delta t)^\text{min}]).

The last number in the list is used for all subsequent clusters;
to keep consistency with the single-rate setups.
In particular using "1" last means that all clusters will be fused.
A rate of 1 is only accepted as the last entry of the list.

More precisely, a trailing 1 truncates the ladder: it does not create a further cluster, and all
elements that would have landed above simply stay in the topmost one.
A consequence worth knowing is that the entry just before the trailing 1 is never used.
:code:`ClusteredLTS = '3 2 5 6 1'` yields four clusters with ratios 3, 2 and 5 between them;
the 6 would only have mattered for a fifth cluster, which the trailing 1 forbids.

.. _maxdiff:

Cluster Smoothing (Maximum Difference Property)
-----------------------------------------------

SeisSol enforces some constraints on the clustering: neighboring elements are always either in the same cluster,
or in clusters that are *adjacent in the ladder*, i.e. whose indices differ by at most one.
Elements that are connected by a dynamic rupture face have to be in the same time-cluster.
This cluster smoothing is called the maximum difference property.

Note that the bound is on the cluster index and not on the time step size.
With uniform rate-2 LTS the two coincide, and neighbors differ by at most a factor of two.
With varying rates (see `Varying-Rate Clustering`_) the admissible ratio between neighbors is
whatever the local rate happens to be, so :code:`ClusteredLTS = '5 2'` allows a factor of five
across the boundary between the first two clusters.
The constraint is a hard requirement of the cluster construction, not a heuristic:
a cluster only ever exchanges data with its immediate neighbors in the ladder.
(the internal reason being: we would need to store more memory buffers if otherwise)

Wiggle factor (experimental)
----------------------------
The *LtsWiggleFactorMin* parameter sets the minimum allowed value for the wiggle factor :math:`0.5 < \lambda \leq 1`.
This wiggle factor can be used to reduce the overall number of time-steps in some cases and hence reduce the cost of the simulation.
Even though it seems counterproductive to reduce the global time step size, it can move the boundaries of the clusters such that
elements move from a cluster with a smaller time step size to clusters with a larger time step size.

SeisSol tries to find the optimal wiggle factor automatically by performing a grid search in the interval

.. math::

    \text{LtsWiggleFactorMin} \leq \lambda \leq 1

The step size of the grid search is controlled by the parameter *LtsWiggleFactorStepsize*.
The optimal wiggle factor is the one that minimizes the cost of updates per unit time.
When the setting *LtsWiggleFactorEnforceMaximumDifference* is set to one, SeisSol enforces the `Maximum Difference Property <_maxdiff>`_.
during the grid search. This leads to a better cost estimate, but computing this cost estimate can be costly for large scale simulations with many MPI ranks.
Hence, this is a trade-off between the time required to initialize the simulation and the time required to perform the actual simulation.
Typically this setting should be activated and only be deactivated when the initialization time becomes a bottleneck.

The wiggle factor was inspired by the implementation in (Breuer, Heinecke, 2022) [1]_

Enforcing maximum number of clusters (experimental)
---------------------------------------------------
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


.. _lts_lattice_search:

Choosing the ladder automatically (experimental)
------------------------------------------------

Everything described so far takes the rate vector from :code:`ClusteredLTS` as given and only tunes
the wiggle factor and the number of clusters.
Setting :code:`LtsClusteringSearch = 'lattice'` lets SeisSol choose the rate vector as well.

.. code-block:: Fortran

    &Discretization
    ...
    ClusteredLTS = 2
    LtsClusteringSearch = 'lattice'

    LtsCostUpdate = 1.0
    LtsCostLaunch = 0.0
    LtsCostFill = 0.0
    LtsMaxRate = 0
    /

The default, :code:`LtsClusteringSearch = 'legacy'`, is the grid search over wiggle factors with
merging from the top that is described in the sections above, and it reproduces the behavior of
earlier SeisSol versions exactly.

The cost model
~~~~~~~~~~~~~~

Both searches minimize a cost per unit of simulated time.
Writing :math:`W_k` for the summed vertex weight of cluster :math:`k` and :math:`n_k` for its update
factor (how many base time steps pass between two of its updates), that cost is

.. math::

    C = \sum_k \frac{c_\text{launch} + c_\text{update} \max(W_k, T_\text{fill})}
                     {n_k \lambda (\Delta t)^\text{min}}

with :math:`c_\text{update}` given by *LtsCostUpdate*, :math:`c_\text{launch}` by *LtsCostLaunch*
and :math:`T_\text{fill}` by *LtsCostFill*.

With the default values this is exactly the update counting that SeisSol has always used, namely
:math:`C = \sum_k W_k / (n_k \lambda (\Delta t)^\text{min})`.
The two additional terms model what counting updates misses, in particular on GPUs:

* :math:`c_\text{launch}` is paid whenever a cluster updates, no matter how little work it carries.
  It stands for kernel launches, the scheduler round trip and the communication handshake.
* :math:`T_\text{fill}` is the amount of work below which an update takes as long as an update
  carrying exactly that much work, because the device never saturates.

Both make a cluster that holds almost nothing expensive, which is what allows the search to decide
that a coarser ladder is worth it.
Sensible values depend on the machine and the backend and have to be measured; there is no
meaningful default beyond zero.

.. warning::

    With the default cost model the lattice search minimizes pure update count, and that objective
    is always best served by the finest ladder it is allowed to build.
    Switching to :code:`'lattice'` without setting *LtsCostLaunch* and/or *LtsCostFill* therefore
    mostly changes which rates are chosen, not how many clusters there are.
    SeisSol prints a note to the log in that case.

How the search works
~~~~~~~~~~~~~~~~~~~~

The update factors of a ladder cannot be arbitrary: a cluster may only advance once all finer
clusters have caught up, so :math:`n_0 = 1` must divide :math:`n_1`, which must divide
:math:`n_2`, and so on.
The set of admissible ladders is therefore the set of divisibility chains in
:math:`[1, n^\text{max}]`, and finding the cheapest one is a shortest path problem on the divisor
lattice, which is solved exactly by dynamic programming.

For every wiggle factor on the same grid the legacy search uses, SeisSol bins the element time step
sizes by :math:`\lfloor \Delta t_e / (\lambda (\Delta t)^\text{min}) \rfloor`, runs the dynamic
program on that histogram and keeps the cheapest candidate.
Only that one candidate is then realized and smoothed, so the expensive part of the search --
enforcing the `Maximum Difference Property <_maxdiff>`_ -- runs once rather than once per candidate.
Because the prediction ignores the demotions that smoothing causes, SeisSol logs the predicted and
the realized cost side by side.

*LtsMaxRate* bounds the ratio between two adjacent clusters; :code:`0` means unbounded.
This can be useful to keep the ladder gradual, at the price of a higher cost.
It has no effect for the legacy search, which does not choose the rates.

Interaction with the other settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* *LtsWiggleFactorMin* and *LtsWiggleFactorStepsize* control the wiggle factor grid for both searches.
* *LtsMaxNumberOfClusters* is honored by both.
* *LtsAutoMergeClusters*, *LtsAllowedRelativePerformanceLossAutoMerge* and *LtsAutoMergeCostBaseline*
  have no effect with :code:`'lattice'`, because choosing the number of clusters is part of what the
  dynamic program does.  SeisSol warns if they are set anyway.
* :code:`ClusteredLTS = 1` still disables LTS entirely, whichever search is selected.
  Otherwise :code:`ClusteredLTS` only provides the fallback ladder for the case where the search
  does not run at all, for instance because the material does not support LTS.

These features should be considered experimental at this point in time.

References
----------

.. [1] Breuer, A., & Heinecke, A. (2022). Next-Generation Local Time Stepping for the ADER-DG Finite Element Method. In 2022 IEEE International Parallel and Distributed Processing Symposium (IPDPS) (pp. 402-413). IEEE.
