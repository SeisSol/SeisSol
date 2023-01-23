Local time-stepping (LTS)
===================================

You can (and should!) enable local time-stepping for your simulations.
This can lead to a better time to solution.
The following settings are relevant:

.. code-block:: Fortran

    &Discretization
    ...
    ClusteredLTS = 2

    LtsWiggleFactorMin = 1.0
    LtsWiggleFactorStepsize = 0.01
    LtsWiggleFactorEnforceMaximumDifference = 1
    LtsMaxNumberOfClusters = 20
    LtsAutoMergeClusters = 0
    LtsAllowedRelativePerformanceLossAutoMerge = 0.1
    /

To enable it, use the setting :code:`ClusteredLTS = 2`.
To disable it, use the value :code:`ClusteredLTS = 1`
Other values are currently not supported.

When using local time-stepping, SeisSol updates elements only as often as needed.
To do this, SeisSol computes a timestep size for each element independently.
The elements are then grouped into so-called time-clusters, which are updated together.
The timestep size of a cluster is the minimum of the timestep sizes of the elements in the cluster.
The first cluster contains elements of timestep sizes in the interval

.. math::

    [\lambda (\Delta t)^\text{min}, 2 \lambda (\Delta t)^\text{min}])

the second cluster elements of size

.. math::

    [2 \lambda (\Delta t)^\text{min}, 4 \lambda (\Delta t)^\text{min}])

and so on.
SeisSol enforces some constraints on the clustering, for example, neighboring elements are always either in the same cluster,
or a cluster which has at most a timestep size difference of 2.
Elements that are connected by a dynamic rupture face have to be in the same time-cluster.
THis is called the maximum difference property.

The factor :math:`\lambda` is typically one.

Wiggle factor (experimental)
----------------------------
The *LtsWiggleFactorMin* parameter sets the minimum allowed value for the wiggle factor :math:`0.5 < \lambda \leq 1`.
This wiggle factor can be used to reduce the overall number of time-steps in some cases.
Even though it seems counterproductive to reduce the global timestep size, it can move the boundaries of the clusters such that
elements move from a cluster with a smaller timestep size to a cluster with a larger timestep size.

SeisSol tries to find the optimal wiggle factor automatically by performing a grid search in the interval

.. math::

    \text{LtsWiggleFactorMin} \leq \lambda \leq 1

The step size of the grid search is controlled by the parameter *LtsWiggleFactorStepsize*.
The optimal wiggle factor is the one that minimizes the cost of updates per unit time.
When the setting *LtsWiggleFactorEnforceMaximumDifference* is set to one, SeisSol enforces the maximum difference property
during the grid search. This leads to a better cost estimate, but can be costly for large scale simulations with many MPI ranks.
The wiggle factor was inspired by the implementation in Breuer, Alexander, and Alexander Heinecke. "Next-Generation Local Time Stepping for the ADER-DG Finite Element Method.".

You can set a maximum number of clusters by the parameter LtsMaxNumberOfClusters.
This can lead to better performance in some cases.

SeisSol also supports the automatic merging of clusters.
For this, you need to set *LtsAutoMergeClusters* to one.
The parameter *LtsAllowedRelativePerformanceLossAutoMerge* controls the allowed relative performance loss when merging clusters compared
to the baseline cost without cluster-merging and wiggle factor.

These features should be considered experimental at this point.