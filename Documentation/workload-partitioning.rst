.. _workload-partitioning:

Workload Partitioning 
~~~~~~~~~~~~~~~~~~~~~~

A computational mesh can be viewed as a graph composed of nodes (the volume cells) and edges (the cell faces), how this graph is initially distributed 
among computing nodes is important for the performance of any application. 
The initial workload partitioning (and the partitioning of the graph representation of the computational mesh)
of SeisSol can be configured by choosing one of the available Node Weight Cost Models and one of the Edge Weight Model Cost Models.

Local Time Stepping (LTS) Node weight Balancing Strategies
----------------------------------------------------------

By default, SeisSol uses a single-constraint node-weight mesh partitioning when LTS is enabled. The node weights are evaluated as follows: 

.. math::

   w_{k} = c_{k} R^{L - l_{k}}


where :math:`w_{k}` - node-weight of the element :math:`k`; :math:`c_{k}` - cost of element :math:`k`, which depends whether 1) a cell is regular, 
2) has :math:`n` dynamic rupture faces and 3) has :math:`m`  free surface with gravity faces; :math:`n, m \in [0, 4]`;  :math:`R` - LTS cluster update ratio;
:math:`L` - total number of LTS clusters; :math:`l_{k} \in [0, L)` - time cluster number, which element :math:`k` belongs to.

Because of the form of the node-weight function, we call this weight balancing strategy as *exponential*. The strategy reflects the computational 
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


Local Time Stepping (LTS) Edge Weight Balancing Strategies
-----------------------------------------------------------

*naive* is the default strategy where every edge is assigned to the same cost 1

.. math::

   w_{(a,b)} = 1

where :math:`(a,b)` denotes the directed edge from node a to b and :math:`w_{(a,b)}` the cost of the edge.

*communication approximation* is the adaptation of the exponential strategy for edges, it can be described as follows:

.. math::

   w_{(a,b)} = R^{L - l_{a}}


where :math:`(a,b)` denotes the directed edge from node a to b and :math:`w_{(a,b)}` the cost of the edge, :math:`R` - LTS cluster update ratio;
:math:`L` - total number of LTS clusters; :math:`l_{a} \in [0, L)` - time cluster number, which element :math:`a` belongs to.

Local Time Stepping (LTS)  Combined Weight Balancing Strategies
-----------------------------------------------------------------

These are specific node weight balancing strategies that have to be used with the communication approximation strategy and they extend the 
exponential-balanced strategy with additional constraints. SeisSol will throw a runtime error if they are not combined with communication approximation.

*communication approximation + minimum messaging* additionally tracks the amount of messages every node sends under the communication approximation cost model.
The edge weights are set to be then:

.. math::

    w_{k} \in \mathbb{R}^{3} \mid
    w_{k} = 
    \begin{bmatrix}
    c_{k} R^{L - l_{k}}\\
    1\\
    \sum_{(k,a) \in E} R^{L - l_{k}}\\
    \end{bmatrix}

It features an additional constraint (the third component of :math:`m_{k}`) that aims to track the amount of messages that are sent from these nodes.
:math:`(a,b)` denotes the directed edge from node a to b, :math:`R` - LTS cluster update ratio;
:math:`L` - total number of LTS clusters; :math:`l_{a} \in [0, L)` - time cluster number, which element :math:`a` belongs to.

*communication approximation + balanced messaging* Here :math:`L` additional constraints are added, where :math:`L` is the total number of LTS clusters. 
Each constraint represents the amount of nodes from the given time cluster the element :math:`k` communicates with.
With these additional constraints the partitioning aims to ensure that the amount of elements that send 
messages to each time cluster are evenly distributed among every partition. 


.. math::

    w_{k} \in \mathbb{R}^{2 + L} \mid
    w_{k} = 
    \begin{bmatrix}
    c_{k} R^{L - l_{k}}\\
    1\\
    |\{(k,a) \in E | l_a = 0\}| \\
    ... \\
    |\{(k,a) \in E | l_a = L - 1\}| \\
    \end{bmatrix}

:math:`R` is the LTS cluster update ratio;
:math:`L` - total number of LTS clusters; :math:`l_{a} \in [0, L)` - time cluster number, which element :math:`a` belongs to. The sets notate the edges that go 
out from the element :math:`k` and go to the time cluster :math:`l`. The cardinality of the set (the size of the set) is used for the value of the constraint. 

Choosing The Cost Models
----------------------

A user can specify a particular partitioning strategy in the *parameters.par* file:

.. code-block:: Fortran

    &Discretization
    ...
    ClusteredLTS = 2
    NodeWeightModelTypeId = 0            ! 0 for Exponential, 1 for Exponential Balanced, 2 for Encoded Balanced, 3 for  Minimum Messaging, 4 for Balanced Messaging
    EdgeWeightModelTypeId = 0            ! 0 for Naive, 1 for Communication Approximation,
                                         ! Node Weight Models 3 (Minimum Messaging) and 4 (Balanced Messaging) require the edge weight strategy to be set to 1 
    /


Note, the *exponential* strategy is going to be used per default if *ClusteredLTS* is :math:`\geq 2` and 
*LtsNodeWeightModelTypeId* is not specified, and the *naive* per default if the *EdgeWeightModelTypeId* is not specified.

*exponential-balanced* can be used when the memory of the compute nodes are limited. On the CPU *communication approximation* and 
*communication approximation + minimum messaging* may provide a slight performance boost over *naive*. The user is suggested to use
*exponential$* combined with *communication approximation* or  *communication approximation + minimum messaging* for CPUs.
Even though *communication approximation + balanced messaging* may provide a slight performance boost for the GPU in general, the results depend 
massively on the mesh, the problem size and the hardware and can often degrade the performance, therefore user is suggested to stick 
to the *naive* and *exponential-balanced* (as memory of the GPUs are often limited) unless the user wants to experiment with the different 
partitioning strategies.