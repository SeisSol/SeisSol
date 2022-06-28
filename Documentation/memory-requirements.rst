.. _memory-requirements:

Memory requirements
~~~~~~~~~~~~~~~~~~~

General
-------

Memory requirements per node are difficult to predict because the elements are not distributed equally when Local Time-Stepping (LTS) is enabled.
However, an adequate rule of thumb to estimate memory requirements is to multiply the number of elements with the number of degrees of freedom per element and then multiply this number with a factor of 10.
Therefore, a run using the viscoelastic wave equation with 100 million elements of order 5 requires about 1.4 terabytes of memory.

When the equivalent distribution of memory among nodes is necessary, the user can choose the *exponential-balanced* node weight cost model for the initial partitioning of the workload.
See :ref:`workload-partitioning` for the available cost models.

