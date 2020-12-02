Memory requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Memory requirements per node are difficult to predict because the elements are not distributed equally when local time-stepping is enabled.
However, an adequate rule of thumb to estimate memory requirements is to multiply the number of elements with the number of degrees of freedom per element and then multiply this number with a factor of 10.
Therefore, a run using the viscoelastic wave equation with 100 million elements of order 5 requires about 1.4 terabytes of memory.
