Optimization for non Intel architectures
========================================

SeisSol is currently only optimized for Intel architectures. However, as
the single-node performance is mostly based on efficient small
matrix-matrix multiplications using
`libxsmm <https://github.com/hfp/libxsmm>`__, we believe that it should
be possible to obtain high performance on other architectures as well by
supplying a generator for efficient small matrix-matrix multiplications.

In this tutorial, we are going to show how to replace the calls to
libxsmm by a custom code generator, which will simply create calls to
Intel's Math Kernel Library (MKL). While this approach is particularly
useless for Intel architectures, it shows how one may replace libxsmm by
an alternative.

A simple code generator in bash
-------------------------------

During the build process, SeisSol calls the program
``libxsmm_gemm_generator`` in order to generate optimized GEMM routines,
where the term GEMM is borrowed from the BLAS level 3 routine with the
same name. This generator application requires up to 17 command line
arguments, which are described in `libxsmm's
documentation <https://libxsmm.readthedocs.io/en/latest/libxsmm_be/#generator-driver>`__.

As a first step, we create a bash script which mimics the interface of
libxsmm and creates calls to cblas_dgemm:

.. code:: bash

   #!/bin/bash

   if ! grep -q "gemm" ${2}; then
     echo "#include <mkl_cblas.h>" >> ${2}
   fi

   echo "void ${3}(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {" >> ${2}
   echo "  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ${4}, ${5}, ${6}, 1.0, A, ${7}, B, ${8}, ${11}, C, ${9});" >> ${2}
   echo "}" >> ${2}

We call this file ``libxsmm_gemm_generator``, we make it executable with
``chmod +x libxsmm_gemm_generator``, and adjust the ``PATH`` environment
variable to wherever this bash script is located with
``export PATH=/path/to/scripts/``. (You can also choose another file
name, then you have to `adjust the SCons
configuration <https://github.com/SeisSol/SeisSol/blob/201703/site_scons/site_tools/LibxsmmTool.py#L58>`__.)

Definition of your architecture
-------------------------------

Let's suppose you want to run SeisSol on Your Awesome Architecture
(YAA). YAA features 256 bit SIMD instructions which requires 32 bytes
memory alignment for optimal performance (aligned loads and stores). In
order to pass this information to SeisSol, we are going to create a new
architecture definition.

First, open the file ``generated_code/gemmgen/Arch.py`` and add
``'yaa': 32`` and ``'yaa': False`` in the ``alignments`` and
``enablePrefetch`` dictionaries.

Second, open the file ``site_scons/arch.py`` and add ``'yaa'`` in
``getArchitecture`` and ``'yaa': 32`` in the ``getAlignment`` function
(yes, it's redundant). In this file, you also have the option to set
architecture specific compiler flags. Here, we will use the same flags
as for the Haswell architecture and replace ``elif cpu == 'hsw':`` with
``elif cpu in ['hsw', 'yaa']:``. Inside the latter if statement, we add

.. code:: Python

   if cpu == 'yaa':
         flags.extend(['-L{0}/lib/intel64/'.format(os.getenv('MKLROOT')), '-lmkl_intel_ilp64', '-lmkl_sequential', '-lmkl_core', '-lm', '-ldl'])
         flags.extend(['-DMKL_ILP64','-I{0}/include'.format(os.getenv('MKLROOT'))])

In order to link against MKL. Don't forget to add ``import os`` at the
top of the file.

Third, as a final tedious step, we open the file
``src/Equations/elastic/generated_code/initialization/precision.h`` and
add

.. code:: C

   ## ifdef DYAA
   ## define DOUBLE_PRECISION
   ## endif

Performance evaluation
----------------------

You can use our performance proxy, whose performance usually closely
matches the performance of SeisSol, to test the single-node performance.
In order to compile it, change to the directory ``auto_tuning/proxy/``.
In there, run

.. code:: bash

   scons equations=elastic order=6 arch=dyaa memLayout=../config/elastic_dknl_O6.xml

The second option specifies the convergence order, feel free to vary the
number from 2 to 8. With arch we select our newly defined architecture.
The last option memLayout allows us to use sparse or block-partitioned
memory layouts (see ``../config/elastic_dhsw_O6.xml`` for an example).
However, as our code generator only supports dense matrix
multiplication, we use a memory layout with only dense matrix
multiplication. Note that in our experience sparse matrix multiplication
becomes less useful with increasing vector width; e.g. on Knights
Landing sparse matrix multiplication does not reduce the time to
solution.

The file ``build/generated_code/gemms.cpp`` contains your generated
kernels and allows you to check if any errors occurred. Furthermore, it
is a good idea to run the unit tests first: From ``auto_tuning/proxy/``,
run ``./build/bin/generated_kernels_test_suite`` on your cluster. The
result should be

.. code:: bash

   Running cxxtest tests (96 tests)..........................................OK!

Now you may finally evaluate your performance by writing a job script
for your cluster. Here, we are going to use the following LoadLeveler
script for `SuperMUC Phase
2 <https://www.lrz.de/services/compute/supermuc/systemdescription/>`__:

.. code:: bash

   #!/bin/bash
   #@ wall_clock_limit = 00:30:00
   #@ job_name = proxy_test
   #@ job_type = parallel
   #@ class = test
   #@ node = 1
   #@ tasks_per_node = 1
   #@ island_count = 1
   #@ node_usage = not_shared
   #@ initialdir = $(home)/seissol_reference/auto_tuning/proxy/
   #@ output = logs/proxy_test.$(schedd_host).$(jobid).out
   #@ error = logs/proxy_test.$(schedd_host).$(jobid).err
   #@ queue

   . /etc/profile
   . /etc/profile.d/modules.sh

   export OMP_NUM_THREADS=56
   export KMP_AFFINITY="compact,granularity=thread"

   ./build/bin/seissol_proxy 100000 100 all

The proxy application delivers the following information:

::

   time for seissol proxy              : 17.757648
   GFLOPS (non-zero) for seissol proxy : 146.865083
   GFLOPS (hardware) for seissol proxy : 422.893843

Hence, we are already at 43 % of peak performance (@ 2.2 GHz). However,
remember that we used dense GEMM routines instead of sparse GEMM
routines. Hence, the actual work done is better measured by time or by
non-zero GFLOPS.

Performance comparison
----------------------

We compare the performance obtained by noarch, yaa, and hsw on a
dual-socket Xeon E5-2697:

====== ======== ========= =========
arch   time [s] NZ-GFLOPS HW-GFLOPS
====== ======== ========= =========
noarch 78.7     33        92
yaa    17.8     147       423
hsw    10.8     240       554
====== ======== ========= =========

To summarise, using MKL already gave us a huge speed-up of 4.4 in
comparison to noarch. However, using a specialized code generator like
libxsmm in combination with auto-tuning we get another speed-up of 1.6.
