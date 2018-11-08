Introduction
============

This is a development version of SeisSol.

**SeisSol is still under heavy development and comes without any
guaranteed funcitonality. At the moment we can only provide very limited
support for general users. Please contact**\ `Alice
Gabriel <http://www.geophysik.uni-muenchen.de/Members/gabriel>`__\ **if
you are interested in a close collaboration.**

Folder Structure
----------------

It contains the following folders:

-  build/ Files for compiling the optimized generated kernel version
-  Maple/ Includes precomputed basis functions and other Maple tools;
   folder is essential to run SeisSol
-  postprocessing/

   -  science Tool box of Matlab, Python scripts for postprocessing
   -  visualization/

      -  DGVisu/ Contains source files for postprocessing the
         high-resolution output, readme included
      -  tools Tool box for all kind of visualization support

-  prepocessing/

   -  seissol_kernels submodule which is used to generate, test and tune
      SeisSol computational backend
   -  partitioning All kind of partitioning support files, including
      re-order approaches
   -  science Tool box of Matlab, Python scripts for preprocessing

-  src/ SeisSol source files are here
-  submodules/ Libraries
