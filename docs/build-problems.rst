..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Frequent build and running issues
=================================

The following issues appear frequently when trying to compile or to run SeisSol.

Building fails due to missing files
-----------------------------------

Your submodules are probably not fully initialized.
Run ``git submodule update --init --recursive`` when inside the Git repository. In particular, that also works if you have your build directory inside the repository.
It is recommended to run this command after each repository update.

Code generation fails
---------------------

There is unfortunately no general remedy to that, as the problem may lie elsewhere.
For a fast, but stable installation, you can use ``PSpaMM``, ``gemmforge`` and ``chainforge``; all three of them being Python packages.

Then, you will need to set ``GEMM_TOOLS_LIST=PSpaMM`` to avoid any other code generators.

If nothing works, you can also try ``GEMM_TOOLS_LIST=Eigen`` without installing any additional generator. However, be aware that using Eigen
usually comes with a substantial performance penalty against small matrix code generators like libxsmm or PSpaMM; hence, we can only recommend it if nothing else works (or you just want to get a minimal build up and running and your meshes and/or simulations are small enough).

Running SeisSol gives ``SIGILL`` (reason 1)
-------------------------------------------

The underlying reason for this problem is usually that your selected ``HOST_ARCH`` in the CMake file and
the architecture you are running on do not match.

That can happen for example, if:

* Some version of the documentation recommends using ``HOST_ARCH=skx`` which enables some basic AVX512 options. However, most CPUs in personal computers, as well as AMD Epyc CPUs of Zen 3 or older do not support AVX512 to date. To avoid running into this problem, use ``HOST_ARCH=hsw``.
* Check again the system you want to run SeisSol on. Sometimes, for example, login nodes can have a different architecture compared to compute nodes.

To get a working build, choose ``noarch``. If that works, you are good on your local machine by using ``hsw`` (for x86_64 machines)
or ``neon`` (for ARM machines).

Running SeisSol gives ``SIGILL`` (reason 2)
-------------------------------------------

This problem also occurs, if you use the ImpalaJIT backend for easi on AARCH64-based CPUs (64-bit ARM), like e.g. in the latest Apple computers or the Nvidia Grace Hopper Superchip.
The crash usually happens then when reading material parameters, i.e. around the log messages of ``Computing LTS weights.``, or after ``Begin init model.``.

In this case (and also all other cases by now), it is best to use the Lua backend for easi instead.
The ImpalaJIT backend is used for ``!Function`` constructs while the Lua backend uses ``!Lua`` constructs. Since easi v1.5.0, a transpiler for ``!Function`` to ``!Lua`` constructs is included.
