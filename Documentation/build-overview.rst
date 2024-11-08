..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause

Building SeisSol
================

Welcome to the instructions for installing and running SeisSol. Therefore, this page serves as a high-level gateway to the individual steps of the build process.
To the current state, SeisSol has many dependencies due to a large amount of code generators, math libraries, and programming models (yes, all three of them).

Thus, we recommend using the Spack package; it can be found in the official Spack registry under ``seissol``.

Otherwise, if the spack package does not work, you can also install SeisSol on your own.

* Start by installing dependencies. See the instructions under :doc:`_build_dependencies`.

* Then, compile SeisSol. See the instructions under :doc:`_build_seissol`.

Once you are done, you can test your installation by running some examples.
