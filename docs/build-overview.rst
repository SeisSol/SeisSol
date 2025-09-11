..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Building SeisSol
================

Welcome to the instructions for installing and running SeisSol. Therefore, this page serves as a high-level gateway to the individual steps of the build process.
To the current state, SeisSol has many dependencies due to a large amount of code generators, math libraries, and programming models (yes, all three of them).

Thus, we recommend using the Spack package; it can be found in the official Spack registry under ``seissol``.

Otherwise, if the spack package does not work, you can also install SeisSol on your own.

* Start by installing all necessary dependencies. As a warning, they may be quite a few in number. See the instructions under :ref:`Build Dependencies <build_dependencies>` on how to start.

* Then, compile SeisSol. See the instructions under :ref:`Build SeisSol <build_seissol>` for possible options.

Once you are done, you can test your installation by :ref:`running <build_run>` the SeisSol proxy, as well as some :ref:`examples <a_first_example>`.
