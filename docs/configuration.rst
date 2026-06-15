..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _Configuration:

Configuration
=============

To set up a SeisSol run, you need to do the following steps (assuming
you are in the root directory of the code repository):

1. Install SeisSol on your machine.

2. Create your :doc:`parameter-file`

3. Copy any additional input files referenced in the parameter file (for
   example file with receiver coordinates) to your launch directory.
   The paths in the parameter file and the easi files can be specified relatively.

4. (For large mesh only) Create symbolic links (``ln -s``) to the mesh
   file(s) in your launch directory.

5. Make sure output and checkpoint directories exist.

6. Optional (highly recommended on clusters): set :doc:`environment-variables` and :doc:`build-run` for tuning.

Checklist for required files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Necessary files
^^^^^^^^^^^^^^^

-  SeisSol executable / SeisSol installation (compiled on the system where the job will run)
-  Parameter file
-  \*.yaml files for setting model parameters

Optional files depending on settings in the parameter file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  receiver files in \*.dat format
-  fault receiver files in \*.dat format
