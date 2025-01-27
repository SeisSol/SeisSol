..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _Configuration:

Configuration
=============

To set up a SeisSol run, you need to do the following steps (assuming
you are in the root directory of the code repository):

1. Create the launch directory:

   ::

       mkdir launch_SeisSol

2. Copy all executables to the launch directory:

   ::

       cp build/SeisSol* launch_SeisSol/

3. Create your :doc:`parameter-file`

4. Copy any additional input files referenced in the parameter file (for
   example file with receiver coordinates) to your launch directory

5. (For large mesh only) Create symbolic links (ln -s) to the mesh
   file(s) in your launch directory

6. Make sure output and checkpoint directories exist

7. Optional: set :doc:`environment-variables` for tuning

Checklist for required files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Necessary files
^^^^^^^^^^^^^^^

-  SeisSol executable (compiled on the system where the job will run)
-  Parameter file
-  \*.yaml files for setting model parameters

Optional files depending on settings in the parameter file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  receiver files in \*.dat format
-  fault receiver files in \*.dat format (in the parameter file)
