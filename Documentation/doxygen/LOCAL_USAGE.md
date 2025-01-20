<!--
    SPDX-FileCopyrightText: 2021-2024 SeisSol Group

    SPDX-License-Identifier: BSD-3-Clause
    SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

    SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
-->

Local Doxygen Documentation
===========================

Install Doxygen
---------------

Fedora:

```bash
sudo dnf install doxygen
```

Ubuntu/Debian:

```bash
sudo apt install doxygen
```

Generate Documentation
----------------------

```bash
doxygen ./Doxyfile
```

Use your browser to view the documentation. For example:

```bash
firefox ./html/index.html
```
