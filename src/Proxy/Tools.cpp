// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Tools.h"

void print_hostname() {
  /*FILE* fp = popen("hostname", "r");
  if (fp > 0) {
    char buffer[256];
    fgets(buffer, 256, fp);
    strtok(buffer, "\n");
    printf("Running on %s.\n", buffer);
  }*/
}

auto sec(struct timeval start, struct timeval end) -> double {
  return (static_cast<double>(
             ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)))) /
         1.0e6;
}
