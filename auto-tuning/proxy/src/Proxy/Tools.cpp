// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Tools.h"
#include <ctime>

auto derive_cycles_from_time(double time) -> double {
  // first try to read proxy env variable with freq
  /*char* p_freq;
  double d_freq;
  double cycles = 1.0;
  p_freq = getenv ("SEISSOL_PROXY_FREQUENCY");
  if (p_freq !=NULL ) {
    d_freq = atof(p_freq);
    printf("detected frequency (SEISSOL_PROXY_FREQUENCY): %f\n", d_freq);
    cycles = time * d_freq * 1.0e6;
  } else {
    FILE* fp;
    fp = popen("lscpu | grep MHz | awk '{print $3}'", "r");
    if(fp > 0) {
      char tmp_buffer[20];
      fread(tmp_buffer, 20, 1, fp);
      d_freq = atof(tmp_buffer);
      printf("detected frequency (lscpu): %f\n", d_freq);
      cycles = time * d_freq * 1.0e6;
      pclose(fp);
    } else {
      cycles = 1.0;
      printf("detected frequency (lscpu) FAILED!\n");
    }
  }
  return cycles;*/
  return 0;
}

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
  return ((double)(((end.tv_sec * 1000000 + end.tv_usec) -
                    (start.tv_sec * 1000000 + start.tv_usec)))) /
         1.0e6;
}
