// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
