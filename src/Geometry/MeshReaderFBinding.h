/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

extern "C" {

// Functions implemented in Fortran
void allocelements(int n);
void allocvertices(int n, int maxElements);
void allocbndobjs(int n);
void allocbndobj(int i, int n);
void allocbndobjfault(int i, int n);
void allocfault(int n);

void hasplusfault();

void getverticesxy(int* size, double** verticesXY);
void getverticesnelements(int* size, int** nelements);
void getverticeselements(int* size, int** elements);

void getelementvertices(int* size, int** vertices);
void getreference(int* size, int** reference);
void getmpireference(int* size, int** mpireference);
void getmpinumber(int* size, int** mpinumber);
void getmpinumberdr(int* size, int** mpinumberdr);
void getboundarytoobject(int* size, int** boundarytoobject);
void getsideneighbor(int* size, int** sideneighbor);
void getlocalneighborside(int* size, int** localneighborside);
void getlocalneighborvrtx(int* size, int** localneighborvrtx);

void getbndsize(int* size);
void getbndrank(int i, int* rank);
void setbndrank(int i, int rank);
void getbndnelem(int i, int* nelem);
void setbndnelem(int i, int nelem);
void getbnddomainelements(int i, int* size, int** domainelements);

void getfaultreferencepoint(double* x, double* y, double* z, int* method);
void getfaultface(int* size, int** faultface);
void getfaultnormals(int* size, double** normals);
void getfaulttangent1(int *size, double** tangent);
void getfaulttangent2(int *size, double** tangent);

void setbndfaultnelem(int i, int nelem);
void getbndfaultelements(int i, int* size, int** faultelements);

}
