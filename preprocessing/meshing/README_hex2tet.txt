%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2005, SeisSol Group
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from this
%    software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%                                                     %%
   %%                   GAMBIT_HEX2TET                    %%
   %%           Creates regular tetrahedral mesh          %%
   %%           on the base of an hexahedral mesh         %%
   %%                                                     %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 Description of usage:

 GAMBIT_HEX2TET works only for domains that are meshed with a regular
 hexahedral mesh, i.e. the numbers nx, ny, and nz of hexahedral elements 
 in each space dimension is constant.

 Steps:

 (1) Create a hexagonal mesh with GAMBIT
 
 (2) Export the mesh as meshfile.neu  

 (3) Run GAMBIT_HEX2TET
     
     - specify the meshfile with the relativ path name 
       with respect to the path, where GAMBIT_HEX2TET.m is stored
       The suffix .neu is assumed as standard, so just use e.g.
       ../directory1/directory2/meshfile
     
     - specify the numbers nx, ny, and nz of hexahedral elements
       in x, y, and z dimension

     - this automatically produces a regular tetrahedral mesh
       where all hexahedral elements are subdivided into 5
       tetrahedrons and the new tetrahedral mesh is stored in 
       GAMBIT format in the same directory as meshfile.neu
       with the new name meshfile_tetra.neu 

 (4) Import the new meshfile_tetra.neu in GAMBIT for further
     processing, e.g. boundary conditions, ...


