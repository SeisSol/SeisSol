%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2008, SeisSol Group
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
%
% @section DESCRIPTION
% Convergence Order computes convergence tables including L_1, L_2, and L_infty errors and the corresponding orders

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                             %%')
disp('     %% Convergence Order computes convergence tables including     %%')
disp('     %% L_1, L_2, and L_infty errors                                %%')
disp('     %% and the corresponding orders                                %%')
disp('     %%                                                             %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')
clear, close all;

%Specify filenames including the errors in the following ordering
% L_1
% L_2
% L_infty
%These three error norms are repeated n times for n Variables 
% (e.g. n=5 for 2D elastic or n=9 for 3D elastic)

load mesh04o3.dat
load mesh08o3.dat
load mesh12o3.dat
load mesh16o3.dat
load mesh20o3.dat
load mesh24o3.dat
load mesh28o3.dat
load mesh32o3.dat

%Combine all errors in one error matrix
error = [mesh04o3,...
         mesh08o3,...
         mesh12o3,...
         mesh16o3,...
         mesh20o3,...
         mesh24o3,...
         mesh28o3,...
         mesh32o3];

%Mesh spacing 
% (a) can be read from the output of the code, e.g.
h2 = [0.37267799, ...
      0.24845199, ...
      0.14907119, ...
      0.12422599, ...
      0.10647943, ...
      9.31694991e-2, ...
      8.281733e-2];

% (a) can be defined analytically for regular meshes, e.g. 
h = 2./[4 8 12 16 20 24 28 32]; 
h = ones(15,1)*h;

%computation of the order
n_row = size(error,1);
n_col = size(error,2);

order = log(error(:,2:n_col)./error(:,1:n_col-1))./log(h(:,2:n_col)./h(:,1:n_col-1));

disp(' '),disp('   The errors are:')
disp(error)
disp(' '),disp('   The orders are:')
disp(order)



