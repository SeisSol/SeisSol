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

home;
disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%               GAMBIT_RECEIVERS_FAST                 %%')
disp('    %%           TO COMPUTE STATION ELEVATIONS             %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(' Uses a precomputed surface triangulation of a mesh-file and station ,')
disp(' coordinates in order to compute')
disp(' the exact elevation of the stations on the particular Gambit mesh,')
disp(' which approximates the topography with piecewise linear polynomials')
disp(' in the surface triangulation!')
disp(' '),disp(' ')
clear, close all;

filename = input('    Filename of mesh(suffix "_tri.mat" and "_vert.mat" are appended):  ','s');

disp('    Surface triangulation with vertices loaded!'); disp(' ')
eval(['load ',filename,'_tri.mat']);
eval(['load ',filename,'_vert.mat']);

rec_filename = input('    Filename of receiver locations(suffix ".dat" is appended):  ','s');
eval(['load ',rec_filename,'.dat']);
eval(['st = ',rec_filename,';']);

%Print mesh information on screen
disp(sprintf('\n\tFound %g surface triangles',size(tri,1)));
disp(sprintf('\tFound %g surface stations with xyz-coordinates:',size(st,1)));
figure; hold on
trisurf(tri,X(:,1),X(:,2),X(:,3),X(:,3)); axis equal,xlabel('x'),ylabel('y'),zlabel('z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          COMPUTE STATION ELEVATION          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = tsearch(X(:,1),X(:,2),tri,st(:,1),st(:,2));
for i = 1:size(st,1)
    %t(i) = tsearch(X(:,1),X(:,2),tri,st(i,1),st(i,2));
    v0 = X(tri(t(i),1),:);
    v1 = X(tri(t(i),2),:);  
    v2 = X(tri(t(i),3),:); 
    A = [1, v0(1), v0(2);1, v1(1), v1(2);1, v2(1), v2(2)];
    f = [v0(3);v1(3);v2(3)];
    c = A\f;
    est(i) = c(1)+c(2)*st(i,1)+c(3)*st(i,2); 
end
depth = input('    Specify depth to burry the receivers under the surface (in meters):  ');
receivers = [st(:,1),st(:,2),est'-depth];
plot3(st(:,1),st(:,2),est'-depth,'r*','MarkerSize',8)

disp('    Receiver coordinates under topography:'); disp(' ')
%disp(receivers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SAVE RECEIVER LOCATIONS TO FILE         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = input(['    Save receiver coordinates as ',filename,'_receivers.dat ? (y/n)'], 's');
if choice=='y'
    dlmwrite([filename ,'_receivers.dat'], receivers, ...
        'delimiter','\t','newline','unix','precision','%16.4f')
    disp('    Receiver coordinates saved!');
else
    disp('    Receiver coordinates NOT saved!');
end

disp('    FINISHED');
