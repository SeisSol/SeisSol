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

clear, close all, home;
disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%                    create Matterhorn                %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Parameters to Edit

nx =41;
ny =41;
load mmma50.xyz
a  = mmma50;
fid = fopen('matterhorn50.jou','w');

x = reshape(a(:,1),nx,ny); 
y = reshape(a(:,2),nx,ny); y = fliplr(y);
z = reshape(a(:,3),nx,ny); z = fliplr(z);

disp(sprintf('Minimum elevation:     %g',min(min(a(:,3)))));
disp(sprintf('Maximum elevation:     %g',max(max(a(:,3)))));
d = input('Give bottom elevation: ');

% Verschiebe die linke untere Ecke in den Koordinatenursprung

x0 = min(min(x(:,:)));
y0 = min(min(y(:,:)));

x(:,:) = x(:,:) - x0;
y(:,:) = y(:,:) - y0;

dx = abs(x(1)-x(2))*1;

fprintf(fid,'/ File opened for write %s.\n',datestr(now));
fprintf(fid,'identifier name "matterhorn" new nosaveprevious \n');

%Create surface vertices
counter = 0;
for j = 1:ny
    for i = 1:nx
        counter = counter + 1;
        fprintf(fid,'vertex create coordinates %g %g %g\n',x(i,j),y(i,j),z(i,j));
        vertices(i,j) = counter;
    end
end  
%Create bottom vertices
fprintf(fid,'vertex create coordinates %g %g %g\n',x(1,1),y(1,1),d);
fprintf(fid,'vertex create coordinates %g %g %g\n',x(nx,1),y(nx,1),d);
fprintf(fid,'vertex create coordinates %g %g %g\n',x(nx,ny),y(nx,ny),d);
fprintf(fid,'vertex create coordinates %g %g %g\n',x(1,ny),y(1,ny),d);

num_surf_edges = (nx-1)*2 + (ny-1)*2; 

%Create bottom edges
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+1,nx*ny+2);
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+2,nx*ny+3);
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+3,nx*ny+4);
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+4,nx*ny+1);

%Create vertical edges
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+1,vertices(1,1));
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+2,vertices(nx,1));
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+3,vertices(nx,ny));
fprintf(fid,'edge create straight "vertex.%i" "vertex.%i" real\n',nx*ny+4,vertices(1,ny));

%Create surface face
counter = 3;
fprintf(fid,'face create vertices');
for j = 1:ny
    for i = 1:nx
        counter = counter+1;
        fprintf(fid,' "vertex.%i"',vertices(i,j));
        if counter>=6
            fprintf(fid,' \\ \n');
            counter = 0;
        end
    end
end
fprintf(fid,'rowdimension %i  tolerance 0 \n',nx);

%Create side faces
fprintf(fid,'face create wireframe "edge.%i" "edge.%i" "edge.%i" "edge.%i" real\n',9,1,5,6);
fprintf(fid,'face create wireframe "edge.%i" "edge.%i" "edge.%i" "edge.%i" real\n',10,2,6,7);
fprintf(fid,'face create wireframe "edge.%i" "edge.%i" "edge.%i" "edge.%i" real\n',11,3,7,8);
fprintf(fid,'face create wireframe "edge.%i" "edge.%i" "edge.%i" "edge.%i" real\n',12,4,8,5);

%Create bottom faces
fprintf(fid,'face create wireframe "edge.%i" "edge.%i" "edge.%i" "edge.%i" real\n',1,2,3,4);

%Create Volume
fprintf(fid,'volume create stitch "face.1" "face.2" "face.3" "face.4" "face.5" "face.6" real\n');
fprintf(fid,'volume mesh "volume.1" tetrahedral size %g\n',dx/1);

%Create free surface boundary
fprintf(fid,'physics create "101" btype "ELEMENT_SIDE" face "face.1"\n');

%Create outflow boundary
fprintf(fid,'physics create "105" btype "ELEMENT_SIDE" face "face.2" "face.3" "face.4" "face.5" "face.6"\n');

fclose(fid);
   
       
