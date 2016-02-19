%%
% @file
% This file is part of SeisSol.
%
% @author Thomas Ulrich (ulrich AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/ulrich)
%
% @section LICENSE
% Copyright (c) 2005-2012, SeisSol Group
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


% Thomas Ulrich 01.2016
% The original routine reads a gambit neu mesh
% This amended version reads a netcdf file
home;
if verLessThan('matlab', '8.5')
    disp('YOUR MATLAB VERSION IS TOO OLD FOR THAT SCRIPT')
    return
end

disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%         PLACE_QUICKLY_FAULTRECEIVERS_NETCDF         %%')
disp('    %%         TO COMPUTE FAULTRECEIVER POSITIONS          %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(' Give netcdf mesh-file and station coordinates in order to compute')
disp(' the exact positions of the stations on the particular mesh,')
disp(' which approximates the fault with piecewise linear polynomials')
disp(' in the fault-surface triangulation!')
disp(' '),disp(' ')
clear, close all;
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
format long

% get input
filename0  =    input('    Filename of mesh(suffix ".nc" is appended)                  :  ','s');
normal       = input('    Normal axis to fault (x,y,z)                                :  ','s');
rec_filename = input('    Filename of 2D receiver locations(suffix ".dat" is appended):  ','s');
plotmesh     = input('    Do you wish to plot the mesh? (y,n)                         :  ','s');
if plotmesh == 'y'
  plotmarker = input('    Do you wish to plot the marker into the mesh? (y,n)         :  ','s');
else
  plotmarker = 'n';
end

% input check
if normal == 'x' || normal == 'y' || normal == 'z'
else
    disp('Wrong input for normal axis!')
    return
end

% load receiver stations
eval(['load ',rec_filename,'.dat']);
eval(['st = ',rec_filename,';']);

filename = sprintf('%s.nc',filename0);


element_size  = ncread(filename,'element_size');
[nPartition, ntrash] = (size(element_size));
element_vertices  = ncread(filename,'element_vertices');
element_boundaries  = ncread(filename,'element_boundaries');
vertex_coordinates  = ncread(filename,'vertex_coordinates');
vertex_size = ncread(filename,'vertex_size');

receivers=[];
if plotmesh == 'y'
   figure; hold on
end
disp(sprintf('\tFound %g faultreceiver with xyz-coordinates:',size(st,1)));

for iPartition=1:nPartition

   %read vertices
   X = vertex_coordinates(:,1:vertex_size(iPartition),iPartition);
   X=X';
   %read elements (+1 because numeration starts at 0 in netcdf)
   tetra = element_vertices(:,1:element_size(iPartition),iPartition)+1;
   tetra = tetra';
   %read elements boudary
   elemBound = element_boundaries(:,1:element_size(iPartition),iPartition);
   elemBound = elemBound';
   %find Boundary elements
   [tetraId, faceId] = find(elemBound==3);
   [NS, ntrash] = size(tetraId);
   if NS==0
      continue
   end

   tri = [];
   tri=zeros(NS,3);
   %GAMBIT convention
   %s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(3,:) = [2,3,4];   s_vert(4,:) = [1,4,3];
   %ICEM convention
   s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(4,:) = [2,3,4];   s_vert(3,:) = [1,4,3];

   for idBound = 1:NS
       tri(idBound,:) = tetra(tetraId(idBound), s_vert(faceId(idBound),:))';
   end



   if plotmesh == 'y'
       %figure; hold on
       trimesh(tri,X(:,1),X(:,2),X(:,3),X(:,3));
       axis equal,xlabel('x'),ylabel('y'),zlabel('z')
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                          COMPUTE STATION ELEVATION          
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Note, fault sides appear twice and therefore two faces will be found that
   % will be cleaned out later. This way, a fault station will be projected
   % also to branches lying behind the main fault.

   if normal == 'x'
      indexes= 2:3;
   elseif normal == 'y'
      indexes= 1:2:3;
   elseif normal == 'z'
      indexes= 1:2;
   end
   TR = triangulation(tri,X(:,indexes));
   t = pointLocation(TR,st(:,1:2));

   indicesInTriangulation = find(~isnan(t));
   ninside  = size(indicesInTriangulation,1);
   
      %Print mesh information on screen
   disp(sprintf('Partition %d:\tFound %g fault triangles,\t%d points inside the triangulation',iPartition, NS, ninside));

   t=t(indicesInTriangulation);
   newst=st(indicesInTriangulation,1:2);
   k=size(t);
   est=[];
   est(1:k) = 0.0;

   % linear interpolation
   for i=1:k  % add new stations generated by branch e.g.

       % load the three vertices for the receiver position
       v0 = X(tri(t(i),1),:);
       v1 = X(tri(t(i),2),:);
       v2 = X(tri(t(i),3),:);

       % generate interpolation matrix
       if normal == 'x'
	   A = [1, v0(2), v0(3);1, v1(2), v1(3);1, v2(2), v2(3)];
	   f = [v0(1);v1(1);v2(1)];
       elseif normal == 'y'
	   A = [1, v0(1), v0(3);1, v1(1), v1(3);1, v2(1), v2(3)];
	   f = [v0(2);v1(2);v2(2)];
       elseif normal == 'z'
	   A = [1, v0(1), v0(2);1, v1(1), v1(2);1, v2(1), v2(2)];
	   f = [v0(3);v1(3);v2(3)];
       end

       % solve matrix and get desired coordinate:
       c = A\f;
       est(i) = c(1)+c(2)*newst(i,1)+c(3)*newst(i,2);
   end

   % plot markers into mesh:
   if normal == 'x'
       xplot = est';
       yplot = newst(:,1);
       zplot = newst(:,2);
   elseif normal == 'y'
       xplot = newst(:,1);
       yplot = est';
       zplot = newst(:,2);
   elseif normal == 'z'
       xplot = newst(:,1);
       yplot = newst(:,2);
       zplot = est';
   end

   % clean double receiver positions (due to double counted fault sides)
   receivers_part = [xplot,yplot,zplot];
   receivers = vertcat(receivers,receivers_part);
end
[tmp,uniquecol]=unique(receivers(:,1:2:end),'rows');
receivers=receivers(uniquecol,:);

disp(sprintf(' %d/%d receiver(s) could be located on the fault', size(receivers,1), size(st,1)));
 

if plotmarker == 'y'
    plot3(receivers(:,1),receivers(:,2),receivers(:,3),'r*','MarkerSize',8)
    axis tight
end

%disp('    Receiver coordinates at fault:'); disp(' ')
%disp(receivers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SAVE RECEIVER LOCATIONS TO FILE         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choice = input(['    Save receiver coordinates as ',filename,'_faultreceivers.dat ? (y/n)'], 's');
disp(['    receiver coordinates saved as ',filename0,'_faultreceivers.dat']);
choice='y';

if choice=='y'
    fid_out  = fopen([filename0 ,'_faultreceivers.dat'],'w');
    fprintf(fid_out,'%20.12f %20.12f %20.12f\n',receivers');
    %fprintf(fid_out,'%.20e %.20e %.20e\n',receivers');
    fclose(fid_out);
    disp('    Receiver coordinates saved!');
else
    disp('    Receiver coordinates NOT saved!');
end



disp('    FINISHED');
 

