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
% the normal vector to the surface can be any vector
home;
if verLessThan('matlab', '8.5')
    disp('YOUR MATLAB VERSION IS TOO OLD FOR THAT SCRIPT')
    return
end

disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%      PLACE_QUICKLY_FAULTRECEIVERS_NETCDF_ANY_NORMAL %%')
disp('    %%      TO COMPUTE FAULTRECEIVER POSITIONS             %%')
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
filename0    = input('    Filename of mesh(suffix ".nc" is appended)                  :  ','s');
snormal      = input('    Normal coordinates \"nx,ny,nz\"                               :  ','s');
rec_filename = input('    Filename of 3D receiver locations(suffix ".dat" is appended):  ','s');
SufaceId     = input('    surface to be considered (1: Free Surface, 3: Fault)        :  ','s');
plotmesh     = input('    Do you wish to plot the mesh? (y,n)                         :  ','s');
if plotmesh == 'y'
  plotmarker = input('    Do you wish to plot the marker into the mesh? (y,n)         :  ','s');
else
  plotmarker = 'n';
end
SufaceId = str2num(SufaceId); 
if SufaceId==1
   depthdecrement     = input('    depth (in m) to decrement for more accurate outputs (ex 0 or 50)        :  ','s');
   depthdecrement = str2double(depthdecrement);
end

coords = strsplit(snormal,',');
normal(1) = str2double(coords(1));
normal(2) = str2double(coords(2));  
normal(3) = str2double(coords(3));  
normal = normal/norm(normal);
%now determine 2 perpendicular vectors:
randomV = rand(3,1);
ux = cross(randomV,normal);
ux = ux/norm(ux);
uy = cross(ux, normal);

% load receiver stations
%eval(['load ',rec_filename,'.dat']);
%eval(['st = ',rec_filename,';']);
st = load([rec_filename,'.dat']);
if size(st,2)~=3
    disp(sprintf('Error in %s: number of columns not equal 3',rec_filename))
    return
end

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
   [tetraId, faceId] = find(elemBound==SufaceId);
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
   x1 = zeros(size(X,1),2);
   for kkk=1:size(X,1)
      x1(kkk,1) = dot(X(kkk,:),ux); 
      x1(kkk,2) = dot(X(kkk,:),uy);
   end
   TR = triangulation(tri,x1);
   x2 = zeros(size(st,1),2);
   for kkk=1:size(st,1)
      x2(kkk,1) = dot(st(kkk,:),ux); 
      x2(kkk,2) = dot(st(kkk,:),uy);
   end
   [t,barcoords] = pointLocation(TR,x2);
   indicesInTriangulation = find(~isnan(t));
   ninside  = size(indicesInTriangulation,1);
   
      %Print mesh information on screen
   disp(sprintf('Partition %d:\tFound %g fault triangles,\t%d points inside the triangulation',iPartition, NS, ninside));

   t=t(indicesInTriangulation);
   barcoords = barcoords(indicesInTriangulation,:);
   
     
   newst=st(indicesInTriangulation,1:2);
   k=size(t);
   
   for i=1:k
        % load the three vertices for the receiver position
       v0 = X(tri(t(i),1),:);
       v1 = X(tri(t(i),2),:);
       v2 = X(tri(t(i),3),:);
       receivers_part = barcoords(i,1)*v0 +  barcoords(i,2)*v1+barcoords(i,3)*v2;
       receivers = vertcat(receivers,receivers_part);       
   end
end
if SufaceId==1
    disp(sprintf('!!!removing %f in normal direction!!!!!!', depthdecrement));
    for i =1:size(receivers,1)
        receivers(i,:) = receivers(i,:)-depthdecrement*normal;
    end
end
disp(sprintf(' %d/%d receiver(s) could be located on the fault', size(receivers,1), size(st,1)));
 

if plotmarker == 'y'
    plot3(receivers(:,1),receivers(:,2),receivers(:,3),'r*','MarkerSize',8)
    axis tight
    view([150 30])
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

%Now writing _sameorder file, in which the initial order of the receiver
%(input file) is kept
disp(['    receiver coordinates saved as ',filename0,'_sameorder.dat']);
%for i=1:2
fid_out  = fopen([filename0 ,'_sameorder.dat'],'w');
for i=1:size(st,1)
    indexes2 = find((abs(receivers(:,1)-st(i,1))<1e-3) & (abs(receivers(:,2)-st(i,2))<1e-3));
    fprintf(fid_out,'%20.12f %20.12f %20.12f\n',st(i,1:2), receivers(indexes2,3));
end
disp('    Receiver coordinates saved!');
fclose(fid_out);
