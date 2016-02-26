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


% Thomas Ulrich 02.2016
home;

disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%      MIRROR_MESH : MAKE SYMMETRICAL A NETCDF        %%')                                %%')
disp('    %%      UNPARTITIONNED HALF-MESH                       %%')
disp('    %%      OUTPUT MESH IN GAMBIT .neu FORMAT              %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp('Make symmetrical an half-mesh.')
disp('Only 3 planes of symmetry are available x=0, y=0, z=0.')
disp('The file was intended to read and write netcdf mesh,')
disp('but because of the complex structures of netcdf mesh,')
disp('it finally writes the mesh in neu format.')
disp('A limitation is that the input netcdf mesh have to be unpartitionned.')
disp('Another limitation is that all elements are written inside the same group.')
disp(' ')
clear, close all;


% get input
filename0  = input('   Filename of mesh(suffix ".nc" is appended)                  :  ','s');
normal     = input('    Normal axis to plane of symmetry (x,y,z)                   :  ','s');
%filename0 = 'halfmodel_planar.1'
%normal = 'y'

%precision for detecting nodes on the plane of symmetry
eps = 1e-5;

% input check
switch normal
    case 'x'
        iaxis=1;
    case 'y'
        iaxis=2;
    case 'z'
        iaxis=3;
    otherwise
        disp('Wrong input for normal axis!')
        return
end


filename_new = '/export/data/ulrich/sym_mesh.neu';
disp(sprintf('output will be written in %s', filename_new))

filename = sprintf('%s.nc',filename0);

element_size  = ncread(filename,'element_size');
nPartition = size(element_size,1);
element_vertices  = ncread(filename,'element_vertices');
element_boundaries  = ncread(filename,'element_boundaries');
vertex_coordinates  = ncread(filename,'vertex_coordinates');
vertex_size = ncread(filename,'vertex_size');
element_group = ncread(filename,'element_group');

vertices  = max(vertex_size);
elements = max(element_size);

%initialize New array
vertex_size_new = zeros(nPartition,1);
vertex_coordinates_new = zeros(3,2*vertices,nPartition);
element_size_new = zeros(nPartition,1);
element_vertices_new = zeros(4,2*elements,nPartition);
element_group_new = zeros(2*elements,nPartition);

for iPartition=1:nPartition

   %read vertices
   X = vertex_coordinates(:,1:vertex_size(iPartition),iPartition);
   
   %initialize New array
   Vertex_Old2New = zeros(vertex_size(iPartition),1);
   Elem_Old2New = zeros(element_size(iPartition),1);
   
   vertex_coordinates_new(:,1:vertex_size(iPartition)) = X;
   currentNewVertex = vertex_size(iPartition);

   %create new_vertices
   for ivertex=1:vertex_size(iPartition)
       if abs(X(iaxis,ivertex))<eps
           Vertex_Old2New(ivertex)=ivertex;
       else
           currentNewVertex = currentNewVertex+1;          
           Vertex_Old2New(ivertex)=currentNewVertex;
           vertex_coordinates_new(1:3,currentNewVertex,iPartition) = X(1:3,ivertex);
           vertex_coordinates_new(iaxis,currentNewVertex,iPartition) = -X(iaxis,ivertex);
       end
   end
   vertex_size_new(iPartition) = currentNewVertex;
   
   element_vertices_new(:,1:element_size(iPartition),iPartition) = element_vertices(:,1:element_size(iPartition),iPartition);
   %create New elements
   currentNewElement = element_size(iPartition);
   for ielem=1:element_size(iPartition)
       %read elements (+1 because numeration starts at 0 in netcdf)       
       tetra = element_vertices(:,ielem,iPartition)+1;
       NewElement =0;
       for k=1:4
           if Vertex_Old2New(tetra(k)) ~= tetra(k)
               NewElement =1;
               break
           end
       end
       Elem_Old2New(currentNewElement) = ielem;
       if NewElement
          currentNewElement = currentNewElement+1;
          Elem_Old2New(ielem) = currentNewElement;
          for k=1:4
             element_vertices_new(k,currentNewElement,iPartition) = Vertex_Old2New(tetra(k))-1;
          end
          %switch nodes 3 and 4 for getting a positive volume
          a = element_vertices_new(3,currentNewElement,iPartition);
          element_vertices_new(3,currentNewElement,iPartition) = element_vertices_new(4,currentNewElement,iPartition);
          element_vertices_new(4,currentNewElement,iPartition) =a;
       else
           Elem_Old2New(ielem) = ielem;
       end 
       
   end
   element_size_new(iPartition) = currentNewElement;

   %create New groups
   element_group_new(1:element_size(iPartition),iPartition) = element_group_new(1:element_size(iPartition),iPartition);
   for ielem=1:element_size(iPartition)
       if Elem_Old2New(ielem) ~= ielem
           element_group_new(Elem_Old2New(ielem),iPartition) = element_group(ielem,iPartition);
       end
   end
end
disp('done reading mesh');

vertices_new = max(vertex_size_new);
elements_new = max(element_size_new);
vertex_coordinates_new = vertex_coordinates_new(1:3,1:vertices_new,1:nPartition);
element_vertices_new = element_vertices_new(1:4,1:elements_new,1:nPartition);

if nPartition~=1
    disp('This script currently only works with 1 partition!')
    return
end
ngroups = max(element_group_new);

fid_out = fopen('/export/data/ulrich/sym_mesh.neu', 'w');
fprintf(fid_out,'%s','        CONTROL INFO 1.2.1');
fprintf(fid_out,'\n%s','** GAMBIT NEUTRAL FILE');
fprintf(fid_out,'\n%s','Simmetrix mesh in GAMBIT neutral file format');
fprintf(fid_out,'\n%s','PROGRAM:               Gambit     1.2.1');
fprintf(fid_out,'\n%s','Mon Feb 15 11:39:00 2016');
fprintf(fid_out,'\n%s','     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL');
%     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
%       606      2538         1         3         3         3

param = zeros(6,1);
param(1)=vertices_new;
param(2) = elements_new;     
param(3) = ngroups;
param(4) = 3;
param(5) = 3;
fprintf(fid_out,'\n%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f',param);

fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','   NODAL COORDINATES 1.2.1');
for i=1:vertices_new
    fprintf(fid_out,'\n%10.0f%20.10e%20.10e%20.10e',i,vertex_coordinates_new(:,i,1));
end
disp('done writing nodes');

fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','      ELEMENTS/CELLS 1.2.1');

for i=1:elements_new
    fprintf(fid_out,'\n%8.0f%3.0f%3.0f%9.0f%8.0f%8.0f%8.0f',i, 6,4, element_vertices_new(1:4,i,1)+1);
end
disp('done writing elements');

fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','       ELEMENT GROUP 2.0.0');
fprintf(fid_out,'\n%s','GROUP:          1 ELEMENTS:');
fprintf(fid_out,'%11.0f%',elements_new);
fprintf(fid_out,'%s',' MATERIAL:          2 NFLAGS:          1');
fprintf(fid_out,'\n%s','                           fluid');
fprintf(fid_out,'\n%s','       0');
fprintf(fid_out,'\n%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f',(1:1:elements_new));
fprintf(fid_out,'\n%s\n','ENDOFSECTION');
disp('done writing groups');


%GAMBIT convention
s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(3,:) = [2,3,4];   s_vert(4,:) = [1,4,3];


%read elements boudary
elemBound = element_boundaries(:,1:element_size(iPartition),iPartition);
elemBound = elemBound';

%find Boundary elements
[tetraId, faceId] = find(elemBound==1);
%switch 3 and 4 (ICEM vs Gambit convention)
ind3 = find(faceId==3);
ind4 = find(faceId==4);
faceId(ind3)=4;
faceId(ind4)=3;

fprintf(fid_out,' BOUNDARY CONDITIONS 2.3.16\n');
nbound=size(tetraId,1);
fprintf(fid_out,'%s%8.0f%8.0f%8.0f%8.0f\n','                             101',1,2*nbound,0,6);
for i = 1:nbound
    fprintf(fid_out,'%10.0f%5.0f%5.0f\n', tetraId(i),6,faceId(i));
end

%switch 1 and 2 (inversion indices due to symmetry)
ind1 = find(faceId==1);
ind2 = find(faceId==2);
faceId(ind1)=2;
faceId(ind2)=1;
for i = 1:nbound
    currentvertices = element_vertices_new(s_vert(faceId(i),:), Elem_Old2New(tetraId(i)), 1);
    if max(currentvertices)>vertices
        %then one vertices has been created
        fprintf(fid_out,'%10.0f%5.0f%5.0f\n', Elem_Old2New(tetraId(i)),6,faceId(i));
    end
end
fprintf(fid_out,'%s\n','ENDOFSECTION');

%find Boundary elements
[tetraId, faceId] = find(elemBound==3);
%switch 3 and 4 (ICEM vs Gambit convention)
ind3 = find(faceId==3);
ind4 = find(faceId==4);
faceId(ind3)=4;
faceId(ind4)=3;

fprintf(fid_out,' BOUNDARY CONDITIONS 2.3.16\n');
nbound=size(tetraId,1);
fprintf(fid_out,'%s%8.0f%8.0f%8.0f%8.0f\n','                             103',1,nbound,0,6);
for i = 1:nbound
    fprintf(fid_out,'%10.0f%5.0f%5.0f\n', tetraId(i),6,faceId(i));
end

%switch 1 and 2 (inversion indices due to symmetry)
ind1 = find(faceId==1);
ind2 = find(faceId==2);
faceId(ind1)=2;
faceId(ind2)=1;
for i = 1:nbound
    currentvertices = element_vertices_new(s_vert(faceId(i),:), Elem_Old2New(tetraId(i)), 1);
    if max(currentvertices)>vertices
        %then one vertices has been created
        fprintf(fid_out,'%10.0f%5.0f%5.0f\n', Elem_Old2New(tetraId(i)),6,faceId(i));
    end
end
fprintf(fid_out,'%s\n','ENDOFSECTION');

[tetraId, faceId] = find(elemBound==5);
ind3 = find(faceId==3);
ind4 = find(faceId==4);
faceId(ind3)=4;
faceId(ind4)=3;
fprintf(fid_out,' BOUNDARY CONDITIONS 2.3.16\n');
nbound=size(tetraId,1);
fprintf(fid_out,'%s%8.0f%8.0f%8.0f%8.0f\n','                             105',1,2*nbound,0,6);
for i = 1:nbound
    fprintf(fid_out,'%10.0f%5.0f%5.0f\n', tetraId(i),6,faceId(i));
end
%switch 1 and 2 (inversion indices due to symmetry)
ind1 = find(faceId==1);
ind2 = find(faceId==2);
faceId(ind1)=2;
faceId(ind2)=1;
for i = 1:nbound
    currentvertices = element_vertices_new(s_vert(faceId(i),:), Elem_Old2New(tetraId(i)), 1);
    if max(currentvertices)>vertices
        %then one vertices has been created
        fprintf(fid_out,'%10.0f%5.0f%5.0f\n', Elem_Old2New(tetraId(i)),6,faceId(i));
    end
end
fprintf(fid_out,'%s\n','ENDOFSECTION');

fclose(fid_out);

disp('    FINISHED');
 

