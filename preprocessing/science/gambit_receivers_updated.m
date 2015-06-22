%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
% @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
%
% @section LICENSE
% Copyright (c) 2005-2013, SeisSol Group
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
disp('    %%                  GAMBIT_RECEIVERS                   %%')
disp('    %%           TO COMPUTE STATION ELEVATIONS             %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(' Give Gambit mesh-file and station coordinates in order to compute')
disp(' the exact elevation of the stations on the particular Gambit mesh,')
disp(' which approximates the topography with piecewise linear polynomials')
disp(' in the surface triangulation!')
disp(' '),disp(' ')
clear, close all;

filename = input('    Filename of mesh(suffix ".neu" is appended):  ','s');
fid   = fopen([filename ,'.neu']);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
param = fscanf(fid,'%i',[6,1]);  NX = param(1);  NT = param(2); NG = param(3); NB = param(4);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
X     = fscanf(fid,'%g',[4,NX]); X = X'; X(:,1) = []; X(:,3)=X(:,3)+rand(size(X(:,3)))*0.0000001;
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
tetra = fscanf(fid,'%g',[7,NT]); tetra = tetra'; tetra(:,1:3) = [];
junk  = fgetl(fid);
for i = 1:NG
    junk  = fgetl(fid);
    junk  = fgetl(fid);
    junk  = fgetl(fid);
    n     = str2num(junk(30:38));
    junk  = fgetl(fid);
    junk  = fgetl(fid);
    ju  = fscanf(fid,'%g',[10,ceil(n/10)]);
end
junk  = fgetl(fid); 
if(length(junk)==0)
   junk = fgetl(fid);
end
for i = 1:NB
    junk = fgetl(fid);
    param = fscanf(fid,'%i',[5,1]);  NS = param(3);
    if (param(1)~=101)
        junk  = fgetl(fid);
        TS    = fscanf(fid,'%g',[3,NS]); LocalFace = TS(3,:)'; TS = TS(1,:)';
        junk  = fgetl(fid); junk = fgetl(fid);
    else    
        junk  = fgetl(fid);
        TS    = fscanf(fid,'%g',[3,NS]); LocalFace = TS(3,:)'; TS = TS(1,:)';
        junk  = fgetl(fid); break
    end
end
fclose(fid);

rec_filename = input('    Filename of receiver locations(suffix ".dat" is appended):  ','s');
eval(['load ',rec_filename,'.dat']);
eval(['st = ',rec_filename,';']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PLOT SURFACE TRIANGULATION         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(3,:) = [2,3,4];   s_vert(4,:) = [1,4,3];

Y = []; tri = [];
j = 0;

Y   = X(tetra(TS(1),s_vert(LocalFace(1),:))',:);
tri = [1,2,3];
for i = 2:NS
    tmp        = X(tetra(TS(i),s_vert(LocalFace(i),:))',:);
    mem        = ismember(Y,tmp);
    mloc       = ismember(tmp,Y);
    a          = find(sum(mem,2)==3);
    b          = find(sum(mloc,2)==3);
    if length(b)>0
        tmp(b,:) = [];
    end
    Y   = [Y;tmp];
     
    trinew = [length(Y):-1:length(Y)-2];
    trinew = [a',trinew];
    tri = [tri;trinew(1:3)];
%    fill3(tmp(:,1),tmp(:,2),tmp(:,3),'r');
    if(mod(i,round(NS/10))==0)
        j = j+1;
        disp(sprintf('     %g percent done!',j*10))
    end
end

%Print mesh information on screen
disp(sprintf('\n\tFound %g surface triangles',NS));
disp(sprintf('\tFound %g surface stations with xyz-coordinates:',size(st,1)));
figure; hold on
trisurf(tri,Y(:,1),Y(:,2),Y(:,3),Y(:,3)); axis equal,xlabel('x'),ylabel('y'),zlabel('z')
%save surface triangulation for use by gambit_receivers_fast
X = Y;
disp('    Surface triangulation with vertices saved!'); disp(' ')
eval(['save ',filename,'_tri.mat tri'])
eval(['save ',filename,'_vert.mat X'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          COMPUTE STATION ELEVATION          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr=triangulation(tri,Y(:,1),Y(:,2));
Cedges = edges(tr);
dt = delaunayTriangulation(Y(:,1),Y(:,2), Cedges);
dtTotrIndex = nan(size(dt,1),1);
ic = incenter(tr);
dtid = pointLocation(dt,ic);
dtTotrIndex(dtid) = 1:size(tr,1);
t = pointLocation(dt,st(:,1),st(:,2));
if isfinite(t)
   t = dtTotrIndex(t);
end
%t = tsearch(Y(:,1),Y(:,2),tri,st(:,1),st(:,2)); !OBSOLETE
for i = 1:size(st,1)
    %t(i) = tsearch(Y(:,1),Y(:,2),tri,st(i,1),st(i,2)); !OBSOLETE
    v0 = Y(tri(t(i),1),:);
    v1 = Y(tri(t(i),2),:);  
    v2 = Y(tri(t(i),3),:); 
    A = [1, v0(1), v0(2);1, v1(1), v1(2);1, v2(1), v2(2)];
    f = [v0(3);v1(3);v2(3)];
    c = A\f;
    est(i) = c(1)+c(2)*st(i,1)+c(3)*st(i,2); 
end
depth = input('    Specify depth to burry the receivers under the surface (in mesh units):  ');
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
