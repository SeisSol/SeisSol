%%
% @file
% This file is part of SeisSol.
% SPDX-License-Identifier: BSD-3-Clause
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
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

home;
disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%                PLACE_FAULTRECEIVERS                 %%')
disp('    %%         TO COMPUTE FAULTRECEIVER POSITIONS          %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(' Give Gambit mesh-file and station coordinates in order to compute')
disp(' the exact positions of the stations on the particular Gambit mesh,')
disp(' which approximates the fault with piecewise linear polynomials')
disp(' in the fault-surface triangulation!')
disp(' '),disp(' ')
clear, close all;

format long

% get input
filename     = input('    Filename of mesh(suffix ".neu" is appended)                 :  ','s');
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

% read neu file
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

% read nodes
X     = fscanf(fid,'%g',[4,NX]);X = X'; 
X(:,1) = []; % delete indices
X(:,3)=X(:,3)+rand(size(X(:,3)))*0.0000001; % random field needed to avoid NaNs. it supports matlab by finding the closest neighbor.
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
% read connectivity information
tetra = fscanf(fid,'%g',[7,NT]); tetra = tetra'; tetra(:,1:3) = [];
junk  = fgetl(fid);

% read groups
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

% read NBSETS
for i = 1:NB
    junk = fgetl(fid);
    param = fscanf(fid,'%i',[5,1]);  NS = param(3);
    % search for fault elements
    if (param(1)~=103)
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
disp(sprintf('\n\tFound %g fault triangles',NS));
disp(sprintf('\tFound %g faultreceiver with xyz-coordinates:',size(st,1)));

if plotmesh == 'y'
    figure; hold on
    trimesh(tri,Y(:,1),Y(:,2),Y(:,3),Y(:,3)); axis equal,xlabel('x'),ylabel('y'),zlabel('z')
end

% saving of data not necessary
%  %save surface triangulation for use by gambit_receivers_fast
%  X = Y;
%  disp('    Fault-surface triangulation with vertices saved!'); disp(' ')
%  eval(['save ',filename,'_tri.mat tri'])
%  eval(['save ',filename,'_vert.mat X'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          COMPUTE STATION ELEVATION          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note, fault sides appear twice and therefore two faces will be found that
% will be cleaned out later. This way, a fault station will be projected
% also to branches lying behind the main fault.

k=0;
for i = 1:size(st,1)
    for j = 1:size(tri)
        % get correct orientation
        if normal == 'x'
            x(1) = Y(tri(j,1),2);
            x(2) = Y(tri(j,2),2);
            x(3) = Y(tri(j,3),2);
            y(1) = Y(tri(j,1),3);
            y(2) = Y(tri(j,2),3);
            y(3) = Y(tri(j,3),3);        
        elseif normal == 'y'
            x(1) = Y(tri(j,1),1);
            x(2) = Y(tri(j,2),1);
            x(3) = Y(tri(j,3),1);
            y(1) = Y(tri(j,1),3);
            y(2) = Y(tri(j,2),3);
            y(3) = Y(tri(j,3),3);        
        elseif normal == 'z'
            x(1) = Y(tri(j,1),1);
            x(2) = Y(tri(j,2),1);
            x(3) = Y(tri(j,3),1);
            y(1) = Y(tri(j,1),2);
            y(2) = Y(tri(j,2),2);
            y(3) = Y(tri(j,3),2);        
        end

        inside = XYinElement(x,y,st(i,1),st(i,2));
        if inside == 1
            k=k+1;
            t(k) = j;
            newst(k,1)=st(i,1);
            newst(k,2)=st(i,2);
        end
    end
end

est(1:k) = 0.0;

% linear interpolation
for i=1:k  % add new stations generated by branch e.g.
    
    % load the three vertices for the receiver position
    v0 = Y(tri(t(i),1),:);
    v1 = Y(tri(t(i),2),:);  
    v2 = Y(tri(t(i),3),:);
    
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
receivers = [xplot,yplot,zplot];
receivers = unique(receivers,'rows');

if plotmarker == 'y'
    plot3(receivers(:,1),receivers(:,2),receivers(:,3),'r*','MarkerSize',8)
    axis tight
end

%disp('    Receiver coordinates at fault:'); disp(' ')
%disp(receivers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SAVE RECEIVER LOCATIONS TO FILE         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = input(['    Save receiver coordinates as ',filename,'_faultreceivers.dat ? (y/n)'], 's');
if choice=='y'
    fid_out  = fopen([filename ,'_faultreceivers.dat'],'w');
    fprintf(fid_out,'%20.12f%20.12f%20.12f\n',receivers');
    fclose(fid_out);
    disp('    Receiver coordinates saved!');
else
    disp('    Receiver coordinates NOT saved!');
end

disp('    FINISHED');
