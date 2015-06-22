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
disp('    %%                 GAMBIT_CHECK_INSPHER                %%')
disp('    %%            TO ANALYSE QUALITY OF 3D_MESH            %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

close all;

filename = input('    Filename(suffix ".neu" assumed):  ','s');

disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Reading data from: \t\t%s',[filename,'.neu']));
disp(' ')
fid   = fopen([filename ,'.neu']);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
param = fscanf(fid,'%i',[6,1]);  NX = param(1);  NT = param(2); 
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
X     = fscanf(fid,'%g',[4,NX]); X = X'; X(:,1) = [];
junk  = fgetl(fid);
junk  = fgetl(fid);
junk  = fgetl(fid);
tetra = fscanf(fid,'%g',[7,NT]); tetra = tetra'; tetra(:,1:3) = [];
junk  = fgetl(fid);
fclose(fid);
disp(sprintf('\t Read \t\t\t\t\t%i vertices',NX));
disp(sprintf('\t Read \t\t\t\t\t%i tetrahedral elements',NT));

disp('-----------------------------------------------------------------------------------')
eval(['load ',filename,'_stencil.dat']);
stencil = cube_005_stencil;
n = input(     'Given number of tetrahedron:   ');
stencil = stencil(n,:);
figure(1)
s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(3,:) = [1,4,3];   s_vert(4,:) = [2,3,4];
for i = 1:length(stencil)
    t = stencil(i);
    tmp = [X(tetra(t,:),1),X(tetra(t,:),2),X(tetra(t,:),3)];
    for j = 1:4
        if i==1
           fill3(tmp(s_vert(j,:),1),tmp(s_vert(j,:),2),tmp(s_vert(j,:),3),'b'); hold on, grid on 
        else
           fill3(tmp(s_vert(j,:),1),tmp(s_vert(j,:),2),tmp(s_vert(j,:),3),'r'); hold on
        end    
    end
    alpha(0.5), axis equal
    pause
end

x_min = min(X(:,1));  x_max = max(X(:,1));
y_min = min(X(:,2));  y_max = max(X(:,2));
z_min = min(X(:,3));  z_max = max(X(:,3));
%plotting tetrahedron with smallest insphere
i = i(1);
s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(3,:) = [1,4,3];   s_vert(4,:) = [2,3,4];
tmp = [X(tetra(i,:),1),X(tetra(i,:),2),X(tetra(i,:),3)];
tx_min = min(tmp(:,1));  tx_max = max(tmp(:,1));
ty_min = min(tmp(:,2));  ty_max = max(tmp(:,2));
tz_min = min(tmp(:,3));  tz_max = max(tmp(:,3));
tMAX = max(abs([tx_min,tx_max,ty_min,ty_max,tz_min,tz_max]));
figure(2)
subplot(121), hold on
for j = 1:4
     fill3(tmp(s_vert(j,:),1),tmp(s_vert(j,:),2),tmp(s_vert(j,:),3),'r');
end
axis equal, grid on, axis([x_min x_max y_min y_max z_min z_max]); 
xlabel('Location in x');ylabel('Location in y');zlabel('Location in z');
subplot(122), hold on
for j = 1:4
     fill3(tmp(s_vert(j,:),1),tmp(s_vert(j,:),2),tmp(s_vert(j,:),3),'r');
end
axis equal, grid on
axis([tx_min-0.1*tMAX tx_max+0.1*tMAX ty_min-0.1*tMAX ty_max+0.1*tMAX tz_min-0.1*tMAX tz_max+0.1*tMAX]);
xlabel('zoomed in x');ylabel('zoomed in y');zlabel('zoomed in z');
