
%%
% @file
% This file is part of SeisSol.
%
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
%
% @section LICENSE
% Copyright (c) 2012, SeisSol Group
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
% Reads in -RF- files.

disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  read_RF reads in the -RF- files that include the            %%')
disp('     %%  rupture front evolution information.                        %%')
disp('     %%  Script expects that all RF files are stored within one      %%')
disp('     %%  folder without any RF from another simulation.              %%')
disp('     %%  Explanation in the struct.                                  %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
workingdir   = input('     workingdir in which are the RF files     :   ','s');
if strcmp(workingdir,'')
    workingdir='.'
end

savedata     = input('     Save data to file? (1=yes, 0=no)         :   ');
if (savedata == 1)
  format     = input('     Data format? (1=.dat, 2=.mat)            :   ');
end;
ploton       = input('     Show plots? (1=yes, 0=no)                :   ');
if (ploton == 1)
  dt_sample  = input('     Contour plot lines every... [in sec]     :   ');
end;
ds_data      = input('     Includes dyn.shear stress arrival? (1=yes, 0=no)         :   ');

disp(' '),  disp('     Read data...' )

% read in data
Ssearch = sprintf('%s/%s',workingdir,'*-RF-*');
liste=dir(Ssearch);

if (ds_data==1)
    ndat = 5; %x,y,z,RF,DS
else
    ndat = 4; %x,y,z,RF
end

files = {liste.name};
istart=1;
for k=1:numel(files)
    % read files
    sfilename = sprintf('%s/%s',workingdir,files{k});
    fid = fopen(sfilename,'r');

    maxnumberlines = fscanf(fid,'%g',[1,1]);
    data_tmp(:,:) = fscanf(fid,'%g',[ndat,maxnumberlines]);
    fclose(fid);
    [junk len] = size(data_tmp);
    iend=(istart-1)+len;
    data(:,istart:iend)=data_tmp;
    istart=iend+1;
    clear data_tmp
end 

% save data
if (savedata == 1)
    if (format == 1)
        % .dat file
        RF=data';
        save RF.dat -ASCII RF;
        clear RF;
        disp(' '),  disp('Data stored unstructured in RF.dat' )
    elseif (format == 2)
        % matlab file
        RF=data;
        save RF.mat -MAT RF;
        clear RF;
        disp(' '),  disp('Data stored unstructured in RF.mat' )
    end
elseif (savedata == 0)
    disp('No data is stored!')
end;

% plot data
if (ploton == 1)
    disp(' '),  disp('     Prepare plot...' )

    % convert data (at the moment for xz plane only!)
    x=data(1,:);
    z=data(3,:);
    t=data(4,:);
    if (ds_data==1)
        ds_t=data(5,:);
    end

    % spatial sampling of data
    nx = 300;
    nz = 150;

    % create strucutred grid
    xlin=linspace(min(x),max(x),nx);
    zlin=linspace(min(z),max(z),nz);
    [X,Z]=meshgrid(xlin,zlin);

    % project unstructured data on structured grid
    T=griddata(x,z,t,X,Z,'cubic');

    % plot in 3D
    mesh(X,Z,T); % interpolated
    axis tight; hold on
    plot3(x,z,t,'.','MarkerSize',15)

    % contour plot
    figure;
    V=0:dt_sample:max(t);
    [c,h]= contour(X,Z,T,V);
    hold on
    axis equal;
    h.LineWidth = 2;
    clabel(c,h,'fontsize',12)

    if (ds_data==1)
        %DS: project unstructured data on structured grid
        index = find(ds_t>0); %plot only values greater than 0
        T_prime=griddata(x(index),z(index),ds_t(index),X,Z,'cubic');

        % contour plot
        V=0:dt_sample:6; %max(ds_t);
        [c1,h1]= contour(X,Z,T_prime,V);
        axis equal;
        h1.LineWidth = 2;
        h1.LineStyle='--';
        clabel(c1,h1,'Color','r','fontsize',12)
    end

elseif (ploton == 0)
    disp('No data is displayed!')
end;

disp('End program!')
