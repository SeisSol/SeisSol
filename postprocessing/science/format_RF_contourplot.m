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
% script for SCEC TPV rupture front contour plot formatting
% assume planar fault plane in the xz plane

clear all

tpv = input('Specify TPV case:  ','s');

disp('Code revision 782???');

% load mat file from read_RF.m script
disp('Load data...')
data=loadmatfile('RF.mat');

%% prepare data

% kick out y component
tmp(1,:) = data(1,:);
tmp(2:3,:) = data(3:4,:);
clear data

% make depth positive
tmp(2,:)=abs(tmp(2,:)); 

disp('Reduce data by a factor of 4 to match the max allowed file size!')
RF=tmp(:,1:4:end);
clear tmp;


% ATTENTION: right-lateral rupture has strike direction -x and has to be
% changed
disp('Prepare data for right-lateral case - switch sign of x coordinate!')
RF(1,:)=-1.*RF(1,:);

[a b] = size(RF);

%% Save file

str3 = date;

% open file
disp('Save file...')
fid=fopen('cplot.dat','w');

fprintf(fid,['# problem=TPV',tpv,'\n']);
fprintf(fid,'# author=Christian Pelties\n');
fprintf(fid,['# date=',str3,'\n']);
fprintf(fid,'# code=SeisSol (ADER-DG)\n');
fprintf(fid,'# code_version=revision 782\n');
fprintf(fid,'# element_size=200 m\n');
fprintf(fid,'# order of approximation in space and time= O4\n');
fprintf(fid,'# Column #1 = horizontal coordinate, distance along strike (m)\n');
fprintf(fid,'# Column #2 = vertical coordinate, distance down-dip (m)\n');
fprintf(fid,'# Column #3 = rupture time (s)\n');
fprintf(fid,'#\n');
fprintf(fid,'# The line below lists the names of the data fields.\n');
fprintf(fid,'# It indicates that the first column contains the horizontal\n');
fprintf(fid,'# coordinate (j), the second column contains the vertical\n');
fprintf(fid,'# coordinate (k), and the third column contains the time (t).\n');
fprintf(fid,'j k t\n');
fprintf(fid,'#\n');
fprintf(fid,'# Here is the rupture history\n');
for j = 1:b
fprintf(fid,'%15.7e %15.7e %15.7e \n', ...
            RF(1,j),RF(2,j),RF(3,j));
end;

% close file
fclose(fid);
disp('Done!')
