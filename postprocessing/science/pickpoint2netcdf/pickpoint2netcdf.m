%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
% @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
%
% @section LICENSE
% Copyright (c) 2006-2013, SeisSol Group
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
% Takes a set of pick point files (in TecPlot format) and produces a netCDF file respecting the COARDS conventions.

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%   Takes a set of pick point files (in TecPlot format) and    %%')
disp('     %%   produces a netCDF file respecting the COARDS               %%')
disp('     %%   conventions.                                               %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
dirname    = input('      Sepcify the directory containing the pickpoint files ', 's');
firstPoint = input('      Specify the index of the first interesting pickpoint (0)');
nPoints    = input('      Specify the number of interesting pickpoints         ');
filename   = input('      Sepcify the output filename ', 's');
writeonsampleevery = 1;
%if alldir=3, x and y displacements are also written
alldir = false;     

% Get all files in the directory
files = dir(dirname);
files = {files(:).name};

% Find the pick point files
file_indices = ~cellfun(@isempty, strfind(files, 'receiver'));
files = files(arrayfun(@(x) x, file_indices));

% Parse files
location = NaN(2, nPoints);
for num = 1:nPoints
     
    in_file = [dirname,'/',files{num+firstPoint}];
    fid   = fopen(in_file);
    
    % title
    fgetl(fid);
    
    % variable names -> find time and z
    vars = regexp(fgetl(fid), '"\s*,\s*"', 'split');
    vars(1) = regexprep(vars(1), '\s*VARIABLES\s*=\s*"', '', 'once');
    vars(end) = regexprep(vars(end), '"\s*', '', 'once');
    timeIndex = find(strcmp(vars, 'Time'));
    zIndex = find(strcmp(vars, 'w'));
    if alldir
        xIndex = find(strcmp(vars, 'u'));
        yIndex = find(strcmp(vars, 'v'));
    end
    
    % Read the coordinate of the pickpoint
    fscanf(fid,'%c',[6,1]);
    x1    = fscanf(fid,'%g',[1,1]);
    fscanf(fid,'%c',[6,1]);
    x2    = fscanf(fid,'%g',[1,1]);
    fscanf(fid,'%c',[6,1]);
    x3    = fscanf(fid,'%g',[1,1]);
    
    location(:, num) = [x2, x1];
    
    % Read the actual data
    fprintf('X: %g   Y: %g  Z: %g\n',x1,x2,x3);
    data = fscanf(fid,'%g',[size(vars, 2),inf]); data = data';
    
    fclose(fid);
    
    %determine and delete double samples from a possible restart
    data(:,1) = (round(data(:,1)*1000000))/1000000;
    [tmp,k,kk] = unique(data(:,1));
    data  = data(k,:);
    
    fprintf('Samples in seismogram nr. %i\t:   %i\n', num,size(data,1));
    
    time = data(:,timeIndex)';
    if alldir
        valuesLinear(:,num,1) = data(:,xIndex);
        valuesLinear(:,num,2) = data(:,yIndex);
        valuesLinear(:,num,3) = data(:,zIndex);
    else
        valuesLinear(:,num) = data(:,zIndex);
    end
end
clearvars data kk tmp

% Generate ycoords, xcoords and values from location and valuesLinear
ycoords = NaN(size(unique(location(1,:))));
xcoords = NaN(size(unique(location(2,:))));
ynext = 1;
xnext = 1;
values = NaN(length(time), length(ycoords), length(xcoords));

%First loop for setting xcoords and ycoords
for num = 1:size(location, 2)    
    % Find the indices of the x and y coordinate
    % Generate new indices of not found
    y = find(ycoords == location(1,num));
    if isempty(y)
        ycoords(ynext) = location(1,num);
        y = ynext;
        ynext = ynext + 1;
    end
    x = find(xcoords == location(2,num));
    if isempty(x)
        xcoords(xnext) = location(2,num);
        x = xnext;
        xnext = xnext + 1;
    end    
end
if ~issorted(ycoords) || ~issorted(xcoords)
    disp('Warning: Pickpoint files are not sorted by x and y coordinates')
    disp('attempting to Reorder')
end

xcoords = sort(xcoords);
ycoords = sort(ycoords);
if  size(xcoords,2)*size(ycoords,2) ~= size(location, 2)
    disp('size(xcoords,2)*size(ycoords,2) ~= number of receivers')
    return
end


for num = 1:size(location, 2)
    
    % Find the indices of the x and y coordinate
    y = find(ycoords == location(1,num));
    x = find(xcoords == location(2,num));
    if alldir
       values(:,y,x,1:3) = valuesLinear(:,num,1:3);
    else
       values(:,y,x) = valuesLinear(:,num);
    end
end




clearvars valuesLinear



disp(' '),  disp('     Finished conversion!'), disp(' ')

% Write netCDF

warning('OFF','MATLAB:DELETE:FileNotFound')
delete(filename)

nccreate(filename,'x','Dimensions',{'x',length(xcoords)},'Datatype','single')
ncwrite(filename,'x',xcoords)
nccreate(filename,'y','Dimensions',{'y',length(ycoords)},'Datatype','single')
ncwrite(filename,'y',ycoords)
nccreate(filename,'time','Dimensions',{'time',floor(length(time)/writeonsampleevery)},'Datatype','single')
ncwriteatt(filename,'time','units','seconds since initial rapture')
timewritten = cat(2,time(1), time(1+writeonsampleevery:writeonsampleevery:writeonsampleevery*floor(length(time)/writeonsampleevery)));
ncwrite(filename,'time',timewritten)

ncwriteatt(filename,'/','Conventions','COARDS')
ncwriteatt(filename,'/','Created','Using SeisSol and pickpoint2netcdf converter')

if alldir
    nccreate(filename,'dx','Dimensions',{'x', 'y', 'time'},'Datatype','single')
    nccreate(filename,'dy','Dimensions',{'x', 'y', 'time'},'Datatype','single')
    nccreate(filename,'dz','Dimensions',{'x', 'y', 'time'},'Datatype','single')
    displOld = zeros(length(ycoords), length(xcoords),3);
    ncwrite(filename,'dx',permute(displOld(:,:,1),[2 1]),[1 1 1])
    ncwrite(filename,'dy',permute(displOld(:,:,2),[2 1]),[1 1 1])
    ncwrite(filename,'dz',permute(displOld(:,:,3),[2 1]),[1 1 1])
else        
    nccreate(filename,'dz','Dimensions',{'x', 'y', 'time'},'Datatype','single')
    displOld = zeros(length(ycoords), length(xcoords));
    ncwrite(filename,'dz',permute(displOld,[2 1]),[1 1 1])
end


for t = 2:length(time)
    % Integrating values
    if alldir
        displ = reshape(values(t-1,:,:,:)+values(t,:,:,:),[length(ycoords) length(xcoords) 3]) * 0.5 * (time(2)-time(1)) + displOld; 
    else
        displ = reshape(values(t-1,:,:)+values(t,:,:),[length(ycoords) length(xcoords)]) * 0.5 * (time(2)-time(1)) + displOld;   
    end
    if mod(t,writeonsampleevery)==0
        fprintf('Writing timestep:   %i\n', t);
        if alldir
           ncwrite(filename,'dx',permute(displ(:,:,1),[2 1]),[1 1 t/writeonsampleevery])
           ncwrite(filename,'dy',permute(displ(:,:,2),[2 1]),[1 1 t/writeonsampleevery])
           ncwrite(filename,'dz',permute(displ(:,:,3),[2 1]),[1 1 t/writeonsampleevery])
        else
           ncwrite(filename,'dz',permute(displ,[2 1]),[1 1 t/writeonsampleevery])    
        end
        

    end
    displOld = displ;
end
