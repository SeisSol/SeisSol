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
firstPoint = input('      Specify the index of the first interesting pickpoint ');
nPoints    = input('      Specify the number of interesting pickpoints         ');
filename   = input('      Sepcify the output filename ', 's');
writeonsampleevery = 1;

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
    valuesLinear(:,num) = data(:,zIndex);
      
end
clearvars data kk tmp

% Generate ycoords, xcoords and values from location and valuesLinear
ycoords = NaN(size(unique(location(1,:))));
xcoords = NaN(size(unique(location(2,:))));
ynext = 1;
xnext = 1;
values = NaN(length(time), length(ycoords), length(xcoords));
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
    
    values(:,y,x) = valuesLinear(:,num);
    
end

clearvars valuesLinear

if ~issorted(ycoords) || ~issorted(xcoords)
    disp('Error: Pickpoint files are not sorted by x and y coordinates')
    disp('Reordering is not supported')
    return
end

disp(' '),  disp('     Finished conversion!'), disp(' ')

% Write netCDF

warning('OFF','MATLAB:DELETE:FileNotFound')
delete(filename)

nccreate(filename,'x','Dimensions',{'x',length(xcoords)},'Datatype','single')
ncwrite(filename,'x',xcoords)
nccreate(filename,'y','Dimensions',{'y',length(ycoords)},'Datatype','single')
ncwrite(filename,'y',ycoords)
nccreate(filename,'time','Dimensions',{'time',length(time/writeonsampleevery)},'Datatype','single')
ncwriteatt(filename,'time','units','seconds since initial rapture')
ncwrite(filename,'time',time)

ncwriteatt(filename,'/','Conventions','COARDS')
ncwriteatt(filename,'/','Created','Using SeisSol and pickpoint2netcdf converter')

nccreate(filename,'z','Dimensions',{'x', 'y', 'time'},'Datatype','single')

displOld = zeros(length(ycoords), length(xcoords));
ncwrite(filename,'z',permute(displOld,[2 1]),[1 1 1])
for t = 2:length(time)
    fprintf('Writing timestep:   %i\n', t);
    
    % Integrating values
    displ = reshape(values(t-1,:,:)+values(t,:,:),[length(ycoords) length(xcoords)]) * 0.5 * (time(2)-time(1)) + displOld;
    if mod(t,writeonsampleevery)==0
       ncwrite(filename,'z',permute(displ,[2 1]),[1 1 t])
    end
    displOld = displ;
end
