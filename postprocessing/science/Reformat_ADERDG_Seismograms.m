%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2006, SeisSol Group
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
% Reads in TECPLOT format of SeisSol's seismogram output and converts it.
home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  Reformat_ADER-DG_Seismograms reads in the TECPLOT format    %%')
disp('     %%  of the typical SEISSOL seismogram output and converts       %%')
disp('     %%  it into a 3-dimensional tensor in MATLAB format.            %%')
disp('     %%                                                              %%')
disp('     %%  Dimension 1 has the length of the number of time samples    %%')
disp('     %%  Dimension 2 has the length of the number of variables + 1   %%')
disp('     %%              (in general this is 10: 1 time     component    %%')
disp('     %%                                      6 stress   components   %%')
disp('     %%                                      3 velocity components   %%')
disp('     %%  Dimension 3 has the length of the number of receivers       %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
workingdir   = input('     Give relative path, where are the data                :  ','s');
if strcmp(workingdir,'')
    workingdir='.'
end
prefix       = input('     Specify root-filename (typically: *file*-receiver-...):  ','s');
nvar         = input('     Give number of variables (including time):            ');


%evaluate nb of files and nb of time samples
%do not consider files with only header
eval(['!find ', workingdir,' -maxdepth 1 -name "',prefix,'*-receiver*" -size +30k | xargs wc -l > tmp.dat']);
tmp_file = 'tmp.dat';
fid   = fopen(tmp_file);
command_result = textscan(fid,'%d %s',[inf 2]);
fclose(fid);
eval(['!rm tmp.dat']);

[nseis,y]=size(command_result{1});
nseis=nseis-1;
disp('assuming 5 lines of header for computing ndt');
ndt=min(command_result{1}(1:end-1))-5;

msg = sprintf('found %d seismogram(s)', nseis);
disp(msg);
msg = sprintf('found %d time sample(s)', ndt);
disp(msg);

if nseis==0
    error('no seismogram found')
end

num = 1;
ind = 0;
location(3,1:nseis) = 0;

for num = 1:nseis
     
    in_file = command_result{2}{num};

    fid   = fopen(in_file);
    junk  = fgetl(fid); junk  = fgetl(fid); 
    junk  = fscanf(fid,'%c',[6,1]);
    x1    = fscanf(fid,'%g',[1,1]);
    junk  = fscanf(fid,'%c',[6,1]);
    x2    = fscanf(fid,'%g',[1,1]);
    junk  = fscanf(fid,'%c',[6,1]);
    x3    = fscanf(fid,'%g',[1,1]);
    disp(sprintf('X: %g   Y: %g  Z: %g',x1,x2,x3));
    location(:,num) = [x1; x2; x3];
    data = fscanf(fid,'%g',[nvar,ndt]); data = data';
    
    %determine and delete double samples from a possible restart
    data(:,1) = (round(data(:,1)*1000000))/1000000;
    [tmp,k,kk] = unique(data(:,1));
    data  = data(k,:);
    
    disp(sprintf('\n Samples in seismogram nr. %i\t:   %i', num,size(data,1)));
    MATLAB_OUTPUT(:,:,num) = data;
    fclose(fid);
      
end

% add to struct
d.data=MATLAB_OUTPUT;
d.location=location;
d.datalegend=char;
d.locationlegend=('X,Y,Z');


disp(' '),  disp('     Finished conversion!'), disp(' ')
out_filename = input('     Give filename to save MATLAB data (adds .mat):  ','s');
out_file = out_filename(1:end-4);
save(out_filename,'-v7.3','d')
disp(sprintf('    \n Saved data as :  %s\n',out_filename));
