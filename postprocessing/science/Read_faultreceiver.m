%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
%
% @section LICENSE
% Copyright (c) 2006-2011, SeisSol Group
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
% Reformat_Faultreceiver reads in the TECPLOT format of the typical SEISSOL seismogram output and converts it into a MATLAB struct format.
% Explanation in the struct.

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  Reformat_Faultreceiver reads in the TECPLOT format          %%')
disp('     %%  of the typical SEISSOL seismogram output and converts       %%')
disp('     %%  it into a MATLAB struct format.                             %%')
disp('     %%                                                              %%')
disp('     %%  Explanation in the struct.                                  %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
workingdir   = input('     Give relative path, where are the data                :  ','s');
if strcmp(workingdir,'')
    workingdir='.'
end

filename = input('     Specify root-filename (typically: **file**-faultreceiver-, without the -):  ','s');

%evaluate nb of files and nb of time samples
%do not consider files with only header
eval(['!find ', workingdir,' -maxdepth 1 -name "',filename,'*faultreceiver*" -size +50k | xargs wc -l > tmp.dat']);



tmp_file = 'tmp.dat';
fid   = fopen(tmp_file);
command_result = textscan(fid,'%d %s',[inf 2]);
fclose(fid);
eval(['!rm tmp.dat']);

[nseis,y]=size(command_result{1});
nseis=nseis-1;

if nseis==0
    error('no seismogram found')
end
%count the number of lines starting with '#'
[status,coutput]=system(['awk ''/#/{a=NR}END{print a}'' ',command_result{2}{1}]);
last_header_line=str2num(coutput);
ndt=min(command_result{1}(1:end-1))-last_header_line;

if (last_header_line>8)
    disp('Warning: more than 8 header line found')
end

msg = sprintf('found %d seismogram(s)', nseis)
disp(msg);

msg = sprintf('found %d time sample(s) (last header line: %d)', ndt, last_header_line);
disp(msg);

%nseis    = input('     Give number of last seismogram:                             ');
%ndt      = input('     Give number of total time samples of the data:              ');
sr       = input('     Data includes - slip rate  - information?  (1= yes, 0= no)  ');
stress   = input('     Data includes - stress     - information?  (1= yes, 0= no)  ');
vel      = input('     Data includes - normal vel - information?  (1= yes, 0= no)  ');
RS       = input('     Data includes - Mu and SV  - information?  (1= yes, 0= no)  ');

 % 1 for time
nvar = 1;
char = 'VARIABLES = "Time"';
if (sr == 1);
    nvar = nvar + 2;
    char = [char,' ,"SRs" ,"SRd"'];
end;
if (stress == 1);
    nvar = nvar + 3;
    char = [char,' ,"T_s" ,"T_d" ,"P_n"'];
end;
if (vel == 1);
    nvar = nvar + 1;
    char = [char,' ,"u_n"'];
end;
if (RS == 1);
    nvar = nvar + 2;
    char = [char,' ,"Mud"','StV'];
end;
disp(' '), disp(sprintf('    \n Total number of variables including time:  %i\n',nvar));


num = 1;
ind = 0;
location(3,1:nseis) = 0;
backgroundstress(3,1:nseis) = 0;

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
    junk  = fscanf(fid,'%c',[7,1]);
    P_0   = fscanf(fid,'%g',[1,1]);
    junk  = fscanf(fid,'%c',[7,1]);
    T_s   = fscanf(fid,'%g',[1,1]);
    junk  = fscanf(fid,'%c',[7,1]);
    T_d   = fscanf(fid,'%g',[1,1]);

    %SKIP extra header lines
    for i=8:last_header_line
        junk  = fgetl(fid);
    end


    data = fscanf(fid,'%g',[nvar,ndt]); data = data';

    disp(sprintf('\n Samples in seismogram nr. %i\t(%s):   %i', num,in_file,size(data,1)));
    MATLAB_OUTPUT(:,:,num) = data;
    location(:,num) = [x1; x2; x3];
    backgroundstress(:,num) = [P_0; T_s; T_d];
    fclose(fid);
      
end

% add to struct
d.data=MATLAB_OUTPUT;d.location=location;d.bg=backgroundstress;
d.datalegend=char;
d.locationlegend=('X,Y,Z');
d.bglegend=('Background Stress: P_0, T_s, T_d');

% output
disp(' '),  disp('     Finished conversion!'), disp(' ')
out_filename = input('     Give filename to save MATLAB data (add .mat):  ','s');
out_file = out_filename(1:end-4);
save(out_filename,'-v7.3','d')
disp(sprintf('    \n Saved data as :  %s\n',out_filename));
