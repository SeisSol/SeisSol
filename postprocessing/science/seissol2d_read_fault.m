function data = seissol2d_read_fault_parallel(name,nCPU,sp_sample)
%%
% @file
% This file is part of SeisSol.
%
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
%
% @section LICENSE
% Copyright (c) 2009, SeisSol Group
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

ndat=4; %Number of variables

disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                             %%')
disp('     %% seissol2d_read_fault reads in the *-flt_header.txt and      %%')
disp('     %% *-flt.dat files produced by SeisSol from the current        %%')
disp('     %% directory and plots the slip rate along a fault caused by   %%')
disp('     %% dynamic rupture.                                            %%')
disp('     %%                                                             %%')
disp('     %% needs:     root name = *                                    %%')
disp('     %%                 nCPU = Total number of MPI domains          %%')
disp('     %%            sp_sample = Number of output points per element  %%')
disp('     %%                                                             %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

disp(' ')

count = 0; data.nx=0;

for iCPU = 0:nCPU-1
    
    %define character string of input files
    if(iCPU<10)
        hdr = [name,'-flt_header.000',num2str(iCPU),'.txt'];
    elseif(iCPU>=10 && iCPU<100)
        hdr = [name,'-flt_header.00',num2str(iCPU),'.txt'];
    else
        hdr = [name,'-flt_header.0',num2str(iCPU),'.txt'];
    end
    
    % check if file exists
    check = exist(hdr,'file');
    if(check~=2)
        continue % step out the loop
    end   
    
    %Read header file
     %hdr = strcat(name,'-flt_header.',iCPU,'.txt');
    [temp.nx] = textread(hdr,'%n',1,'headerlines',1);
    [data.dt] = textread(hdr,'%n',1,'headerlines',3);
    [data.ord] = textread(hdr,'%n',1,'headerlines',5);

    if(sp_sample==1)
        [temp.p,temp.x,temp.y] = textread(hdr,'%f%f%f','headerlines',7);
    else
        [tmp.p,tmp.x,tmp.y] = textread(hdr,'%f%f%f','headerlines',7);
        temp.p=tmp.p(1:sp_sample:end);
        temp.x=tmp.x(1:sp_sample:end);
        temp.y=tmp.y(1:sp_sample:end);
        clear tmp;
    end

    %Read matrix of fault data
    if(iCPU<10)
        dat = [name,'-flt-000',num2str(iCPU),'.dat'];
    elseif(iCPU>=10 && iCPU<100)
        dat = [name,'-flt-00',num2str(iCPU),'.dat'];
    else
        dat = [name,'-flt-0',num2str(iCPU),'.dat'];
    end
     %dat  = strcat(name,'-flt.dat');
    fid=fopen(dat);
    raw=fscanf(fid,'%g',inf);
    %raw = fread(fid,[data.nx+2,inf],'single') ;
    fclose(fid);
    %raw = reshape(raw,[],ndat,data.nx);
    raw=shiftdim(raw,1);
    raw = reshape(raw,ndat,temp.nx,[]);

    temp.nx=length(temp.p);

    % Reformat each field [nx,nt]
    temp.v  = squeeze(raw(1,1:sp_sample:end,:));
    temp.st = squeeze(raw(2,1:sp_sample:end,:));
    temp.sn = squeeze(raw(3,1:sp_sample:end,:));
    temp.u  = squeeze(raw(4,1:sp_sample:end,:));
    
    % add to data struct
    count = count + 1;
    if(count == 1)
        data.p  = [temp.p];
        data.x  = [temp.x];
        data.y  = [temp.y];
        data.v  = [temp.v];
        data.st = [temp.st];
        data.sn = [temp.sn];
        data.u  = [temp.u];
    else
        data.p  = vertcat(data.p, temp.p);
        data.x  = vertcat(data.x, temp.x);
        data.y  = vertcat(data.y, temp.y);
        data.v  = vertcat(data.v, temp.v);
        data.st = vertcat(data.st, temp.st);
        data.sn = vertcat(data.sn, temp.sn);
        data.u  = vertcat(data.u, temp.u);
    end
    
    data.nx = data.nx + temp.nx;
    
    clear temp;
end

[b data.nt]=size(data.v(1,:));
data.time = [data.dt:data.dt:data.dt*data.nt];

disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                             %%')
disp('     %%                         OUTPUT                              %%')
disp('     %%                                                             %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                             %%')
disp('     %%                 Matlab Structure Type                       %%')
disp('     %%                                                             %%')
disp('     %%      sample points                          : nx            %%')
disp('     %%      timestep                               : dt            %%')
disp('     %%      ADER DG order                          : ord           %%')
disp('     %%      Positions within element               : p             %%')
disp('     %%      Position in domain                     : x,y           %%')
disp('     %%      Slip rate                              : v             %%')
disp('     %%      Shear traction                         : st            %%')
disp('     %%      normal stress                          : sn            %%')
disp('     %%      normal velocity                        : u             %%')
disp('     %%      Nr of timesteps                        : nt            %%')
disp('     %%      time vector                            : time          %%')
disp('     %%                                                             %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

disp(' ')
