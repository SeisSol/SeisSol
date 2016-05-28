%%
% @file
% This file is part of SeisSol.
%
% @author Stephanie Wollherr
%
% @section LICENSE
% Copyright (c) 2016, SeisSol Group
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
% Reads -ENERGY- files.
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  read_EN reads in the -EN- files that include the            %%')
disp('     %%  energies released in each MPI domain.                       %%')
disp('     %%                                                              %%')
disp('     %%                                                              %%')
disp('     %%  INFO:                                                       %%')
disp('     %%  The numbers in the EN file for each MPI domain              %%')
disp('     %%  is a time series of the corresponding energies              %%')
disp('     %%  kinetic energy, plastic energy and elastic strain energy    %%')
disp('     %%  summation over each element in the domain                   %%')
disp('     %%                                                              %%')
disp('     %%  Script expects that all EN files are stored within one      %%')
disp('     %%  folder without any EN file from another simulation.         %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;

plast = input('Simulation with plasticity? 1=yes, 0=no ');
pltOption = input('Plot all energies in one plot? 1=yes, 0=no ');

if plast==0
    nvar = 2; %time kineticEnergy
elseif plast ==1
    nvar = 4; %time kinetic Energy plasticEnergy estrainEnergy
end
 
% read in data
liste=dir('*-EN*');
files = {liste.name};
ntotal = numel(files);
disp(['Found ',num2str(ntotal),' -EN- files.']),disp(' ');


fid   = fopen(files{1},'r');
fgets(fid);
data_tmp = fscanf(fid,'%g',[nvar,inf] );
fclose(fid);
%number of timesteps
[tmp ndt] = size(data_tmp);

%
disp(['Found ', num2str(ndt), ' timesteps.']);
energy_total = zeros(nvar-1,ndt);

%structure of -EN-files
%1=time; 2=kinetic; (3=plastic; 4=elastic strain)

for k=1:ntotal
    % read files
    fid = fopen(files{k},'r');
    fgets(fid);  % Ignore first line
    data_tmp = fscanf(fid,'%g',[nvar,ndt] );
    fclose(fid);
    energy_total = energy_total + data_tmp(2:nvar,:); %sum for every timestep
    if k==1
        time = data_tmp(1,1:ndt);
    end
    clear data_tmp
end

dt = time(4)-time(3);

if plast==1
    acc_energy = zeros(1,ndt);
    %integral over time to get the total plastic energy change
    for i=2:ndt
        acc_energy(i)= acc_energy(i-1) + dt*energy_total(2,i);
    end
end

%plotting section

if plast==1
   if pltOption==1 %all values in one plot

    plot(time,energy_total(1,:)); hold on;
    plot(time,acc_energy);
    plot(time,energy_total(3,:));
    xlabel('time');
    legend('kinetic', 'plastic','elastic strain', 'Location', 'northwest');
    title('Energies')

   elseif pltOption==0 %all values seperately

    subplot(4,1,1)
    plot(time,energy_total(1,:));
    title('kinetic energy');
    xlabel('time')
    
    subplot(4,1,2)
    plot(time,energy_total(2,:));
    title('dissipated plastic energy change');
    xlabel('time')

    subplot(4,1,3)
    plot(time,acc_energy);
    title('dissipated plastic energy');
    xlabel('time')

    subplot(4,1,4)
    plot(time,energy_total(3,:));
    title('elastic strain energy');
    xlabel('time')
   end
    
elseif plast==0
    plot(time,energy_total(1,:));
    title('kinetic energy');
    xlabel('time')
end

