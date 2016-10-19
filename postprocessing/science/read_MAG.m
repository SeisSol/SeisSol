%%
% @file
% This file is part of SeisSol.
%
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
%
% @section LICENSE
% Copyright (c) 2014, SeisSol Group
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
% Reads -MAG- files.
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  read_MAG reads in the -MAG- files that include the          %%')
disp('     %%  scalar seismic moment information of each MPI domain        %%')
disp('     %%  that includes + fault elements.                             %%')
disp('     %%                                                              %%')
disp('     %%  INFO:                                                       %%')
disp('     %%  The number in the MAG file for each MPI domain              %%')
disp('     %%  with + fault elements is:                                   %%')
disp('     %%  scalar seismic moment = shear modulus * area * average slip %%')
disp('     %%  whereas shear modulus : from the + element [Pa]             %%')
disp('     %%  average slip : accumulated slip for every GP divided by     %%')
disp('     %%                 number of GPs for this element [m]           %%')
disp('     %%          area : triangular surface [m^2]                     %%')
disp('     %%                                                              %%')
disp('     %%  Script expects that all MAG files are stored within one     %%')
disp('     %%  folder without any MAG file from another simulation.        %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;

% read in data
liste=dir('*-MAG-*');
files = {liste.name};
ntotal = numel(files);
disp(['Found ',num2str(ntotal),' -MAG- files.']),disp(' ')
seismicmoment = 0.0;
for k=1:ntotal
    % read files
    fid = fopen(files{k},'r');
    data_tmp = fscanf(fid,'%g');
    fclose(fid);
    seismicmoment=seismicmoment + data_tmp;
    clear data_tmp
end

seismicmoment_dynecm = seismicmoment*10^7;

%result:
disp(['Total scalar seismic moment is: ',num2str(seismicmoment),' [Nm]'])
disp(['Total scalar seismic moment is: ',num2str(seismicmoment_dynecm),' [dyne cm]'])
disp(' ')

%Calculate moment magnitude
mag = (2/3)*log10(seismicmoment_dynecm)-10.7;
disp(['The moment magnitude M_w is:  ',num2str(mag)])
