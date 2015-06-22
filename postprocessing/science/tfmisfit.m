function tfmisfit(file,file_ref,keep)
%%
% @file
% This file is part of SeisSol.
%
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
%
% @section LICENSE
% Copyright (c) 2010, SeisSol Group
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
% function to calculate TF misfit (Kristekova, Kristek, Moczo, Day, 2006)
%
% call: tfmisfit('file.mat','file_ref.mat',keep)
%
% file     : file contains all seismograms (mat format - file(samples, [time,u,v,w], nr_index)))
% file_ref : file contains reference seismograms
% keep     : keep all files (1), keep only MISFIT-GOF.DAT (0)
%
% %%%
% TF Misfit parameters have to be defined in HF_TF-MISFIT_GOF!
% this file together with the executable (tf_misfits_gof) should be in the current folder
% %%%

% load data
X1=load(file);
X2=load(file_ref);
var=fieldnames(X1);
d=X1.(var{1});
var=fieldnames(X2);
d_ref=X2.(var{1});
clear X1 X2 var;

[time_samples var nr]=size(d);
[time_samples_r var_r nr_r]=size(d_ref);

% check consistency
if (time_samples ~= time_samples_r)
    disp('ERROR: files contain differtent numbers of time samples!')
    return
elseif(nr ~= nr_r)
    disp('ERROR: files contain differtent numbers of seismograms!')
    return
elseif(var ~=var_r)
    disp('ERROR: files contain differtent numbers of variables!')
    return
end;

for i=1:nr
    % pick actual traces and save them in ascii
    data=d(:,:,i);save data.dat data -ascii;
    ref=d_ref(:,:,i);save ref.dat ref -ascii;
    clear data ref
    
    % calculate TF Misfit
    eval('!./tf_misfits_gof HF_TF-MISFIT_GOF');
    
    % possible cleanings and renamings
    delete('*.dat');delete('S1.DAT');delete('S2.DAT');
    if (keep == 0)
        if (i == 1)
            mkdir('MISFITFILES')
        end;
        % remove all unnecassary files and rename MISFIT-GOF.DAT
        movefile('MISFIT-GOF.DAT',['MISFITFILES/MISFIT_p',num2str(i)]);
        delete('*.DAT');
    elseif(keep == 1)
        % rename all files
        mkdir(['FILES_',num2str(i)]);
        movefile('*.DAT',['FILES_',num2str(i),'/.']);
    end;
    disp('Seismogram done:'),disp(i);
end;
disp('Succesfully finished TF misfit calculation!')
