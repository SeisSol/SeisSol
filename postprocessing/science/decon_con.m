function signal2 = decon_con(wavelet1,signal1,wavelet2,dt,cf_left,cf_right)
%%
% @file
% This file is part of SeisSol.
%
% @section LICENSE
% Copyright (c) SeisSol Group
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
% DECON_CON 
% function signal2 = decon_con(wavelet1,signal1,wavelet2,dt,cf_left,cf_right)
% Computes time signal2 (with input avelet2) out of signal1 (with input wavelet1)
% The arguments wavelet1, wavelet2 and signal1 have to be column vectors with sampling rate dt.
% A lowpass filter is applied on signal2.
%
% cf_left  is the left corner frequency up to which the filter function is 1.
% cf_right is the right corner frequency beyond which the filter is zero.
% Between cf_left and cf_right the filter function decays with as the cosine-function.

% Fourier Transformation
F1=fft(wavelet1);
F2=fft(signal1);
G1=fft(wavelet2);

% Deconvolution 
amg=F2./F1; 

% Convolution 
G2= G1.*amg;
signal_length = length(G2);

% check for odd number of samples and make it even 
% by skipping the last sample
if(mod(signal_length,2)==1)
    signal_length = signal_length-1;
    G2 = G2(1:end-1);
end

% Frequency sampling rate
df = 1/(signal_length*dt);

freq = df * (0:signal_length-1)';

filter_func = 0.5*( 1 + cos(pi*(freq-cf_left)/(cf_right-cf_left)) );
filter_func(freq<cf_left)  = 1;
filter_func(freq>cf_right) = 0;

% multiply spectrum G2 with filter function 
G2 = G2.*filter_func;

% make filtered signal symmetric
G2(signal_length/2+2:end) = conj( G2(signal_length/2:-1:2) );

% Apply inverse Fourier transform and apply baseline correction
signal2 = real(ifft(G2));
signal2=signal2(:)-signal2(1);
