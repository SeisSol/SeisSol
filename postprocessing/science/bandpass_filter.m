function s_filt = bandpass_filter(s,dt,cf_lowleft,cf_lowright,c_highcut,deg)
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
% BANDPASS_FILTER applies a band-pass filter to a time signal 
% given by the "colums" vector s.
% dt is the sampling rate of the signal s.
% cf_lowleft  is the left corner frequency up to which the filter function is zero.  
% cf_lowright is the right corner frequency beyond which the fileter is 1.
% c_highcut   is the corner frequency beyond which the high frequencies are tapered by the
%             degree deg butterworth filter.
% deg         is the degree of the butterworth filter
% Between cf_lowleft and cf_lowright the filter function increases with as the
% cosine-function.

s_length = length(s);

% check for odd number of samples and make it even 
% by skipping the last sample
if(mod(s_length,2)==1)
    s_length = s_length-1;
    s = s(1:end-1);
end

% Fourier transform the signal
F = fft(s);
% Frequency sampling rate
df = 1/(s_length*dt);

freq = df * (0:s_length-1)';

filter_func1 = 0.5*( 1 + cos(pi*(freq-cf_lowright)/(cf_lowright-cf_lowleft)) );
filter_func1(freq<cf_lowleft)  = 0;
filter_func1(freq>cf_lowright) = 1;

filter_func2 = 1./sqrt(1+(freq/c_highcut).^(2*deg));
%combination of high- and low-pass filter
for i=1:length(filter_func1)
   filter_func(i,1) = min(filter_func1(i),filter_func2(i)); 
end
%plot(freq,filter_func)

% multiply spectrum F with filter function 
s_filt = F.*filter_func;

% make filtered signal symmetric
s_filt(s_length/2+2:end) = conj( s_filt(s_length/2:-1:2) );

% Apply inverse Fourier transform
s_filt = real(ifft(s_filt));
