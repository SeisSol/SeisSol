function [am,f,ph] = sfft(y,dt,whole) 
%%
% @file
% This file is part of SeisSol.
%
% @author Martin Mai (Martin.Mai AT kaust.edu.sa, http://www.kaust.edu.sa/faculty/mai.html)
%
% @section LICENSE
% Copyright (c) 1998, SeisSol Group
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
%  function [am,f,ph] = sfft(y,dt,whole) 
%  calculates the Fourier Transform of a signal
% 	
%  Input: y  - signal; can be array of many signals of equal
%	       length and sampling, arranged in columns
%	  dt - sampling interval or time axis
%	  whole - optional: flag to get entire spectrum ('y')
%
%  Output: am - amplitude spectrum
%	   f  - frequency axis for pos. frequencies only
%	   ph - phase spectrum

if nargin == 1; dt = 0.01; whole = 'n';
elseif nargin == 2; whole = 'n'; 
end;

m = size(y,1);

if size(y,1) == 1, y = y(:); end
if mod(m,2)  == 1
    m=m-1; y = y(1:m,:);
end


if length(dt) > 1
   dt = dt(2)-dt(1);		% assuming dt is regular spaced t-axis
end

W = 1/dt;			% sampling freq in Hz 
f = (W/m)*linspace(-m/2,m/2,m)';% frequency axis


%%% calculate absolute FFT, scaled by sampling interval 
%%% and shifted symmetrically

FT = (fft(y))*dt;
for i = 1:size(y,2);
  am(:,i) = abs(fftshift(FT(:,i)));
  ph(:,i) = angle(FT(:,i));
end

if whole == 'n';
  AM = am; PH = ph;
  clear am ph;
  f = f(m/2+1:end);
  for i = 1:size(y,2)
    am(:,i) = AM(m/2+1:end,i);
    ph(:,i) = PH(m/2+1:end,i);
  end
end
