function inside = XYinElement(x,y,xS,yS)
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
% XYinElement
%  XYInElement returns 1 if the point (xS/yS) is inside the
%  element defined by XY, else it returns 0
%  adapted from SeisSol2Ds XYinElement function

epsilon = 0.00001; % tolerance 

refFactor = 0.5;
VOL = volume(x,y);

xi  = refFactor/VOL*( ( x(3)*y(1) - x(1)*y(3) ) + xS*(y(3)-y(1)) + yS*(x(1)-x(3)) );
eta = refFactor/VOL*( ( x(1)*y(2) - x(2)*y(1) ) + xS*(y(1)-y(2)) + yS*(x(2)-x(1)) );

      
if (xi <(0.0-epsilon)) || (eta <(0.0-epsilon)) || (eta > (1.0-xi+epsilon))
  inside = 0;
else
  inside = 1;
end

function VOL = volume(x,y)
% Sarrus
VOL = 0.5*( (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)) );
