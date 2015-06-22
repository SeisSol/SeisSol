function [Vertex,Conectivity]=File_ts2STL(name)
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

disp(' '),disp(' ')
disp('File_ts2STL')
disp(' '),disp(' ')
disp('converts GoCads ts files in STL files read by ICEM')

% Open file to read
hdr = strcat(name,'.ts');
fid = fopen(hdr);
% Read the full file character by character and store the ascii code in one
% single array
Ascii = fread(fid);
% Close file
fclose(fid);



nascii = length(Ascii);
icolumn =1; irow = 1;
buildnumber = 0; assignnumber = 0; number = 0; power = 0;
sc = 0; neg=1;
for iascii = 1:nascii
    data = Ascii(iascii);
    if data==10         %new line
        if sc==0
            irow      = irow+1;
            icolumn   = 1;
        else
            assignnumber=1;
        end
    elseif data==32     %space
        assignnumber = 1;
    elseif data==45     %-
        neg = -1;
    elseif data==46     %.
        sign = -1;
    elseif data>=48 && data<=57 %decimal numbers
        buildnumber = 1;
    elseif data==76 && sc == 0  %L
        sc = 1;
        icolumn =1; irow = 1;
        buildnumber = 0; assignnumber = 0; number = 0; power = 0;
    end

    if buildnumber==1
        digit = neg*str2num(char(data));
        if sign==1
            number = number*10 + digit;
        elseif sign==-1
            power  = power+1;
            number = number + digit/(10^power);
        end
        buildnumber = 0;
    end
    if assignnumber==1
        if sc==0
            Vertex(icolumn,irow) = number;
        else
            Conectivity(icolumn,irow) = number;
        end
        icolumn      = icolumn+1;
        sign         = 1;
        power        = 0;
        number       = 0;
        assignnumber = 0;
        neg          = 1;
        if data==10 && sc==1 
            irow      = irow+1;
            icolumn   = 1;
        end        
    end

end

file = strcat(name,'.stl');
fid=fopen(file,'w');
fprintf(fid,'solid ');
fprintf(fid,name);
fprintf(fid,'\n');

for ielem=1:length(Conectivity)   
    x=Conectivity(2,ielem);
    y=Conectivity(3,ielem);
    z=Conectivity(4,ielem);
    fprintf(fid,'facet normal 0 0 0\n');
    fprintf(fid,'outer loop\n');
    fprintf(fid,'vertex%20.11e%20.11e%20.11e\n',Vertex(3,x),Vertex(4,x),Vertex(5,x));
    fprintf(fid,'vertex%20.11e%20.11e%20.11e\n',Vertex(3,y),Vertex(4,y),Vertex(5,y));
    fprintf(fid,'vertex%20.11e%20.11e%20.11e\n',Vertex(3,z),Vertex(4,z),Vertex(5,z));
    fprintf(fid,'endloop\n');
    fprintf(fid,'endfacet\n');
end
fprintf(fid,'endsolid ');
fprintf(fid,name);
fclose(fid)
return




     





