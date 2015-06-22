%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2006, SeisSol Group
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
home;
disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%               MATLAB2GAMBIT                         %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

clear, close all;

filename = input(' MATLAB filename:  ','s');

%%%% READ IN MATLAB MESH FILE
fid   = fopen([filename,'.msh']);
tmp   = fscanf(fid,'%i',[2,1]);   n_vert = tmp(1); n_tri = tmp(2);
tri   = fscanf(fid,'%i',[3,n_tri]);   tri = tri';
X     = fscanf(fid,'%g',[3,n_vert]);  X   = X';
fclose(fid);

fid_out  = fopen([filename,'.neu'],'w');
fprintf(fid_out,'%s','        CONTROL INFO 2.0.0');
fprintf(fid_out,'\n%s','** GAMBIT NEUTRAL FILE');
fprintf(fid_out,'\n%s',filename);
fprintf(fid_out,'\n%s','PROGRAM:                Gambit     VERSION:  2.0.0');
fprintf(fid_out,'\n%s',date);
fprintf(fid_out,'\n%s','     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL');
fprintf(fid_out,'\n%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f',n_vert,n_tri,1,1,2,2);
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','   NODAL COORDINATES 2.0.0');
counter = (1:n_vert)'; 
fprintf(fid_out,'\n%10.0f%20.10e%20.10e',[counter,X(:,1:2)]');
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','      ELEMENTS/CELLS 2.0.0');
counter = (1:n_tri)';
fprintf(fid_out,'\n%8.0f%3.0f%3.0f%9.0f%8.0f%8.0f',[counter,ones(size(counter))*3,ones(size(counter))*3,tri]');
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','       ELEMENT GROUP 2.0.0');
fprintf(fid_out,'\n%s','GROUP:          1 ELEMENTS:');
fprintf(fid_out,'%11.0f%',n_tri);
fprintf(fid_out,'%s',' MATERIAL:          2 NFLAGS:          1');
fprintf(fid_out,'\n%s','                           solid');
fprintf(fid_out,'\n%s','       0');
fprintf(fid_out,'\n%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f',(1:1:n_tri));
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s',' BOUNDARY CONDITIONS 2.0.0');

trimesh(tri,X(:,1),X(:,2),zeros(n_vert,1)), view(0,90), axis square, hold on
b = find(X(:,3)==101);
plot(X(b,1),X(b,2),'r.')

Btri   = ismember(tri,b);
sumtri = sum(Btri,2);
b = find(sumtri==2);
sumedge(:,1) = sum(Btri(b,[1,2]),2);
sumedge(:,2) = sum(Btri(b,[2,3]),2);
sumedge(:,3) = sum(Btri(b,[1,3]),2);
[tmp,edge]   = find(sumedge==2);
tmp          = sortrows([tmp,edge]);
edge         = tmp(:,2); 
CX = ( X(tri(:,1),1)+X(tri(:,2),1)+X(tri(:,3),1) ) / 3;
CY = ( X(tri(:,1),2)+X(tri(:,2),2)+X(tri(:,3),2) ) / 3;
plot(CX(b),CY(b),'b.')

fprintf(fid_out,'\n%32.0f%8.0f%8.0f%8.0f%8.0f',101,1,length(b),0,6);
fprintf(fid_out,'\n%10.0f%5.0f%5.0f',[b,ones(size(b))*3,edge]');
fprintf(fid_out,'\n%s\n','ENDOFSECTION');


fclose(fid_out);
