%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2005, SeisSol Group
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
disp('    %%                 GAMBIT_RED_REFINE (2D)              %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

clear, close all;

filename = input('    Filename(suffix ".neu" assumed):  ','s');
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Reading data from: \t\t\t%s',[filename,'.neu']));
disp(' ')
fid_in   = fopen([filename ,'.neu']);
fid_out  = fopen([filename ,'_refined.neu'],'w');
junk  = fgetl(fid_in);      fprintf(fid_out,'%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
param = fscanf(fid_in,'%i',[6,1]);  
NX = param(1);  NT = param(2); NG = param(3); NBSETS = param(4); NDFCD = param(5); NDFVL = param(6); 
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);      
junk  = fgetl(fid_in);    
X     = fscanf(fid_in,'%g',[3,NX]); X = X'; X(:,1) = [];
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
tri   = fscanf(fid_in,'%g',[6,NT]); tri = tri'; tri(:,1:3) = [];
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fscanf(fid_in,'%g',[10,ceil(NT/10)]);
%junk  = fgetl(fid_in);
%junk  = fgetl(fid_in);
%param = fscanf(fid_in,'%i',[5,1]);
%bc    = fscanf(fid_in,'%g',[3,param(3)]); bc = bc';
%junk  = fgetl(fid_in);
%junk  = fgetl(fid_in);
fclose(fid_in);

disp(sprintf('\t Read \t\t\t\t\t%i vertices',NX));
disp(sprintf('\t Read \t\t\t\t\t%i triangular elements',NT));

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Start refinement ...'));

subplot(121)
h = trimesh(tri,X(:,1),X(:,2),zeros(length(X),1));view(0,90),axis square
set(h,'EdgeColor','k')
xlabel('x','FontWeight','bold','FontSize',12)
ylabel('y','FontWeight','bold','FontSize',12,'Rotation',0)
axis([-1 1 -1 1 0 1]), view(0,90), axis square, drawnow;

tin=cputime;

MX = [(X(tri(:,1),1)+X(tri(:,2),1)),(X(tri(:,2),1)+X(tri(:,3),1)),(X(tri(:,3),1)+X(tri(:,1),1))]/2;
MY = [(X(tri(:,1),2)+X(tri(:,2),2)),(X(tri(:,2),2)+X(tri(:,3),2)),(X(tri(:,3),2)+X(tri(:,1),2))]/2;
X_ref  = [];
bc_ref = [];
c=0;
for t = 1:NT
    for n = 1:4
        c = c+1;        
        switch n
            case 1
                x = [X(tri(t,1),:);MX(t,1),MY(t,1);MX(t,3),MY(t,3)];
            case 2
                x = [X(tri(t,2),:);MX(t,2),MY(t,2);MX(t,1),MY(t,1)];
            case 3
                x = [X(tri(t,3),:);MX(t,3),MY(t,3);MX(t,2),MY(t,2)];
            case 4
                x = [MX(t,1),MY(t,1);MX(t,2),MY(t,2);MX(t,3),MY(t,3)];
        end
        [tmp,io,in] = intersect(X_ref,x,'rows');
        len = size(X_ref,1);
        if length(in) > 0;
            x(in,:)        = [];
            ind            = [io',len+1,len+2]; 
            tri_ref(c,1:3) = ind(1:3);
            X_ref          = [X_ref;x];
        else
            tri_ref(c,1:3) = (len+1:len+3);
            X_ref          = [X_ref;x];  
        end
        % compute counter-clockwise oriented edges
        p = [X_ref(tri_ref(c,1),:),X_ref(tri_ref(c,2),:),X_ref(tri_ref(c,3),:)];
        % compute cross product 
        cr = (p(3)-p(1)).*(p(6)-p(2))-(p(4)-p(2)).*(p(5)-p(1));
        % reverse order of vertices where necessary
        tri_ref(c,:) = fliplr(tri_ref(c,:));

    end
%    bcind = find(bc(:,1)==t);
%    if(length(bcind)>0)
%        bedge = bc(bcind,3);
%        switch bedge
%            case 1
%                bc_new = [c-3;c-2];
%                evec   = X(tri(t,2),:)-X(tri(t,1),:); 
%            case 2
%                bc_new = [c-2;c-1];
%                evec   = X(tri(t,3),:)-X(tri(t,2),:); 
%            case 3
%                bc_new = [c-1;c-3];
%                evec   = X(tri(t,1),:)-X(tri(t,3),:); 
%        end
%        
%        bc_ref = [bc_ref; [bc_new,[3;3],[1;3]]];
%    end 
    
end

tri = tri_ref; X = X_ref;  
[NX,tmp]      = size(X);
[NT,tmp]      = size(tri);
disp(sprintf('\t Writing output-file: \t\t\t%s',[filename,'_refined.neu']));
disp(sprintf('\n\t Write \t\t\t\t\t%i vertices',NX));
disp(sprintf('\t Write \t\t\t\t\t%i triangular elements',NT));

fprintf(fid_out,'\n%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f',NX,NT,NG,NBSETS,NDFCD,NDFVL);
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','   NODAL COORDINATES 2.0.0');
fprintf(fid_out,'\n%10.0f%20.10e%20.10e',[(1:NX)',X]');
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','      ELEMENTS/CELLS 2.0.0');
fprintf(fid_out,'\n%8.0f%3.0f%3.0f%9.0f%8.0f%8.0f',[(1:NT)',ones(NT,1)*3,ones(NT,1)*3,tri]');
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','       ELEMENT GROUP 2.0.0');
fprintf(fid_out,'\n%s','GROUP:          1 ELEMENTS:');
fprintf(fid_out,'%11.0f%',NT);
fprintf(fid_out,'%s',' MATERIAL:          2 NFLAGS:          1');
fprintf(fid_out,'\n%s','                           fluid');
fprintf(fid_out,'\n%s','       0');
fprintf(fid_out,'\n%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f',(1:1:NT));
fprintf(fid_out,'\n%s\n','ENDOFSECTION');
%fprintf(fid_out,'\n%s',' BOUNDARY CONDITIONS 2.0.0');
%fprintf(fid_out,'\n%32.0f%8.0f%8.0f%8.0f%8.0f',param(1),param(2),param(3)*2,param(4),param(5));
%fprintf(fid_out,'\n%10.0f%5.0f%5.0f',bc_ref');
%fprintf(fid_out,'\n%s\n','ENDOFSECTION');
fclose(fid_out);
subplot(122)
h = trimesh(tri,X(:,1),X(:,2),zeros(length(X),1));view(0,90),axis square
%set(gca,'XTick',[-0.5 0 0.5]),set(gca,'YTick',[-0.5 0 0.5])
%set(gca,'XTickLabel',[-0.5 0 0.5],'FontSize',11,'FontWeight','bold')
%set(gca,'YTickLabel',[-0.5 0 0.5],'FontSize',11,'FontWeight','bold')
%xlabel('x','FontSize',12,'FontWeight','bold')
%ylabel('y','FontSize',12,'FontWeight','bold','Rotation',0)
%save tri.mat tri; save X.mat X;
set(h,'EdgeColor','k')
xlabel('x','FontWeight','bold','FontSize',12)
ylabel('y','FontWeight','bold','FontSize',12,'Rotation',0)
%set(gca,'XTick',[-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5],...
%  'XTickLabels',[-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5])
%set(gca,'YTick',[-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5],...
%  'YTickLabels',[-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5])
axis([-1 1 -1 1 0 1]), view(0,90), axis square, drawnow;

tt = cputime-tin;

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Refinement finished successfully!  (%g CPU sec)\n',tt));
disp('-----------------------------------------------------------------------------------')
disp(' ')
