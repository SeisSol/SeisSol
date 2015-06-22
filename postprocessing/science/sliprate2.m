%%
% @file
% This file is part of SeisSol.
%
% @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
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

clear;clc;home;

disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                             %%')
disp('     %% SlipRate.m reads in the .SR. files produced by SeisSol from %%')
disp('     %% the current directory and plots the slip rate along a fault %%')
disp('     %% caused by dynamic rupture.                                  %%')
disp('     %%                                                             %%')
disp('     %%    NOTE: Make sure only .SR. files from ONE simulation are  %%')
disp('     %%          present in the current directory!                  %%')
disp('     %%                                                             %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

disp(' ')

clear; close all;

disp(' ')

name=input(' Enter -flt_header name:');
%Read header file
hdr = strcat(name,'-flt_header.txt');
[data.nx] = textread(hdr,'%n',1,'headerlines',1);
[data.dt] = textread(hdr,'%n',1,'headerlines',3);
[data.ord] = textread(hdr,'%n',1,'headerlines',5);
[data.p,data.x,data.y] = textread(hdr,'%f%f%f','headerlines',7);

FileEx=input(' Enter first file name:');
disp(' ')
filelength=length(FileEx);

%Finding where is the ordering information in the file names (between dots
%second and third of file names)
count=0;
for i=1:filelength
    if(FileEx(i)=='.')
        count=count+1;
        if count==2
            num_ini=i+1;
        elseif count==3;
            num_end=i-1;
        end
    end
end

disp(' Reading all .SR. files in directory...')
[R]=evalc('ls *SR*');

j=0;
filecount=0;
FileList=cell(1);
for i=1:length(R)
    j=j+1;
    if(R(j)~='  ')
        filecount=filecount+1;
        filenow=R(j:j+filelength);
        filenow=strtrim(filenow);
        FileList(filecount)=cellstr(filenow);
        FileList=[FileList;cell(1)];
        ord(filecount)=str2num(filenow(num_ini:num_end));
        j=j+filelength;
    else
        continue
    end
    if (abs(length(R)-j)<=filelength)
        break
    end
end

disp([' Found ',num2str(filecount),' files of Slip Rate in the directory'])

%Reordering of files
[B indx]=sort(ord);
FileList2=FileList;
for i=1:filecount
    FileList(i)=FileList2(indx(i));
end

disp(' Files successfully found and reordered!')
disp(' ')

maxtot=0;
sgntmp=1;
txtadd=' ';
for i=1:filecount    
    fid = fopen(char(FileList(i)),'r');
    
    for k=1:6
        tmp=fgets(fid);
    end

    % Reading file header
    tmp=fgets(fid);time=str2num(tmp(1:22));it=str2num(tmp(23:end));
    tmp=fgets(fid);tmp=fgets(fid);npoints=str2num(tmp);    
    tmp=fgets(fid);tmp=fgets(fid);order=str2num(tmp(1:22));npoints=str2num(tmp(23:end));    
    tmp=fgets(fid);

    S=fscanf(fid,'%f');
    
    for j=1:length(S)/2
        X(j) = S(2*j-1);
        SR(j)= S(2*j);
    end
    
    fclose('all');

    if(time>0)    
        maxtmp=max(abs(SR));
        if(max(SR)>=-min(SR))
            sgntmp=1;
        else
            sgntmp=-1;
        end
        if(maxtmp>maxtot)
            maxtot=maxtmp;
        end
        if(sgntmp==-1)
            txtadd='  sign of slip rate changed!';
        end
    end
    disp(['Plotting file ',char(FileList(i)),'...',txtadd])        
    subplot(211)
    plot(data.x,SR*sgntmp,'LineWidth',2),title(['Time=',num2str(time),' sec'],'FontSize',12),ylabel('Slip Rate (m/s)','FontSize',12)
    set(gca,'FontSize',12)        
    if(time>0)            
        ylim([-1.2*maxtot 1.2*maxtot]);
    end
 
    subplot(212)
    plot(data.x,SR*sgntmp,'LineWidth',2),xlabel('Infault position (m)','FontSize',12),ylabel('Slip Rate (m/s)','FontSize',12)
    set(gca,'FontSize',12)     
%    print -depsc SlipRate
        
    pause;

end

clear; close all;
disp(' ')
disp(' SlipRate output successfully finished!')
