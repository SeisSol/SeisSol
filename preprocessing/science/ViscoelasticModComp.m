%%
% @file
% This file is part of SeisSol.
%
% @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
%
% @section LICENSE
% Copyright (c) 2007, SeisSol Group
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BRIEF INTRODUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In viscoelastic materials there exists phase velocity dispersion,
%% meaning that the velocities of the waves vary with the frequency. In the
%% viscoelastic model used in SeisSol (Generalized Maxwell Body, GMB), the
%% velocities are maximum for infinity frequency and get lower at lower
%% frequencies, but other models behave differently. The shape of this
%% curves depends on the attenuation chosen (quality factor QP and QS) as
%% well as the bandwidth and central frequency of the attenuating
%% mechanisms.
%%   For most practical applications, we know QP and QS as well as the
%% seismic wave velocities at some particular frequency, often the source's
%% peak frequency. The input of SeisSol, however, is QP and QS and the
%% seismic wave velocities at infinite frequency. In order to transform
%% from one value set to the other we just need to run this script.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The script's input is the following:
%%   QPval,QSval: quality factor of the P and S waves
%%   n:           number of attenuating mechanisms
%%   freq:        frequency at which the phase velocities are defined
%%   f_ratio:     frequency bandwidth of attenuating mechanisms
%%   rho:         the material's density, which is frequency independent
%%   vp_t, vs_t:  target wave velocities at frequency "freq"
%%   vp_0, vs_0:  wave velocities at infinite frequency (unrelaxed
%%                velocities)
%%   In a first run, vp_0 and vs_0 are unknown so we can set them at the same
%% value as vp_t and vs_t. The program will then output in the screen the
%% recommended values for vp_0 and vs_0 given the proposed setup. For
%% double-checking one can then re-run the script now assigning the
%% previously obtained vp_0 and vs_0 values so that we are certain that the
%% values computed are correct. 
home
clear all;
close all;
clc
disp('======================================')
disp('COMPUTATION OF PHASE VELOCITIES AT PEAK')
disp(' FREQUENCY FOR VISCOELASTIC MODELLING  ')
disp('======================================')
disp('   ((c) Josep de la Puente Alvarez  ')
disp('        LMU Geophysics 2007)')
disp('======================================')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attenuation Q-Solver parameters:
QPval=20.0;          % Desired QP constant value
QSval=10.0;          % Desired QS constant value
n=5;                 % Number of mechanisms we want to have

%Frequency parameters:
freq=1;            %  Central frequency of the absorption band (in Hertz). Good to center it
                     %  at the source's central frequency
f_ratio=100;         % The ratio between the maximum and minimum frequencies of our bandwidth
                     % (Usually between 10^2 and 10^4)

%Desired P and S wave velocities at peak frequency:
vp_t=2;           % P-wave velocity
vs_t=1;           % S-wave velocity

%Material parameters at infinity frequency (unrelaxed moduli):
rho=1;            % Density
vp_0=2.120858074089156;    % P-wave velocity
vs_0=1.122397868862394;    % S_wave velocity

disp('--------------------------------------')
disp(' MATERIAL VALUES AT INFINTE FREQUENCY:');
disp('--------------------------------------')
disp(strcat('Density=  ',num2str(rho), 'Kg/m^3'));
disp(strcat('VP= ',num2str(vp_0), 'm/s'));
disp(strcat('VS= ',num2str(vs_0), 'm/s'));
% disp(strcat('mu= ',num2str(vs_0^2*rho), 'Pa'));
% disp(strcat('lambda= ',num2str(vp_0^2*rho-2*vs_0^2*rho), 'Pa'));
disp(' ');

disp('--------------------------------------')
disp(' ATTENUATION SETTINGS:');
disp('--------------------------------------')
disp(strcat('Number of mechanisms=  ',num2str(n)));
disp(strcat('Q for P-wave=  ',num2str(QPval)));
disp(strcat('Q for S-wave=  ',num2str(QSval)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PROBLEM SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derived quantities initialization
w_max=10001;          % Maximum frequency (minimum freq = -w_max). Is better to keep an odd number 
w_inc=0.1;           % Increment in frequency
w0=freq*(2*pi);        % Conversion to angular velocity
w=[0:w_inc:w_max];
w=2*w-w_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% GMB-EK COMPLEX M SOLUTION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equation system initialization %%
kmax=2*n-1;             %Number of equations to solve (the system is overdetermined)
AP=zeros(kmax,n);       %Initialization of system matrix (P wave)
AS=zeros(kmax,n);       %Initialization of system matrix (S wave)
QP=ones(kmax,1)/QPval;  %Desired values of Q for each mechanism inverted (P wave)
QS=ones(kmax,1)/QSval;  % " " (S wave)
YP=zeros(n,1);
YS=zeros(n,1);

%% Selection of the logarithmically equispaced frequencies
wmean=2*pi*freq;
wmin_disc=wmean/sqrt(f_ratio); 

if n==1
    w_disc=w0;
else
    for j=1:kmax
        w_disc(j)=exp(log(wmin_disc)+(j-1)/(kmax-1)*log(f_ratio));
    end
end

%% Filling of the linear system matrix %%
for m=1:kmax
    for j=1:n
        AP(m,j)=(w_disc(2*j-1).*w_disc(m)+w_disc(2*j-1).^2/QPval)./(w_disc(2*j-1).^2+w_disc(m).^2);
    end
end

for m=1:kmax
    for j=1:n
        AS(m,j)=(w_disc(2*j-1).*w_disc(m)+w_disc(2*j-1).^2/QSval)./(w_disc(2*j-1).^2+w_disc(m).^2);
    end
end

%% Solving of the system %%
YP=AP\QP;
YS=AS\QS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% VISUALIZATION OF Q IN FREQUENCY %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting values for the continuous representation of Q %%
wmin=wmean/sqrt(f_ratio); 
wmax=wmean*sqrt(f_ratio); 
xfrq=[wmin:(wmax-wmin)/(10000-1):wmax];

%% P-wave Q continuous values
numP=0;
denP=1;
for j=1:n
    numP=numP+(YP(j)*w_disc(2*j-1)*xfrq(:))./(w_disc(2*j-1)^2+xfrq(:).^2);
    denP=denP-(YP(j)*w_disc(2*j-1).^2)./(w_disc(2*j-1)^2+xfrq(:).^2);
end
Q_contP=denP./numP;

%% S-wave Q continuous values
numS=0;
denS=1;
for j=1:n
    numS=numS+(YS(j)*w_disc(2*j-1)*xfrq(:))./(w_disc(2*j-1)^2+xfrq(:).^2);
    denS=denS-(YS(j)*w_disc(2*j-1).^2)./(w_disc(2*j-1)^2+xfrq(:).^2);
end
Q_contS=denS./numS;

%% Computing fitting quality (RMS and maximum difference)
maxPdif=0;
maxSdif=0;

for j=1:length(Q_contP)
    tempP=abs(Q_contP(j)-QPval);
    if tempP >= maxPdif
        maxPdif=tempP;
    end
    tempS=abs(Q_contS(j)-QSval);
    if tempS >= maxSdif
        maxSdif=tempS;
    end
end

subplot(1,2,1),
semilogx(xfrq/(2*pi),Q_contP,xfrq/(2*pi),QPval*ones(length(xfrq),1)),
title('Q for the P-wave'),legend('computed','desired')
xlabel('frequency (Hz)'),ylabel('Q'),axis([wmin/(2*pi) wmax/(2*pi) 0 QPval*1.25])
subplot(1,2,2),
semilogx(xfrq/(2*pi),Q_contS,xfrq/(2*pi),QSval*ones(length(xfrq),1)),
title('Q for the S-wave'),legend('computed','desired'),
xlabel('frequency (Hz)'),ylabel('Q'),axis([wmin/(2*pi) wmax/(2*pi) 0 QSval*1.25])

disp(strcat('Attenuation bandwidth= [',num2str(wmin/(2*pi)),',',num2str(wmax/(2*pi)),'] Hz'));
disp('');
disp(strcat('Maximum QP fitting error=  ',num2str(maxPdif/QPval*100),' %'));
disp('');
disp(strcat('Maximum QS fitting error=  ',num2str(maxSdif/QSval*100),' %'));
disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% VELOCITIES COMPUTATION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complex modulus
MUP=vp_0^2*rho;    % Elastic P modulus
MUS=vs_0^2*rho;    % Elastic S modulus

MP=MUP*ones(length(w),1);
MS=MUS*ones(length(w),1);
vp=zeros(length(w),1);
vs=zeros(length(w),1);

for j=1:n
    MP(:)=MP(:)-MUP*(YP(j)*(w_disc(2*j-1)./(w_disc(2*j-1)+i*w(:))));
    MS(:)=MS(:)-MUS*(YS(j)*(w_disc(2*j-1)./(w_disc(2*j-1)+i*w(:))));
end

% Complex wave velocities (NOTE: Substitute by vp_0 and vs_0 to get the
% elastic solution)
vp(:)=sqrt(MP(:)/rho);
vs(:)=sqrt(MS(:)/rho);

%Computing which Munrelaxed values we need to get desired attenuation at
%frequency peak
PSI1P=1;
PSI2P=0;
for j=1:n
    PSI1P(:)=PSI1P(:)-YP(j)./(1+(w0/w_disc(2*j-1))^2);
    PSI2P(:)=PSI2P(:)+YP(j).*w0/w_disc(2*j-1)./(1+(w0/w_disc(2*j-1))^2);
end
R=sqrt(PSI1P.^2+PSI2P.^2);
PINF=rho*vp_t^2*(R+PSI1P)./(2*R.^2);

PSI1S=1;
PSI2S=0;
for j=1:n
    PSI1S(:)=PSI1S(:)-YS(j)./(1+(w0/w_disc(2*j-1))^2);
    PSI2S(:)=PSI2S(:)+YS(j).*w0/w_disc(2*j-1)./(1+(w0/w_disc(2*j-1))^2);
end
R=sqrt(PSI1S.^2+PSI2S.^2);
SINF=rho*vs_t^2*(R+PSI1S)./(2*R.^2);

LINF=PINF-2*SINF;
VPINF=sqrt(PINF/rho);
VSINF=sqrt(SINF/rho);

limw=w_max;
for j=1:length(w)
    if (w(j)-w0)^2<=limw
        wj=j;
        limw=(w(j)-w0)^2;
    end
end

disp('--------------------------------------')
disp(' RESULTING PHASE VELOCITIES AT PEAK FREQ:');
disp('--------------------------------------')
disp(strcat('Approx peak freq=  ',num2str(w(wj)/2/pi),' Hz'));
disp(strcat('P velocity at peak freq=  ',num2str(1./real(1./vp(wj))),' m/s'));
disp(strcat('S velocity at peak freq=  ',num2str(1./real(1./vs(wj))),' m/s'));
disp(' ');
disp('--------------------------------------')
disp(' RECOMMENDED UNRELAXED MODULI:');
disp('--------------------------------------')
disp(strcat('mu=  ',num2str(SINF),' Pa'));
disp(strcat('lambda=  ',num2str(LINF),' Pa'));
disp(strcat('VP associated=  ',num2str(VPINF),' m/s'));
disp(strcat('VS associated=  ',num2str(VSINF),' m/s'));

disp(' ');

figure,semilogx(w((end+1)/2:end)/2/pi,1./real(1./vp((end+1)/2:end))),ylabel('P-wave velocity (m/s)'),xlabel('Frequency (Hz)'),title('P-wave dispersion')
figure,semilogx(w((end+1)/2:end)/2/pi,1./real(1./vs((end+1)/2:end))),ylabel('S-wave velocity (m/s)'),xlabel('Frequency (Hz)'),title('S-wave dispersion')
