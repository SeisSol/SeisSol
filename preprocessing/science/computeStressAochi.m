function computeStressAochi()

strikeslip=false;
%sumatra hypocenter: strike 309 dip: -19
strike_rad = 309.*pi/180.
dip_rad= -19.*pi/180. 
%Material Parameter
rho=2670-1000;
g=9.8;
z=30e3;
cohesion = 0.4e6;
%intial level of stress (A5, Aochi and Madariaga 2003)
R=0.7;
%friction coefficients
mu_st=0.18;
mu_dy=0.12;


sigmazz=rho * g *z;

%most favorable direction (A4, AM2003)
Phi = pi/4.-0.5*atan(mu_st);

s2=sin(2.*Phi);
c2=cos(2.*Phi);

if (strikeslip==true)
    P=sigmazz;
    %ds (delta_sigma) is computed assuming that P=sigzz (A6, Aochi and Madariaga 2003)
    ds =  (mu_dy * P + R*(cohesion + (mu_st-mu_dy)*P)) / (s2 + mu_dy*c2 + R*(mu_st-mu_dy)*c2);
    sm=P;
    sii(1)= P + ds
    sii(2)= P
    sii(3)= P - ds
    
else % reverse fault

    c2bis = c2 - cos(2.*(Phi+dip_rad));
    %ds (delta_sigma) is deduced from R (A5, Aochi and Madariaga 2003), 
    %assuming that sig1 and sig3 are in the yz plane
    %sigzz and sigma_ini are then related by a phi+dip rotation (using A3, AM03)
    %sigmazz = sm  - ds * cos(2.*(Phi+dip_rad))
    %Now we have to assume that P = sm (not any more equal to sigmazz) 
    %and we can obtain the new expression of ds: 
    ds =  (mu_dy * sigmazz + R*(cohesion + (mu_st-mu_dy)*sigmazz)) / (s2 + mu_dy*c2bis + R*(mu_st-mu_dy)*c2bis);
    sm =  sigmazz + ds * cos(2.*(Phi+dip_rad));
    
    sii(1)= sm + ds
    %could be any value between sig1 and sig3
    sii(2)= sm
    sii(3)= sm - ds
end

%stress in the fault coordinate system
twoTheta_rad=2*Phi;
%different signs when comparing with the theory (tensor rotation)
% but it is because sigma1 on 2nd coordinate and sigma 3 on 1st coordinate
sigma_ini = sm - ds * cos(twoTheta_rad)
tau_ini = ds * sin(twoTheta_rad)
    
if (strikeslip==true)
    
    %stress in the cartesian coordinate system
    %pi/2- comes from the awkward convention on strike
    phi_xyz=-pi/2+ strike_rad + Phi;
    
%     twoTheta_rad = 2*phi_xyz;
%     sigmaxx_prime = sm - ds * cos(twoTheta_rad) 
%     sigmayy_prime = sm + ds * cos(twoTheta_rad) 
%     tauxy_prime = ds * sin(twoTheta_rad)

%     b11=sigmaxx_prime/P
%     b33=sigmayy_prime/P
%     b13=tauxy_prime/P
    
    c=cos(phi_xyz);
    s=sin(phi_xyz);
    R=[[c -s 0];[s c 0];[0 0 1]];
    Stress=[[sii(1) 0 0];[0 sii(3) 0];[0 0 sii(2)]];
    %minus sign to put in seissol convention
    Stress_cartesian = -R*Stress*(R')
    seissolInput = sprintf('%e,%e,%e,%e,%e,%e\n', Stress_cartesian(1,1),Stress_cartesian(2,2),Stress_cartesian(3,3),Stress_cartesian(1,2),Stress_cartesian(2,3),Stress_cartesian(1,3))
    %just for checking that everything is correct (if strike=0)
    Ratio_strike0 = (-Stress_cartesian(1,2)-mu_dy*Stress_cartesian(1,1))/(-cohesion + (mu_st-mu_dy)*Stress_cartesian(1,1))
    Stress_cartesian_norm=Stress_cartesian/sigmazz
    
else %faille inverse
    %stress in the cartesian coordinate system
    phi_xyz=dip_rad+Phi;
    twoTheta_rad = 2.*phi_xyz;
    %sigmaxx_prime = sa + sb * cos(twoTheta_rad) + tauxy * sin(twoTheta_rad)
    %sigmayy_prime = sa + sb * cos(twoTheta_rad) - tauxy * sin(twoTheta_rad)
    %tauxy_prime = -sb * sin(twoTheta_rad) + tauxy * cos(twoTheta_rad)

%     sigmayy_prime = sm + ds * cos(twoTheta_rad) 
%     sigmazz_prime = sm - ds * cos(twoTheta_rad) 
%     tauxz_prime = ds * sin(twoTheta_rad)
%     fprintf('sigmaxx_prime to be taken in range :%e %e\n',sii(1),sii(3))
%     %par exemple
%     sigmaxx_prime = sm
%     Stress=[[sigmaxx_prime 0 0];[0 sigmayy_prime tauxz_prime];[0 tauxz_prime sigmazz_prime]]
%     b11=sigmayy_prime/sigmazz
%     b33=sigmazz_prime/sigmazz
%     fprintf('b22 to be taken in range :%f %f\n',sii(1)/sigmazz,sii(3)/sigmazz)    
%     b13=tauxz_prime/sigmazz

    Stress=-[[sii(1) 0 0];[0 sii(2) 0];[0 0 sii(3)]]

    %first rotation: in xz plane
    phi_xyz=-(dip_rad+Phi);
    c=cos(phi_xyz);
    s=sin(phi_xyz);
    %R1=[[1 0 0];[0 c -s];[0 s c]]
    R1=[[c 0 s];[0 1 0];[-s 0 c]]
    
    %I cant explain the minus sign...
    c=cos(-strike_rad);
    s=sin(-strike_rad);
    %c=cos(pi/2-strike_rad);
    %s=sin(pi/2-strike_rad);  
    R2=[[c -s 0];[s c 0];[0 0 1]]
    %R1*Stress*(R1')

    %here it seems not right
    Stress_cartesian = R2*R1*Stress*(R1')*(R2')
    seissolInput = sprintf('%e,%e,%e,%e,%e,%e\n', Stress_cartesian(1,1),Stress_cartesian(2,2),Stress_cartesian(3,3),Stress_cartesian(1,2),Stress_cartesian(2,3),Stress_cartesian(1,3))
    Stress_cartesian_norm=-Stress_cartesian/sigmazz
end
    


%Centers of Mohr circles
C(1)=0.5*(sii(2)+sii(3));
C(2)=0.5*(sii(1)+sii(3));
C(3)=0.5*(sii(1)+sii(2));
%radius of Mohr circles
R(1)=abs(0.5*(sii(2)-sii(3)));
R(2)=abs(0.5*(sii(1)-sii(3)));
R(3)=abs(0.5*(sii(1)-sii(2)));


circle(C(2),0.,R(2));

% for i=1:3,
%    circle(C(i),0.,R(i));
% end

hold on

%current stress state
plot(sigma_ini,tau_ini,'+')

%static friction
x=linspace(min(sii),max(sii),50);
y=cohesion+mu_st*x;
plot(x,y)
%dynamic friction
y=mu_dy*x;
plot(x,y)
axis equal



function h = circle(x,y,r)

    hold on

    th = 0:pi/50:pi;

    xunit = r * cos(th) + x;

    yunit = r * sin(th) + y;

    h = plot(xunit, yunit);

    hold off
    return

