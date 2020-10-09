% % % This script is written for SeisSol, to view the fault output for a list
% % % of dynamic rupture models in the same structure. Currently the parameter_viewlist={'ASl' 'Vr' 'stressdrop' 'magnitude' 'sliprate' 'PSR'};
% % %Change the directory and step_time to your model output
% % %@author: boli
clear all
close all
clc
model_dir='/import/freenas-m-05-seissol/bo/CHEESE/PSHA/HFF/simple_model/h5test/smoothed_close/';
model_name=dir([model_dir '*SH*']); %here SH is my model common name characters
post_dir=strcat(model_dir,'postanalysis');
output_dir_prefix='/output_o4/';
output_model_prefix='HFFtest';
step_time=0.25; %output timestep in second
mkdir(post_dir);

%%
tic
for ist=1:length(model_name)
    savedir=strcat(post_dir,'/',model_name(ist).name);
      mkdir(savedir);
namestr= strsplit(model_name(ist).name, '.');
   fcell=strcat(model_dir,model_name(ist).name,output_dir_prefix,output_model_prefix,'-fault_cell.h5'); %%all the model outputs have same structure

fvertex= strcat(model_dir,model_name(ist).name,output_dir_prefix,output_model_prefix,'-fault_vertex.h5');

vex = hdf5read(fvertex,'/mesh0/geometry');
x=vex(1,:);
y=vex(2,:);
z=vex(3,:);
ele = hdf5read(fcell,'/mesh0/connect');

ele_x=x(ele(:)+1);
ele_y=y(ele(:)+1);
ele_z=z(ele(:)+1);
loc_x=mean(reshape(ele_x,[size(ele,1),size(ele,2)]));
loc_y=mean(reshape(ele_y,[size(ele,1),size(ele,2)]));
loc_z=mean(reshape(ele_z,[size(ele,1),size(ele,2)]));

length_x=max(loc_x)-min(loc_x);
length_y=max(loc_y)-min(loc_y);
length_z=max(loc_z)-min(loc_z);


%%
ele=double(ele);
vex=vex';
ele = ele'+1;
tri = triangulation(ele,vex(:,1),vex(:,2),vex(:,3)); % prepare to plot trisurf

ASl=hdf5read(fcell,'/mesh0/ASl');
ndt=size(ASl,2); %the last time step
N_threshold=find(ASl(:,ndt)>=0.1); %you can change the threshold according to your model
N_below_threshold=find(ASl(:,ndt)<0.1); 

%%
display(['start plotting the output for model ',model_name(ist).name]);
figure5=figure('Position',[100,100,1500,1500],'visible','off')
subplot(3,2,1)
%%plot the moment rate figure and calculate magnitude
    mr_name=dir([model_dir model_name(ist).name output_dir_prefix output_model_prefix '-EnF*']); 
    M1=importdata([model_dir model_name(ist).name output_dir_prefix mr_name(1).name],' ',1);
    Mo_rate=zeros(size(M1.data,1),1);
    for isttn=1:length(mr_name);
        M_trans=importdata([model_dir model_name(ist).name output_dir_prefix mr_name(isttn).name],' ',1);
        Mo_rate=Mo_rate+M_trans.data(:,2);
        clear M_trans
    end
    dt=mean(diff(M1.data(:,1)));
    Mo=Mo_rate*dt;
    Mw=2/3*log10(sum(Mo))-6.07;
    plot(M1.data(:,1),Mo_rate,'k','LineWidth',2)
    set(gca,'Fontsize',10,'Fontweight','bold')
set(gca,'LineWidth',2)
xlabel('Tme (s)', 'fontsize', 10)
ylabel('Moment Rate', 'fontsize', 10)
title(strcat('Mw=',num2str(Mw)), 'fontsize', 10)
% print(gcf,'-djpeg',[savedir,'/MomentRate_plot']) %%you can change if you
% want each parameter output in one figure
clear Mo dt Mw M1 mr_name
%%
%plot the slip without threshold filter, fault view
subplot(3,2,2)
trisurf(tri,ASl(:,ndt),'edgecolor','none');
hold on
pbaspect([length_x length_y length_z]/min([length_x length_y length_z]))
colorbar
colorbar_limit_refer=ASl(N_threshold,ndt);
caxis([0 prctile(colorbar_limit_refer(:),95)])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'LineWidth',2)
xlabel('EW (m)', 'fontsize', 10)
ylabel('NS (m)', 'fontsize', 10)
zlabel('Depth (m)', 'fontsize', 10)
title('Total Slip (m)', 'fontsize', 10)

%%
%plot the rupture velocity with threshold filter, fault view
Vr=hdf5read(fcell,'/mesh0/Vr');
subplot(3,2,3)
Vr1=Vr;
Vr1(N_below_threshold,ndt)=NaN;
trisurf(tri,Vr1(:,ndt),'edgecolor','none');

hold on
pbaspect([length_x length_y length_z]/min([length_x length_y length_z]))
colorbar
caxis([0 prctile(Vr1(:),99)])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'LineWidth',2)
xlabel('EW (m)', 'fontsize', 10)
ylabel('NS (m)', 'fontsize', 10)
zlabel('Depth (m)', 'fontsize', 10)
title('Rupture Velocity (m/s)', 'fontsize', 10)

%plot the rupture velocity with threshold filter, histogram
subplot(3,2,4)
Vr_threshold=Vr(N_threshold,ndt);
bins=200:200:8000;
hist(Vr_threshold,bins)
xlim([0 8000])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'LineWidth',2)
xlabel('Rupture Velocity (m/s)', 'fontsize', 10)
ylabel('Element counts', 'fontsize', 10)
title(strcat('mean Vr=',num2str(mean(Vr_threshold))), 'fontsize', 10)
clear Vr

%%
%%stress drop
subplot(3,2,5)
T_s=hdf5read(fcell,'/mesh0/T_s');
T_d=hdf5read(fcell,'/mesh0/T_d');
stress_drop=sqrt(T_s(:,ndt).^2+T_d(:,ndt).^2)/(1e6);
stress_drop_filter=stress_drop;
stress_drop_filter(N_below_threshold)=NaN;
trisurf(tri,stress_drop_filter,'edgecolor','none');

hold on
pbaspect([length_x length_y length_z]/min([length_x length_y length_z]))
colorbar
caxis([0 prctile(stress_drop_filter(:),90)])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'LineWidth',2)
xlabel('EW (m)', 'fontsize', 10)
ylabel('NS (m)', 'fontsize', 10)
zlabel('Depth (m)', 'fontsize', 10)
title(['Stress Drop mean=',num2str(mean(stress_drop(N_threshold))),'Mpa'], 'fontsize', 10)
clear T_s T_d stress_drop stress_drop_filter

%%
%plot the peak slip rate with threshold filter, fault view
PSR=hdf5read(fcell,'/mesh0/PSR');
subplot(3,2,6)
PSR1=PSR;
PSR1(N_below_threshold,ndt)=NaN;
trisurf(tri,PSR1(:,ndt),'edgecolor','none');

hold on
pbaspect([length_x length_y length_z]/min([length_x length_y length_z]))
colorbar
caxis([0 prctile(PSR1(:),98)])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'LineWidth',2)
xlabel('EW (m)', 'fontsize', 10)
ylabel('NS (m)', 'fontsize', 10)
zlabel('Depth (m)', 'fontsize', 10)
title('Peak Slip Rate (m/s)', 'fontsize', 10)

sgtitle(['Model2closed ',model_name(ist).name],'Interpreter','none');
print(gcf,'-djpeg',[post_dir,'/Model2closed_',model_name(ist).name,'-postview'])

%if the model_name contains ".", then you can use the following line to
%save figures
% print(gcf,'-djpeg',[post_dir,'/Model2closed_',namestr{1},'_',namestr{2},'-postview'])
%%
%%slip rate
SRd=hdf5read(fcell,'/mesh0/SRd');
SRs=hdf5read(fcell,'/mesh0/SRs');
SR_window=round(linspace(1,ndt,10)); %select the middle 2-7 time step to show slip rate
SR_plot=sqrt(SRd(:,SR_window(3)).^2+SRs(:,SR_window(3)).^2); %use the 3 window to help set the color bar limit
SR_plot(SR_plot<0.0001)=NaN;
SR_colorlimit=prctile(SR_plot,98);

display(['start ploting the slip rate for model ',model_name(ist).name]);
figure10=figure('Position',[100,100,1000,1500])
for isr=1:6
    subplot(3,2,isr)
SR_plot=sqrt(SRd(:,SR_window(isr+1)).^2+SRs(:,SR_window(isr+1)).^2);
SR_filter=SR_plot;
SR_filter(N_below_threshold)=NaN;
trisurf(tri,SR_filter,'edgecolor','none');
title(['Time ',num2str(SR_window(isr+1)*step_time-step_time),'s'], 'fontsize', 10)

pbaspect([length_x length_y length_z]/min([length_x length_y length_z]))
colorbar
caxis([0 SR_colorlimit])
axis off
clear SR_plot SR_filter
end
sgtitle(['Model2 closed ',model_name(ist).name,' slip rate'],'Interpreter','none');
print(gcf,'-djpeg',[post_dir,'/Model2closed_',model_name(ist).name,'-sliprate'])
%print(gcf,'-djpeg',[post_dir,'/Model2closed_',namestr{1},'_',namestr{2},'-sliprate'])
clear SR_colorlimit SRs SRd
%%
clear ASl loc_* ele_* x y z N_threshold
 close all
end
toc

