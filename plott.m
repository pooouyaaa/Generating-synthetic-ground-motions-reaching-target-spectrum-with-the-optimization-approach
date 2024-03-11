%Author Pouya Tavakoli
%Article: Generating synthetic ground motions reaching target spectrum with
%the optimization approach (2023)
clc;clear;close all
%in this file the detail information of the generated earthquake is
%extracted
load sandiegoAspectra.txt; % loading the target spectra.
targetspectra=sandiegoAspectra;
TT=load('eathquakeperiod.txt'); %loading the period sampling
T=TT';
load xmin.txt;%the optimum values for the optimization variables
xmin=xmin';
NN=100;                 % number of the Cosine Function which sum togheter for generating articifical earthquake
t=0:0.01:50;            % duration of the earthquake
Td=max(t);
DR=0.05;                % the viscous damping value for calculating acceleration spectra
dt=0.01;                % time intervals for calculation respanse spectrum with Duhamel's integral       
fq=1./T;
%------------------------------------------------------------------------------------------------------
for k=1:1:NN
%      u means acceleration on the gorund
u(k,:)=xmin(1,k)*cos(((2*pi)/(TT(k)))*t+xmin(1,k+100));
% if k<12
%     figure(1) 
%     plot(t,u)
%     xlim([0 4])
% %     ylim([-0.035 0.035])
%     xlabel('Time(s)') 
%     ylabel('acceleration(g)') 
%     title('value of each frequency between 0.15HZ to 1.65HZ') 
%     grid on
%     grid on
%     drawnow
%     pause(0.3)
% end
%figure one shows the amplitude and phase angles of the first 12 Cosine
%Functions
end
uu=sum(u);
% ENV =(t.^2/16).*(t<4)+1.*(4<t & t<35)+exp(-0.0357*(t-35)).*(35<t & t<80)+(0.05+0.0000938*(120-t).^2).*(t>80);
% envelope type B
ENV =(t.^2/16).*(t<4)+1.*(4<t & t<15)+exp(-0.0992*(t-15)).*(15<t & t<30)+(0.05+0.0005*(50-t).^2).*(t>30);
uddg=uu.*ENV;   % generation non-stationary earthquake

% figure(2) % show the type B of Envelope function 
% plot(t,ENV) 
% xlabel('Time(s)') 
% ylabel('ENVELOPE FUNCTION') 
% title('ENVELOPE') 
% grid on

% figure(3)
% plot(t,uu)
% xlabel('Time(s)') 
% ylabel('acceleration(g)') 
% title('earthquake record with stationary waves') 
% grid on
%----------------------------------------------------------------------------------
 % modification of base line and velocity

A1=trapz(t,uddg);
A2=trapz(t,uddg.*t);
A3=trapz(t,uddg.*(t.^2));

KK=[1 Td/2 (Td^2)/3 ; 0.5 Td/3 (Td^2)/4 ; 1/3 Td/4 (Td^2)/5 ];
FF=[A1./Td; A2/(Td^2); A3/(Td^3)];
CC=inv(KK)*FF;    
EE=CC(1)+CC(2)*t+CC(3)*(t.*t);
 uddg_corrected=uddg-EE;
 velocity_corrected=cumtrapz(t,uddg_corrected*9.81);  %product9.81 because acceleration calculated based on g
 displacement_coreected=cumtrapz(t,velocity_corrected);
 
  figure(4)
    plot(t,uddg_corrected) %the acceleration of corrected earthquake is plotted
    xlabel('Time (s)')
ylabel('Acceleration (g)')
maxuddg=max(abs(uddg_corrected));
 title(['Maximum acceleration= ', num2str(maxuddg),'g'])
    grid on
    set(gca,'XTick',[], 'YTick', [])
 figure(5)
    plot(t,velocity_corrected) %the velocity of corrected earthquake is plotted
    xlabel('Time (s)')
ylabel('Velocity (m/s)')
maxvelocity=max(abs(velocity_corrected));
 title(['Maximum velocity= ', num2str(maxvelocity),'m/s'])
 grid on
    
    figure(6)
    plot(t,displacement_coreected) %the displacement of corrected earthquake is plotted
    xlabel('Time (s)')
ylabel('Displacement (m)')
maxdisp=max(abs(displacement_coreected));
 title(['Maximum displacement= ', num2str(maxdisp),'m'])
    grid on
    uddg_corrected=uddg_corrected';
    save uddgfinal.txt uddg_corrected -ascii  % the accelration of earthquake is saved
    %-----------------------------------------------------
    PlotPos = [1 1 20 20];
hFig1 = figure(30); %plotting displacement, velocity and acceleration of earthquake
set(hFig1,'units','centimeters','position',PlotPos)
hold on
   subplot(3,1,1)
     plot(t,uddg_corrected)
    xlabel('Time (s)','FontSize',14)
ylabel('Acceleration (g)','FontSize',14)
maxuddg=max(abs(uddg_corrected));
 title(['Maximum acceleration= ', num2str(maxuddg,'%.2f'),'g'],'FontSize',11)
 subplot(3,1,2)
  plot(t,velocity_corrected)
    xlabel('Time (s)','FontSize',14)
ylabel('Velocity (m/s)','FontSize',14)
maxvelocity=max(abs(velocity_corrected));
 title(['Maximum velocity= ', num2str(maxvelocity,'%.2f'),'m/s'],'FontSize',11)
 subplot(3,1,3)
    plot(t,displacement_coreected)
    xlabel('Time (s)','FontSize',14)
ylabel('Displacement (m)','FontSize',14)
maxdisp=max(abs(displacement_coreected));
 title(['Maximum displacement= ', num2str(maxdisp,'%.2f'),'m'],'FontSize',11)
   print(gcf,'1-article.png','-dpng','-r500');  
% -----------------------------------------------------------------------------------------------------------
nn=length(uddg);
us(1)=0;
usdot(1)=0;
A=zeros(1,numel(T));V=zeros(1,numel(T));D=zeros(1,numel(T));
for i=1:length(T);
    ws=(2*pi)/T(i);
    wD=ws*sqrt(1-(DR)^2);
    A1=exp(-DR*ws*dt)*(DR*sin(wD*dt)/(1-DR^2)^0.5+cos(wD*dt));
    B=exp(-DR*ws*dt)*(sin(wD*dt)/wD);
    C=-(1/ws)^2*(2*DR/ws/dt+exp(-DR*ws*dt)*(((1-2*DR^2)/wD/dt-DR/(1-DR^2)^0.5)*sin(wD*dt)-(1+2*DR/ws/dt)*cos(wD*dt)));
    D1=-(1/ws)^2*(1-2*DR/ws/dt+exp(-DR*ws*dt)*(((-1+2*DR^2)/wD/dt)*sin(wD*dt)+(2*DR/ws/dt)*cos(wD*dt)));
    a=-exp(-DR*ws*dt)*(ws*sin(wD*dt)/(1-DR^2)^0.5);
    b=exp(-DR*ws*dt)*(-DR*sin(wD*dt)/(1-DR^2)^0.5+cos(wD*dt));
    c=-(1/ws)^2*(-1/dt+exp(-DR*ws*dt)*((ws/(1-DR^2)^0.5+DR/(1-DR^2)^0.5/dt)*sin(wD*dt)+cos(wD*dt)/dt));
    d=-(1/ws)^2/dt*(1-exp(-DR*ws*dt)*(DR*sin(wD*dt)/(1-DR^2)^0.5+cos(wD*dt)));
    for j=1:nn-1;
        us(j+1)=A1*us(j)+B*usdot(j)+C*uddg(j)+D1*uddg(j+1);
        usdot(j+1)=a*us(j)+b*usdot(j)+c*uddg(j)+d*uddg(j+1);
    end
     D(i)=max(abs(us));
     V(i)=ws.*D(i);
     A(i)=(ws^2).*D(i);
end
T=T';,D=D';,V=V';,A=A';
%--------------------------------------------------------------------------------------------------------
figure(8)
plot(T,A,'r','Linewidth',2)
hold on
plot(T,targetspectra,'b','Linewidth',1)
xlabel('Period (s)')
xlim([0 8.5])
ylabel('Acceleration specrtra (g)') 
title('compare target spectra vs artificial earthquake')
legend('Artificial earthquake','Target DRS')
grid on

figure(9)
plot(T,V*9.81,'b','Linewidth',1.5)   %product to 9.81 because acceleration was calculated based on g
xlim([0 8.5])
xlabel('Period (s)') 
ylabel('Velocity specrtra (m/s)') 
title('velocity spectra of Artificial Earthquake')
grid on

figure(10)
plot(T,D*9.81,'b','Linewidth',1.5) %product to 9.81 because acceleration was calculated based on g
xlim([0 8.5])
xlabel('Period (s)') 
ylabel('Displacement specrtra (m)') 
title('Displacement spectra of Artificial Earthquake')
grid on
%--------------------------------
   PlotPos = [1 1 20 20];
hFig1 = figure(40); 
set(hFig1,'units','centimeters','position',PlotPos)
hold on
subplot(4,1,[1 2])
plot(T,A,'r','Linewidth',2)
hold on
plot(T,targetspectra,'b','Linewidth',1)
xlabel('Period (s)','FontSize',14)
xlim([0 8.5])
ylabel('Acceleration specrtra (g)','FontSize',14) 
 title('compare target spectra vs artificial earthquake')
legend('Artificial earthquake','Target DRS','FontSize',12)
axes('position',[0.3, 0.7, 0.3, 0.15]);
plot(T,A,'r',T,targetspectra,'b','Linewidth',2)
xlim([1.2 2.5])
subplot(4,1,3)
plot(T,V*9.81,'b','Linewidth',1.5)   %product9.81 because acceleration calculated based on g
xlim([0 8.5])
xlabel('Period (s)','FontSize',14) 
ylabel('S_{V} (m/s)','FontSize',14)
ylim([0,0.5])
subplot(4,1,4)
plot(T,D*9.81,'b','Linewidth',1.5) %product9.81 because acceleration calculated based
xlim([0 8.5])
xlabel('Period (s)','FontSize',14) 
ylabel('S_{D} (m)','FontSize',14)
print(gcf,'2-article.png','-dpng','-r500');  
%---------------------------------------


figure(11)
na=length(uddg_corrected) ;
n=nextpow2(na) ;
ffta=fft(uddg_corrected,2^n);
ampffta=abs(ffta(1:length(ffta)/2))*0.01;
f=1/(2*0.01)*linspace(0,1,length(ampffta)) ;
plot(f,ampffta) ;
 xlim([0 6])
xlabel('Frequency(Hz)');
ylabel('Fourier amplitude');
 title('Frequency content of Artificial Earhquake ');
 grid on;


%% Spectogram of Artificial Earthquake(Sandiego-site:A)
tslide=0:0.1:max(t);
spect=[];
for j=1:length(tslide)  
    %wavelet func
w=0.5*exp(-0.05*pi*(t-tslide(j)).^2).*sin(pi*(t-tslide(j)));
%  figure(12)
%   subplot(3,1,1), plot(t,uddg_corrected,'b',t,w,'r')
%   xlabel('time(s)');
%  ylabel('acceleration');
%  title('earthquake record ');
%  grid on;
agw=w.*uddg_corrected';
%   subplot(3,1,2), plot(t,agw)
%   xlabel('time(s)');
%  ylabel('fraction of signal');
%  title('wavelet ');
%  grid on;
 fftwaveltet=fft(agw,2^n);
 ampfftwavelet=abs(fftwaveltet(1:length(fftwaveltet)/2))*0.01;
 f=1/(2*0.01)*linspace(0,1,length(ampfftwavelet));
 spect(j,:)=ampfftwavelet ;
%  subplot(3,1,3), plot(f,ampfftwavelet/max(ampfftwavelet)) ;
%  xlabel('Frequency(Hz)');
%  ylabel('fourier amplitude');
%  title('FFT');
%  xlim([0 16]);
%  grid on;
%  drawnow
%  pause(0.000001)
end
figure(13)
a=spect';
pcolor(tslide,f,a), shading interp;
set(gca);
colormap(hsv);
xlabel('Time (S)');
 ylabel('Frequency (HZ)');
 ylim([0 6]);
  title('Spectogram of Artificial Earthquake(Sandiego-site:A)');
% ----------------------------------------------------------------------
 PlotPos = [5 5 38 12];
hFig1 = figure(50); 
set(hFig1,'units','centimeters','position',PlotPos)
hold on
subplot(1,3,1)
plot(f,ampffta) ;
 xlim([0 6])
xlabel('Frequency(Hz)','FontSize',14);
ylabel('Fourier amplitude','FontSize',14);
subplot(1,3,[2 3])
pcolor(tslide,f,a), shading interp;
set(gca);
colormap(hsv);
colorbar
xlabel('Time (s)','FontSize',14);
 ylabel('Frequency (Hz)','FontSize',14);
 ylim([0 6]);
  print(gcf,'3-article.png','-dpng','-r150');  
 %% relative error (using equation 30 in the article)
 delta= A-targetspectra;
 norm_delta=norm (delta);
 norm_target=norm (targetspectra);
 relative_error= (norm_delta/norm_target)*100








