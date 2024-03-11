%Author Pouya Tavakoli
%Article: Generating synthetic ground motions reaching target spectrum with
%the optimization approach (2023)
function f=ObjFF(cc)
% in this file the obejctive funciton of the optimization is defined

load sandiegoAspectra.txt; % the file of the target acceleration spectra has to be loaded 
%in this file the target is the design acceleration spectrum of SanDiego
%city soil class A. the values of the target spectra was scaled to the g
%(gravity of earth). this file only contains the accelration values and the
%corresponding periods are added in the another seperate file named"earthquakeperiod"
targetspectra=sandiegoAspectra;
t=0:0.01:50; % duration of the earthquake is 50 seconds with 0.01s interval
NN=100;
TT=load('eathquakeperiod.txt'); % the period sampling according to the eq (22). 
%the earthquake will consist of these periods. also these periods are
%corresponding to the values of the acceleration response spectrum.
r=size(cc,1);f=NaN*ones(r,1);
for ii=1:r
for k=1:1:NN
%      u means acceleration on the gorund
u(k,:)=cc(ii,k)*cos(((2*pi)/(TT(k)))*t+cc(ii,k+100)); %equation (10)
end
uu=sum(u);%equation (10)

%% jeanings-honser-tai envlope for produce non stationary waves
% ENV =(t.^2/16).*(t<4)+1.*(4<t & t<35)+exp(-0.0357*(t-35)).*(35<t &
% t<80)+(0.05+0.0000938*(120-t).^2).*(t>80); %this envelope for generating
% the 80 seconds earthquake
ENV =(t.^2/16).*(t<4)+1.*(4<t & t<15)+exp(-0.0992*(t-15)).*(15<t & t<30)+(0.05+0.0005*(50-t).^2).*(t>30); %equation (12)
uu1=uu.*ENV; %equation(11)
uddg=uu1; %the generated earthquake
  
%% calculating the acceleration response spectrum of the generated earthquake (interpolation of excitation method)
DR=0.05; %5 percent damping ratio
dt=0.01; %interval time for the numerical solution
nn=length(uddg);
us(1)=0;
usdot(1)=0;
TT=TT';
T=TT; % at the vlaues of the sampling period, the values of displacement response spetrum will be calulated.
A=zeros(1,numel(T));V=zeros(1,numel(T));D=zeros(1,numel(T));
for i=1:length(T);
    %interpolation of excitation method 
    %readers can refer to Chopra structrual dynamics book.
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
     D(i)=max(abs(us)); %displacement response spectrum
     V(i)=ws.*D(i);     %velocity response spectrum
     A(i)=(ws^2).*D(i); %acceleration response spectrum
end
T=T';,D=D';,V=V';,A=A';

     f(ii)=norm(targetspectra(1:100)-A(1:100)); %equation (23)
         
end

