clear 
clc
close all

% User defined parameters 
A=0.68;            % Atmospheric emissivity coefficient (calibrated)
c_Uw=0.98;         % Cross wind scaling coefficient (calibrated)
c_MC=0.91;         % Microclimate coefficient (calibrated)
D_BG=5.8*10^(-6);  % Background diffusion coefficient (m2/s) (calibrated)
Ae=300*300;        % Circulator application areas (m2) per unit
Qac=0.2;           % Circulator flowrate (m3/s)

% Allocate imported array to column variable names
data= xlsread('/Inputdata.xlsx','Input','B2:K19753','',@convertSpreadsheetExcelDates);

Ta_Input = data(:,1);     % Reported air temperature (C)
Uwi_Input = data(:,2);    % Reported wind speed (m/s)
Jsn_Input = data(:,3);    % Reported solar radiation (W/m2)
RH_Input = data(:,4);     % Reported relative humidity
Depth_Input = data(:,5);  % Reported lake depth (m)
sd_Input = data(:,6);     % Field measured secchi depth (m)
CL_Input = data(:,7);     % Fraction of the sky covered with clouds
theta_Input = data(:,8);  % Reporated wind direction (degree)
ux_Input = data(:,9);     % Cross wind component of wind velocity (ux)
uy_Input = data(:,10);    % Longitudinal component of wind velocity (uy)  

% Interpolate hourly input data to 20 seconds time step data
dt=20;           % Time step (second)
Ti=25;           % Initial temperature, here, using mean air temperature of the first day
Multiplier=3600/20;  % time steps per hour  
newdata=zeros((size(data,1)-1)*Multiplier, size(data,2));

for i=1:size(data,1)-1
    for j=1:Multiplier       
    newdata(Multiplier*(i-1)+j,:)=data(i,:)+(j-1)*(data(i+1,:)-data(i,:))/Multiplier;
    end
end

% Interpolated 20-second input data
Ta = newdata(:,1);
Uwi = newdata(:,2);
Jsn = newdata(:,3);
RH = newdata(:,4);
Depth = newdata(:,5);
sd = newdata(:,6);
CL = newdata(:,7);
ux = newdata(:,9);
uy = newdata(:,10);

% Model setup
nt=size(newdata,1);   % Total number of time increments 
ns=15;             % Number of vertical segments of water column
Fa=0.4;            % Proportion of shortwave radiation absorbed in water surface                   
Abo=0.06;          % Albedo (unitless)
c=4186;            % Specific heat of water (J/Kg/C)
dz=Depth./ns;      % Length of segment (m)
c_ke=1.72;         % Parameter of extinction coefficient
B=0.1;             % Parameter in calculating diffusion (denominator)
n=-1;              % Parameter in calculating diffusion (power) 
qc=Qac/Ae;         % Artificial circulation rate (m/s)

% Adding cross wind scaling coefficient (c_Uw)
Uw=sqrt(c_Uw.*ux.^2+uy.^2);
 
% Matrix dimensions of ouput variables
Tw=zeros(nt,ns);      % Water temperature (C) 
ro=zeros(nt,ns);      % Density of fresh water (kg/m3) 
N2=zeros(nt,ns-1);    % Brunt-vaisala frequency
Ri=zeros(nt,ns-1);    % Richardson number
D=zeros(nt,ns-1);     % Turbulent diffusion (m2/s)
D0=zeros(nt,1);       % Diffusion at neutral stability (m2/s)
 
% Set initial condition
Tw(1,1:ns)=Ti*ones(1,ns); 

% Wind shear velocity (womega) (m/s); note explanation of equation below.
womega=0.00079.*Uw.^1.22;   
%{ 
 womega=sqrt(tao/ro_w), where tao is shear stress at the air-water surface (m/s), and ro_w is water density of 997 kg/m3 
 tao=(ro_a)*(Cd)*(Uw^2), where ro_a is air density of 1.2 kg/m3, Cd is drag coefficient (m/s)
 Cd=0.00052*(Uw^0.44), where drag coefficient is a function of wind speed Uw (m/s)
 Simplifying the above equations, we get womega=0.00079*Uw^1.22 
%}

% Other derived inputs
D0=Depth./34.*womega;     % Calculation of diffusion at neutral stability (D0) (m2/s)
ke=c_ke./sd;              % Calculation of extinction coefficient (ke) (/m)
ea=RH.*4.596.*exp(17.27.*Ta./(Ta+237.3));  % Calculation of air vapor pressure (ea)(mmHg)     
Jar=11.7.*10.^(-8).*(Ta+273).^4.*(1+0.17.*(CL).^2).*(A+0.031.*(ea).^0.5).*(1-0.03); % Calculation of air longwave radiation (Jar) (Ly/d)      
fUw=19+0.95.*Uw.^2; % Calculation of fUw, which represents the dependence of wind speed transfer on water surface  

% Calculating variables for every time step
for i=2:nt+1

ro(i-1,:)=999.842594+6.793952*10.^(-2)*(Tw(i-1,:))-9.09529*10.^(-3)*(Tw(i-1,:)).^2+1.001685*10.^(-4)*(Tw(i-1,:)).^3-...
    1.120083*10.^(-6)*(Tw(i-1,:)).^4+6.536332*10.^(-9)*(Tw(i-1,:)).^5;   % Calculating density of water (ro) (kg/m3) 
es(i-1)=4.596*exp(17.27*Tw(i-1,1)/(237.3+Tw(i-1,1)));  % Calculation of saturated vapor pressure (es) (mmHg) 
Jbr(i-1)=0.97*11.7*10^(-8)*(Tw(i-1,1)+273)^4;  % Calculation of waterbody longwave back radiation (Jbr) (Ly/d)
Je(i-1)=c_MC*fUw(i-1)*(es(i-1)-ea(i-1));   % Calculation of  evaporative flux (Je) (Ly/d)
Jc(i-1)=c_MC*0.47*fUw(i-1)*(Tw(i-1,1)-Ta(i-1));   % Calculation of  convective flux (Jc) (Ly/d)
Js(i-1)=(Fa*(1-Abo)*Jsn(i-1)+0.4846*(Jar(i-1)-Jbr(i-1)-Je(i-1)-Jc(i-1)));    % Calculation of total fluxes (Js) (W/m2), note 0.4846 is the unit converter from Ly/d to W/m2

% Calculating N2, Ri and D   
for p=1:ns-1
  N2(i-1,p)=(9.81/ro(i-1,p)*((ro(i-1,p+1)-ro(i-1,p)))/dz(i-1)); 
  Ri(i-1,p)=N2(i-1,p)/(womega(i-1)^2/(p*dz(i-1))^2);  
  D(i-1,p)=D_BG+D0(i-1)*(1+B*max((Ri(i-1,p)),0))^(n);  
end      

% Looping through both time and space and mix FULL water column
Tw(i,1)= Tw(i-1,1)+dt*((D(i-1,1)*(Tw(i-1,2)-Tw(i-1,1))/(dz(i-1))^2)+qc*0.5*((Tw(i-1,ns-1)-Tw(i-1,1))+...
    (Tw(i-1,ns)-Tw(i-1,2)))/2/dz(i-1)+((1-Fa)*Jsn(i-1)*(exp(-ke(i-1)*0)-exp(-ke(i-1)*(2-1)*dz(i-1))))/...
    ro(i-1,1)/c/(dz(i-1))+Js(i-1)/ro(i-1,1)/c/(dz(i-1)));     % Calculation of water temperature at surface     
Tw(i,2)=Tw(i-1,2)+dt*((D(i-1,1)*(Tw(i-1,1)-Tw(i-1,2))/(dz(i-1))^2+D(i-1,2)*(Tw(i-1,3)-Tw(i-1,2))/...
    (dz(i-1))^2)+qc*0.5*((Tw(i-1,ns)-Tw(i-1,2))+(Tw(i-1,1)-Tw(i-1,3)))/2/dz(i-1)+((1-Fa)*Jsn(i-1)*...
    (exp(-ke(i-1)*(1)*dz(i-1))-exp(-ke(i-1)*(2)*dz(i-1))))/ro(i-1,2)/c/(dz(i-1)));    % Calculation of water temperature at 2nd top segment 
Tw(i,ns)=Tw(i-1,ns)+dt*((D(i-1,ns-1)*(Tw(i-1,ns-1)-Tw(i-1,ns))/(dz(i-1))^2)+qc*0.5*((Tw(i-1,ns-2)-...
    Tw(i-1,ns))+(Tw(i-1,ns-1)-Tw(i-1,1)))/2/dz(i-1)+((1-Fa)*Jsn(i-1)*(exp(-ke(i-1)*(ns-1)*dz(i-1))))/...
    ro(i-1,ns)/c/(dz(i-1)));   % Calculation of water temperature at bottom
for q=3:ns-1
 Tw(i,q)=Tw(i-1,q)+dt*((D(i-1,q-1)*(Tw(i-1,q-1)-Tw(i-1,q))/(dz(i-1))^2+D(i-1,q)*(Tw(i-1,q+1)-Tw(i-1,q))/...
     (dz(i-1))^2)+qc*0.5*((Tw(i-1,q-2)-Tw(i-1,q))+(Tw(i-1,q-1)-Tw(i-1,q+1)))/2/dz(i-1)+((1-Fa)*Jsn(i-1)*...
     (exp(-ke(i-1)*(q-1)*dz(i-1))-exp(-ke(i-1)*(q)*dz(i-1))))/ro(i-1,q)/c/(dz(i-1)));     % Calculation of water temperature at middle
end

% Dealing with water column instability (at night)
delT=Tw(i,2:ns)-Tw(i,1:ns-1);
max_delT=max(delT);
if max_delT>0.05
 A1=find(delT==max_delT); % Returns linear indices corresponding to the entries of delT that equals max_delT
  for mix=A1:ns 
     if mix<ns
      Tw_top(i)=mean(Tw(i,1:mix));
       if Tw_top(i)>Tw(i,mix+1)
        Tw(i,1:mix)=Tw_top(i);
       break
       end
     end
     if mix==ns 
       if Tw(i,1)<mean(Tw(i,1:end))
        Tw(i,1:end)=mean(Tw(i,1:end));
       end
     end
   end
end
end %end loop through time

% Plot of hourly mean water temperature at various depths
Tw_output00= mean(reshape(Tw(2:end,1),Multiplier,[]))';   % Calculated water temperature at surface
Tw_output20= mean(reshape(Tw(2:end,3),Multiplier,[]))';   % Calculated water temperature at 20% of total depth
Tw_output40= mean(reshape(Tw(2:end,6),Multiplier,[]))';   % Calculated water temperature at 40% of total depth
Tw_output60= mean(reshape(Tw(2:end,9),Multiplier,[]))';   % Calculated water temperature at 60% of total depth
figure()
subplot(4,1,1);
plot(1:length(Tw_output00),Tw_output00)
title('Water temperature at surface')
subplot(4,1,2);
plot(1:length(Tw_output20),Tw_output20)
title('Water temperature at 20% depth')
subplot(4,1,3);
plot(1:length(Tw_output40),Tw_output40)
title('Water temperature at 40% depth')
subplot(4,1,4);
plot(1:length(Tw_output60),Tw_output60)
title('Water temperature at 60% depth')
xlabel('Hour')

% Plot of hourly mean diffusion at various depths
D_output00= mean(reshape(D(:,1),Multiplier,[]))';   % Calculated diffusion at surface
D_output20= mean(reshape(D(:,3),Multiplier,[]))';   % Calculated diffusion at 20% of total depth
D_output40= mean(reshape(D(:,6),Multiplier,[]))';   % Calculated diffusion at 40% of total depth
D_output60= mean(reshape(D(:,9),Multiplier,[]))';   % Calculated diffusion at 60% of total depth
figure()
subplot(4,1,1);
plot(1:length(D_output00),D_output00)
title('Diffusion at surface')
subplot(4,1,2);
plot(1:length(D_output20),D_output20)
title('Diffusion at 20% depth')
subplot(4,1,3);
plot(1:length(D_output40),D_output40)
title('Diffusion at 40% depth')
subplot(4,1,4);
plot(1:length(D_output60),D_output60)
title('Diffusion at 60% depth')
xlabel('Hour')

% Plot of hourly effective diffusion (m2/s) (i.e., turbulent diffusion + artificial circulation) at various depths
De_output00=D_output00+Qac.*14.*Depth_Input(2:end,1)./ns./Ae;   % Calculated effective diffusion at surface
De_output20=D_output20+Qac.*14.*Depth_Input(2:end,1)./ns./Ae;   % Calculated effective diffusion at 20% of total depth
De_output40=D_output40+Qac.*14.*Depth_Input(2:end,1)./ns./Ae;   % Calculated effective diffusion at 40% of total depth
De_output60=D_output60+Qac.*14.*Depth_Input(2:end,1)./ns./Ae;   % Calculated effective diffusion at 60% of total depth
figure()
subplot(4,1,1);
plot(1:length(De_output00),De_output00)
title('Effective diffusion at surface')
subplot(4,1,2);
plot(1:length(De_output20),De_output20)
title('Effective diffusion at 20% depth')
subplot(4,1,3);
plot(1:length(De_output40),De_output40)
title('Effective diffusion at 40% depth')
subplot(4,1,4);
plot(1:length(De_output60),De_output60)
title('Effective diffusion at 60% depth')
xlabel('Hour')
