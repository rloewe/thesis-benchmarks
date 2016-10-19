%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_all=Stommel2box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a Stommel-like 2box model, for APM 115 workshop
% Eli Tziperman
% Modified by Laure Zanna, for EPS 131
close all; clear all;

p=set_parameters;
ode_options=odeset('RelTol',1.e-7,'abstol',1.e-8);

% find steady state solutions
tspan=[0:50:1500]*p.year;
for ii=1:3
  if ii==1
    X0=[36 35]'; p.Temperature=[25 5]; 
  elseif ii==2
    X0=[36 34.255]; p.Temperature=[25 5];
  else 
    X0=[36 34.5]'; % (S1 S2) - Initial conditions
    p.Temperature=[10 5]; % (T1 T2)
  end
[time,X] = ode45(@(t,X) rhs(t,X,p),tspan,X0,ode_options);

% plots:
figure
subplot(2,1,1)
plot(time/p.year,X(:,1),'r','linewidth',2)
hold on
plot(time/p.year,X(:,2),'b','linewidth',2)
ylabel('S','fontsize',14);
legend('S box 1','S box 2')
if ii==1
  title(['Thermally driven solution'],'fontsize',14,'fontweight','b')
elseif ii==2
  title(['Thermally driven solution starting near the unstable state'],'fontsize',14,'fontweight','b')
elseif ii==3 
  title(['Salinity driven solution'],'fontsize',14,'fontweight','b')
end
  % Initial DeltaS
X(1,2)-X(1,1)

% Final DeltaS
X(end,2)-X(end,1)

subplot(2,1,2)
i=0;
for t=tspan
  i=i+1;
  q12plot(i)=q12(X(i,:),p)/p.Sverdrup;
end
plot(time/p.year,q12plot,'r','linewidth',2)
hold on
xlabel('time (yr)','fontsize',14);
ylabel('Transport [Sv]','fontsize',14);

FW_tmp=p.S0*p.FW_amp0*p.meter/p.year;
% Thermally Driven Flow 
SteadyState1=(0.5*p.alpha*(p.Temperature(2)-p.Temperature(1))+...
    0.5*sqrt((p.alpha*(p.Temperature(2)-p.Temperature(1)))^2-(4*p.beta/p.k)*FW_tmp*p.area))/p.beta
%figure(1); plot(FW_amp0,SteadyState1/p.beta,'b+'); hold on;
SteadyState2=(0.5*p.alpha*(p.Temperature(2)-p.Temperature(1))-...
    0.5*sqrt((p.alpha*(p.Temperature(2)-p.Temperature(1)))^2-(4*p.beta/p.k)*FW_tmp*p.area))/p.beta
%plot(FW_amp0,SteadyState2/p.beta,'ro'); hold on;
% Salt driven flow
SteadyState3=(0.5*p.alpha*(p.Temperature(2)-p.Temperature(1))+...
    0.5*sqrt((p.alpha*(p.Temperature(2)-p.Temperature(1)))^2+(4*p.beta/p.k)*FW_tmp*p.area))/p.beta
%plot(FW_amp0,SteadyState3/p.beta,'k*'); hold on;
SteadyState4=(0.5*p.alpha*(p.Temperature(2)-p.Temperature(1))-...
    0.5*sqrt((p.alpha*(p.Temperature(2)-p.Temperature(1)))^2+(4*p.beta/p.k)*FW_tmp*p.area))/p.beta
%plot(FW_amp0,SteadyState4/p.beta,'g.'); hold on;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flow=q12(X,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DeltaT=p.Temperature(2)-p.Temperature(1);
DeltaS=X(2)-X(1);
flow=-p.k*(p.alpha*DeltaT-p.beta*DeltaS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=rhs(time,X,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaS12=X(2)-X(1);
Fs=FW(time,p);

r=[Fs(1)*p.area+abs(q12(X,p))*deltaS12, ...
   Fs(2)*p.area-abs(q12(X,p))*deltaS12]'/p.volume;

if p.optimization==1
  r=r*1.e12;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fs=FW(time,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface salt flux in northern box

ramp=1;
amp=p.FW_amp0;

Fs(1)=   ramp*p.S0*amp*p.meter/p.year;
Fs(2)=  -ramp*p.S0*amp*p.meter/p.year;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=set_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.meter=1;
p.sec=1;
p.km=1000*p.meter;
p.area=(4000*p.km)^2;
p.depth=4*p.km;
p.volume=p.area*p.depth;
p.year=365*24*3600;
p.time_max=500000*p.year;
p.alpha=4*10^-4; % K^-1
p.beta=4*10^-3; % ppt^-1
p.Sverdrup=10^6*p.meter^3/p.sec;
p.k=10*p.Sverdrup/0.005; % Sv/(delta sigma)
p.S0=35.0; % Salinity reference
p.FW_amp0=0.2; % meter/year for calculating steady state that
               % initializes the hysteresis.
p.Nrand=50;
p.optimization=0;
